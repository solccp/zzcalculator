/*
 * This software is provided under the terms of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * Copyright (c) 2007 Tom Portegys, All Rights Reserved.
 * Permission to use, copy, modify, and distribute this software
 * and its documentation for NON-COMMERCIAL purposes and without
 * fee is hereby granted provided that this copyright notice
 * appears in all copies.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.
 */

/*
 * Compute an MD5 hash for a graph.
 */

#include "graphHash.hpp"

// Constructor.
GraphHash::GraphHash(Graph *graph, bool hashLabels)
{
   int         i, i2;
   VertexCoder *child;

   map<pair<Graph::Vertex *, Graph::Edge *>,
       VertexCoder *, VertexCoder::ltcmpConnection> *vertexMap;

   // Save graph.
   this->graph = graph;

   // Sort vertices by label.
   for (i = 0, i2 = (int)graph->vertices.size(); i < i2; i++)
   {
      vertices.push_back(graph->vertices[i]);
   }
   sort(vertices.begin(), vertices.end(), ltcmpVertices);

   // Generate code for graph using vertex coders.
   vertexCoder = new VertexCoder();
   assert(vertexCoder != NULL);
   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      vertexMap = new map<pair<Graph::Vertex *, Graph::Edge *>,
                          GraphHash::VertexCoder *, VertexCoder::ltcmpConnection>;
      assert(vertexMap != NULL);
      child = new VertexCoder(vertices[i], NULL, vertexCoder, 1, vertexMap);
      assert(child != NULL);
      vertexCoder->children.push_back(child);
   }
   vertexCoder->generateCode(hashLabels);
   memcpy(code, vertexCoder->code, MD5_SIZE);
}


// Destructor.
GraphHash::~GraphHash()
{
   delete vertexCoder;
   vertices.clear();
}




// Less-than comparison of vertices by label.
bool GraphHash::ltcmpVertices(Graph::Vertex *a, Graph::Vertex *b)
{
   return(a->label < b->label);
}


// Vertex coder constructors.
GraphHash::VertexCoder::VertexCoder()
{
   vertex     = NULL;
   parentEdge = NULL;
   codeValid  = false;
   creator    = NULL;
   expanded   = false;
   generation = 0;
   vertexMap  = NULL;
}


GraphHash::VertexCoder::VertexCoder(Graph::Vertex *vertex,
                                    Graph::Edge *parentEdge, VertexCoder *creator, int generation,
                                    map<pair<Graph::Vertex *, Graph::Edge *>,
                                        VertexCoder *, ltcmpConnection> *vertexMap)
{
   this->vertex     = vertex;
   this->parentEdge = parentEdge;
   codeValid        = false;
   this->creator    = creator;
   expanded         = false;
   this->generation = generation;
   this->vertexMap  = vertexMap;
}


// Vertex coder destructor.
GraphHash::VertexCoder::~VertexCoder()
{
   int i, i2;

   // Purge children created by other coders to avoid duplicate deletions.
   if (vertex == NULL)
   {
      for (i = 0, i2 = (int)children.size(); i < i2; i++)
      {
         if (children[i]->vertexMap != NULL)
         {
            delete children[i]->vertexMap;
         }
         children[i]->purgeChildren();
      }
   }
   for (i = 0, i2 = (int)children.size(); i < i2; i++)
   {
      if (children[i] != NULL)
      {
         delete children[i];
      }
   }
   children.clear();
}


// Purge children created by other coders.
void GraphHash::VertexCoder::purgeChildren()
{
   int i, i2;

   for (i = 0, i2 = (int)children.size(); i < i2; i++)
   {
      if (children[i]->creator == this)
      {
         children[i]->purgeChildren();
      }
      else
      {
         children[i] = NULL;
      }
   }
}


// Generate code.
void GraphHash::VertexCoder::generateCode(bool hashLabels)
{
   int               i, j, s;
   struct MD5Context md5c;
   unsigned char     *input;

   // Root?
   s = (int)children.size();
   if (vertex == NULL)
   {
      // Iteratively expand children until fully expanded.
      for (i = 0; i < s; i++)
      {
         children[i]->expand();
      }
      while (true)
      {
         // Generate codes without labels to encode structure.
         for (i = 0; i < s; i++)
         {
            children[i]->generateCode(false);
         }

         // All vertices are fully expanded?
         for (i = 0; i < s; i++)
         {
            if (!children[i]->expanded)
            {
               break;
            }
         }
         if (i == s)
         {
            break;
         }

         // Expand vertices.
         for (i = 0; i < s; i++)
         {
            if (!children[i]->expanded)
            {
               children[i]->expand();
            }
         }
      }
   }

   // Recursively generate code.
   unsigned char oldCode[MD5_SIZE];
   memcpy(oldCode, code, MD5_SIZE);
   for (i = 0; i < s; i++)
   {
      children[i]->generateCode(hashLabels);
   }
   sort(children.begin(), children.end(), ltcmpCode);
   MD5Init(&md5c);
   j = 0;
   if (vertex != NULL)
   {
      if (hashLabels)
      {
         j += sizeof(vertex->label);
      }
      if (parentEdge != NULL)
      {
         if (hashLabels)
         {
            j += sizeof(parentEdge->label);
         }
         j++;
      }
   }
   input = new unsigned char[(s * MD5_SIZE) + j];
   assert(input != NULL);
   if (vertex != NULL)
   {
      if (hashLabels)
      {
         memcpy(input, &vertex->label, sizeof(vertex->label));
      }
      if (parentEdge != NULL)
      {
         if (hashLabels)
         {
            memcpy(&input[sizeof(vertex->label)], &parentEdge->label, sizeof(parentEdge->label));
         }
         if (parentEdge->directed)
         {
            if (parentEdge->source == vertex)
            {
               input[j - 1] = 0;
            }
            else
            {
               input[j - 1] = 1;
            }
         }
         else
         {
            input[j - 1] = 2;
         }
      }
   }
   for (i = 0; i < s; i++)
   {
      memcpy(&input[j], children[i]->code, MD5_SIZE);
      j += MD5_SIZE;
   }
   MD5Update(&md5c, input, j);
   MD5Final(code, &md5c);
   delete [] input;

   // Mark vertex as expanded if code not changed.
   if (codeValid)
   {
      if (memcmp(oldCode, code, MD5_SIZE) == 0)
      {
         expanded = true;
      }
   }
   codeValid = true;
}


// Expand coder one level deeper.
void GraphHash::VertexCoder::expand()
{
   int         i, i2;
   VertexCoder *child;

   pair<Graph::Vertex *, Graph::Edge *> key;
   map<pair<Graph::Vertex *, Graph::Edge *>,
       VertexCoder *, ltcmpConnection>::iterator itr;

   if (expanded)
   {
      return;
   }
   if (children.empty())
   {
      // Expand this vertex.
      for (i = 0, i2 = (int)vertex->edges.size(); i < i2; i++)
      {
         if (vertex == vertex->edges[i]->source)
         {
            key.first  = vertex->edges[i]->target;
            key.second = vertex->edges[i];
            if ((itr = vertexMap->find(key)) == vertexMap->end())
            {
               child = new VertexCoder(vertex->edges[i]->target, vertex->edges[i], this, generation + 1, vertexMap);
               assert(child != NULL);
               children.push_back(child);
               (*vertexMap)[key] = child;
            }
            else
            {
               child = itr->second;

               // Share next generation.
               if (child->generation == generation + 1)
               {
                  children.push_back(child);
               }
            }
         }
         else
         {
            key.first  = vertex->edges[i]->source;
            key.second = vertex->edges[i];
            if ((itr = vertexMap->find(key)) == vertexMap->end())
            {
               child = new VertexCoder(vertex->edges[i]->source, vertex->edges[i], this, generation + 1, vertexMap);
               assert(child != NULL);
               children.push_back(child);
               (*vertexMap)[key] = child;
            }
            else
            {
               child = itr->second;

               // Share next generation.
               if (child->generation == generation + 1)
               {
                  children.push_back(child);
               }
            }
         }
      }
   }
   else
   {
      // Expand deeper.
      for (i = 0, i2 = (int)children.size(); i < i2; i++)
      {
         if (children[i]->creator == this)
         {
            children[i]->expand();
         }
      }
   }
}


// Print code.
void GraphHash::VertexCoder::printCode(FILE *fp)
{
   if (codeValid)
   {
      for (int i = 0; i < MD5_SIZE; i++)
      {
         fprintf(fp, "%d ", code[i]);
      }
      fprintf(fp, "\n");
   }
   else
   {
      fprintf(fp, "Code invalid\n");
   }
}


// Recursive dump in Graphviz "dot" format.
void GraphHash::VertexCoder::dump(FILE *fp)
{
   int i, i2;

   fprintf(fp, "\t\"%p\" [", this);
   if (vertex != NULL)
   {
      if (vertex->label != Graph::NULL_LABEL)
      {
         fprintf(fp, "label=\"%d\",", vertex->label);
      }
      else
      {
         fprintf(fp, "label=\"\",");
      }
   }
   else
   {
      fprintf(fp, "label=\"Root\",");
   }
   fprintf(fp, "shape=ellipse];\n");
   for (i = 0, i2 = (int)children.size(); i < i2; i++)
   {
      children[i]->dump(fp);
   }
   for (i = 0, i2 = (int)children.size(); i < i2; i++)
   {
      fprintf(fp, "\t\"%p\" -> \"%p\" [label=\"", this, children[i]);
      if (children[i]->parentEdge != NULL)
      {
         if (children[i]->parentEdge->label != Graph::NULL_LABEL)
         {
            fprintf(fp, "%d", children[i]->parentEdge->label);
         }
         if (children[i]->parentEdge->directed)
         {
            if (children[i]->parentEdge->source == vertex)
            {
               fprintf(fp, "(f)");
            }
            else
            {
               fprintf(fp, "(b)");
            }
         }
      }
      fprintf(fp, "\"];\n");
   }
}


// Less-than comparison by code.
bool GraphHash::VertexCoder::ltcmpCode(VertexCoder *a, VertexCoder *b)
{
   if (memcmp(a->code, b->code, MD5_SIZE) < 0)
   {
      return(true);
   }
   else
   {
      return(false);
   }
}
