// For conditions of distribution and use, see copyright notice in graphHash.hpp

/*
 * Graph.
 */

#include "graph.hpp"

// Constructor.
Graph::Graph()
{
}


// Destructor.
Graph::~Graph()
{
   int i, i2;

   vector<Edge *> edgeList;

   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      vertices[i]->listEdges(edgeList);
   }
   for (i = 0, i2 = (int)edgeList.size(); i < i2; i++)
   {
      delete edgeList[i];
   }
   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      delete vertices[i];
   }
   vertices.clear();
}


// Add vertex.
Graph::Vertex *Graph::addVertex(unsigned short label)
{
   Vertex *vertex = new Vertex(label);

   assert(vertex != NULL);
   vertices.push_back(vertex);
   return(vertex);
}


// Connect vertices.
Graph::Edge *Graph::connectVertices(Vertex *source, Vertex *target,
                                    bool directed, unsigned short label)
{
   Edge *edge = new Edge();

   assert(edge != NULL);
   edge->source   = source;
   edge->target   = target;
   edge->directed = directed;
   edge->label    = label;
   source->edges.push_back(edge);
   if (source != target)
   {
      target->edges.push_back(edge);
   }
   return(NULL);
}


// Get vertex by label (returns first found).
Graph::Vertex *Graph::getVertex(unsigned short label)
{
   int i, i2;

   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      if (vertices[i]->label == label)
      {
         return(vertices[i]);
      }
   }
   return(NULL);
}


// Load graph.
void Graph::load(FILE *fp)
{
   int            i, i2, j, j2, k;
   unsigned short label;
   Vertex         *vertex;
   Edge           *edge;

   FREAD_INT(&i2, fp);
   for (i = 0; i < i2; i++)
   {
      FREAD_SHORT((short *)&label, fp);
      vertex = new Vertex(label);
      assert(vertex != NULL);
      vertices.push_back(vertex);
   }
   for (i = 0; i < i2; i++)
   {
      FREAD_INT(&j2, fp);
      for (j = 0; j < j2; j++)
      {
         FREAD_SHORT((short *)&label, fp);
         edge = new Edge(label);
         assert(edge != NULL);
         vertices[i]->edges.push_back(edge);
         FREAD_INT(&k, fp);
         edge->source = vertices[k];
         FREAD_INT(&k, fp);
         edge->target = vertices[k];
         FREAD_BOOL(&edge->directed, fp);
      }
   }
}


// Save graph.
void Graph::save(FILE *fp)
{
   int i, i2, j, j2, k;

   i2 = (int)vertices.size();
   FWRITE_INT(&i2, fp);
   for (i = 0; i < i2; i++)
   {
      FWRITE_SHORT((short *)&vertices[i]->label, fp);
   }
   for (i = 0; i < i2; i++)
   {
      j2 = (int)vertices[i]->edges.size();
      FWRITE_INT(&j2, fp);
      for (j = 0; j < j2; j++)
      {
         FWRITE_SHORT((short *)&vertices[i]->edges[j]->label, fp);
         for (k = 0; k < i2; k++)
         {
            if (vertices[k] == vertices[i]->edges[j]->source)
            {
               break;
            }
         }
         FWRITE_INT(&k, fp);
         for (k = 0; k < i2; k++)
         {
            if (vertices[k] == vertices[i]->edges[j]->target)
            {
               break;
            }
         }
         FWRITE_INT(&k, fp);
         FWRITE_BOOL(&vertices[i]->edges[j]->directed, fp);
      }
   }
}


// Print graph.
void Graph::print(FILE *fp)
{
   print(NULL, fp);
}


void Graph::print(char *label, FILE *fp)
{
   int i, i2, j, j2;

   if (label == NULL)
   {
      fprintf(fp, "Graph:\n");
   }
   else
   {
      fprintf(fp, "Graph: %s\n", label);
   }
   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      if (vertices[i]->label != NULL_LABEL)
      {
         fprintf(fp, "vertex %p, label = %d\n", vertices[i], vertices[i]->label);
      }
      else
      {
         fprintf(fp, "vertex %p\n", vertices[i]);
      }
      for (j = 0, j2 = (int)vertices[i]->edges.size(); j < j2; j++)
      {
         if (vertices[i]->edges[j]->directed)
         {
            if (vertices[i]->edges[j]->source == vertices[i])
            {
               if (vertices[i]->edges[j]->label != NULL_LABEL)
               {
                  fprintf(fp, "\t%p -> edge %p, label = %d -> %p\n", vertices[i]->edges[j]->source,
                          vertices[i]->edges[j], vertices[i]->edges[j]->label, vertices[i]->edges[j]->target);
               }
               else
               {
                  fprintf(fp, "\t%p -> edge %p -> %p\n", vertices[i]->edges[j]->source,
                          vertices[i]->edges[j], vertices[i]->edges[j]->target);
               }
            }
            else
            {
               if (vertices[i]->edges[j]->label != NULL_LABEL)
               {
                  fprintf(fp, "\t%p <- edge %p, label = %d <- %p\n", vertices[i]->edges[j]->target,
                          vertices[i]->edges[j], vertices[i]->edges[j]->label, vertices[i]->edges[j]->source);
               }
               else
               {
                  fprintf(fp, "\t%p <- edge %p <- %p\n", vertices[i]->edges[j]->target,
                          vertices[i]->edges[j], vertices[i]->edges[j]->source);
               }
            }
         }
         else
         {
            if (vertices[i]->edges[j]->label != NULL_LABEL)
            {
               fprintf(fp, "\t%p <- edge %p, label = %d -> %p\n", vertices[i]->edges[j]->source,
                       vertices[i]->edges[j], vertices[i]->edges[j]->label, vertices[i]->edges[j]->target);
            }
            else
            {
               fprintf(fp, "\t%p <- edge %p -> %p\n", vertices[i]->edges[j]->source,
                       vertices[i]->edges[j], vertices[i]->edges[j]->target);
            }
         }
      }
   }
}


// Dump graph in Graphviz "dot" format.
void Graph::dump(FILE *fp)
{
   dump(NULL, fp);
}


void Graph::dump(char *label, FILE *fp)
{
   fprintf(fp, "digraph {\n");
   if (label == NULL)
   {
      fprintf(fp, "\tgraph [size=\"8.5,11\",fontsize=14];\n");
   }
   else
   {
      fprintf(fp, "\tgraph [size=\"8.5,11\",fontsize=14,label=\"%s\"];\n", label);
   }
   dumpSub(fp);
   fprintf(fp, "};\n");
}


// Dump "guts" of graph.
void Graph::dumpSub(FILE *fp)
{
   int i, i2, j, j2;

   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      if (vertices[i]->label != NULL_LABEL)
      {
         fprintf(fp, "\t\"%p\" [label=%d,shape=ellipse];\n", vertices[i], vertices[i]->label);
      }
      else
      {
         fprintf(fp, "\t\"%p\" [label=\"\",shape=ellipse];\n", vertices[i]);
      }
   }
   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      for (j = 0, j2 = (int)vertices[i]->edges.size(); j < j2; j++)
      {
         if (vertices[i]->edges[j]->directed)
         {
            if (vertices[i]->edges[j]->source == vertices[i])
            {
               if (vertices[i]->edges[j]->label != NULL_LABEL)
               {
                  fprintf(fp, "\t\"%p\" -> \"%p\" [label=%d];\n",
                          vertices[i]->edges[j]->source, vertices[i]->edges[j]->target,
                          vertices[i]->edges[j]->label);
               }
               else
               {
                  fprintf(fp, "\t\"%p\" -> \"%p\";\n",
                          vertices[i]->edges[j]->source, vertices[i]->edges[j]->target);
               }
            }
         }
         else
         {
            if (vertices[i]->edges[j]->source == vertices[i])
            {
               if (vertices[i]->edges[j]->label != NULL_LABEL)
               {
                  fprintf(fp, "\t\"%p\" -> \"%p\" [label=%d,dir=both];\n",
                          vertices[i]->edges[j]->source, vertices[i]->edges[j]->target,
                          vertices[i]->edges[j]->label);
               }
               else
               {
                  fprintf(fp, "\t\"%p\" -> \"%p\" [dir=both];\n",
                          vertices[i]->edges[j]->source, vertices[i]->edges[j]->target);
               }
            }
         }
      }
   }
}


// Vertex constructor.
Graph::Vertex::Vertex(unsigned short label)
{
   this->label = label;
}


// List vertex edges.
void Graph::Vertex::listEdges(vector<Graph::Edge *>& edgeList)
{
   int i, i2;

   for (i = 0, i2 = (int)edges.size(); i < i2; i++)
   {
      if (edges[i]->source == this)
      {
         edgeList.push_back(edges[i]);
      }
   }
   edges.clear();
}


// Edge constructor.
Graph::Edge::Edge(unsigned short label)
{
   this->label = label;
   source      = target = NULL;
   directed    = true;
}
