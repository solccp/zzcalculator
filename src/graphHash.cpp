// For conditions of distribution and use, see copyright notice in graphHash.hpp

/*
 * Compute an MD5 hash for a graph.
 */

#include "graphHash.hpp"

// Constructor.
GraphHash::GraphHash()
{
   memset(code, 0, MD5_SIZE);
   graph       = NULL;
   vertexCoder = NULL;
}


// Destructor.
GraphHash::~GraphHash()
{
   if (vertexCoder != NULL)
   {
      delete vertexCoder;
   }
   vertices.clear();
}


// Hash graph.
#ifdef THREADS
bool GraphHash::hash(Graph *graph, int numThreads, bool hashLabels)
#else
bool GraphHash::hash(Graph *graph, bool hashLabels)
#endif
{
   int         i, i2;
   VertexCoder *child;

   vector<int> vertexList;

   map<pair<Graph::Vertex *, Graph::Edge *>,
       VertexCoder *, VertexCoder::ltcmpConnection> *vertexMap;

   // Clear code.
   memset(code, 0, MD5_SIZE);

   // Save graph.
   this->graph = graph;

   // Sort vertices by label.
   vertices.clear();
   for (i = 0, i2 = (int)graph->vertices.size(); i < i2; i++)
   {
      vertices.push_back(graph->vertices[i]);
   }
   sort(vertices.begin(), vertices.end(), ltcmpVertices);

   // Generate code for graph using vertex coders.
   if (vertexCoder != NULL)
   {
      delete vertexCoder;
   }
   vertexCoder = new VertexCoder();
   assert(vertexCoder != NULL);
   if (vertexCoder == NULL)
   {
      return(false);
   }
   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      vertexMap = new map<pair<Graph::Vertex *, Graph::Edge *>,
                          GraphHash::VertexCoder *, VertexCoder::ltcmpConnection>;
      assert(vertexMap != NULL);
      if (vertexMap == NULL)
      {
         return(false);
      }
      child = new VertexCoder(vertices[i], NULL, vertexCoder, 1, vertexMap);
      assert(child != NULL);
      if (child == NULL)
      {
         delete vertexMap;
         return(false);
      }
      vertexCoder->children.push_back(child);
   }
#ifdef THREADS
   bool result = vertexCoder->generateCode(numThreads, hashLabels, &vertexList);
#else
   bool result = vertexCoder->generateCode(hashLabels, &vertexList);
#endif
   if (result)
   {
      memcpy(code, vertexCoder->code, MD5_SIZE);
   }
   return(result);
}


// Print graph and its hash.
void GraphHash::print(FILE *fp)
{
   if (graph == NULL)
   {
      return;
   }
   graph->print(fp);
   fprintf(fp, "Hash: ");
   for (int i = 0; i < MD5_SIZE; i++)
   {
      fprintf(fp, "%d ", code[i]);
   }
   fprintf(fp, "\n");
}


// Dump graph and graph hash in Graphvis "dot" format.
void GraphHash::dump(FILE *fp)
{
   char hash[MD5_SIZE * 5];

   if (graph == NULL)
   {
      return;
   }
   fprintf(fp, "digraph {\n");
   fprintf(fp, "\tgraph [size=\"8.5,11\",fontsize=14];\n");
   fprintf(fp, "\tsubgraph cluster_0 {\n");
   fprintf(fp, "\tlabel=\"Graph\";\n");
   graph->dumpSub(fp);
   fprintf(fp, "\t};\n");
   fprintf(fp, "\tsubgraph cluster_1 {\n");
   sprintf(hash, "Hash: ");
   for (int i = 0; i < MD5_SIZE; i++)
   {
      sprintf(&hash[strlen(hash)], "%d ", code[i]);
   }
   fprintf(fp, "\tlabel=\"%s\";\n", hash);
   vertexCoder->dump(fp);
   fprintf(fp, "\t};\n");
   fprintf(fp, "};\n");
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
#ifdef THREADS
bool GraphHash::VertexCoder::generateCode(int numThreads, bool hashLabels,
                                          vector<int> *vertexList)
{
   // Start additional update threads.
   this->numThreads = numThreads;
   terminate        = false;
   int  n      = 0;
   bool result = true;
   if (numThreads > 1)
   {
      if (pthread_barrier_init(&expandBarrier, NULL, numThreads) != 0)
      {
         return(false);
      }
      if (pthread_mutex_init(&expandMutex, NULL) != 0)
      {
         pthread_barrier_destroy(&expandBarrier);
         return(false);
      }
      threads = new pthread_t[numThreads - 1];
      assert(threads != NULL);
      if (threads == NULL)
      {
         pthread_mutex_destroy(&expandMutex);
         pthread_barrier_destroy(&expandBarrier);
         return(false);
      }
      struct ThreadInfo *info;
      for (int i = 0; i < numThreads - 1; i++)
      {
         info = new struct ThreadInfo;
         assert(info != NULL);
         if (info == NULL)
         {
            result = false;
            break;
         }
         info->coder      = this;
         info->threadNum  = i + 1;
         info->vertexList = vertexList;
         if (pthread_create(&threads[i], NULL, expandThread, (void *)info) != 0)
         {
            delete info;
            result = false;
            break;
         }
         n++;
      }
   }

   // Generate code.
   if (result)
   {
      result = generateCode(hashLabels, vertexList);
   }

   // Terminate threads.
   if (numThreads > 1)
   {
      // Unblock threads waiting on update barrier.
      terminate = true;
      expandVertices(0, vertexList);
      for (int i = 0; i < n; i++)
      {
         pthread_join(threads[i], NULL);
         pthread_detach(threads[i]);
      }
      pthread_mutex_destroy(&expandMutex);
      pthread_barrier_destroy(&expandBarrier);
      delete threads;
   }
   return(result);
}


// Vertex expander thread.
void *GraphHash::VertexCoder::expandThread(void *arg)
{
   struct ThreadInfo      *info  = (struct ThreadInfo *)arg;
   GraphHash::VertexCoder *coder = info->coder;
   int threadNum = info->threadNum;

   vector<int> *vertexList = info->vertexList;
   delete info;
   while (true)
   {
      coder->expandVertices(threadNum, vertexList);
   }
   return(NULL);
}


#endif

bool GraphHash::VertexCoder::generateCode(bool hashLabels, vector<int> *vertexList)
{
   int               i, j, s;
   struct MD5Context md5c;
   unsigned char     *input;

   // Root?
   s = (int)children.size();
   if (vertex == NULL)
   {
      // Fully expand vertices.
      while (true)
      {
         vertexList->clear();
         for (i = 0; i < s; i++)
         {
            if (!children[i]->expanded)
            {
               vertexList->push_back(i);
            }
         }
         if (vertexList->empty())
         {
            break;
         }
         else
         {
#ifdef THREADS
            if (!expandVertices(0, vertexList))
#else
            if (!expandVertices(vertexList))
#endif
            {
               return(false);
            }
         }
      }
   }

   // Recursively generate code.
   for (i = 0; i < s; i++)
   {
      if (!children[i]->generateCode(hashLabels))
      {
         return(false);
      }
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
   if (input == NULL)
   {
      return(false);
   }
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
            memcpy(&input[sizeof(vertex->label)],
                   &parentEdge->label, sizeof(parentEdge->label));
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
   codeValid = true;
   return(true);
}


// Expand vertices.
#ifdef THREADS
bool GraphHash::VertexCoder::expandVertices(int threadNum, vector<int> *vertexList)
#else
bool GraphHash::VertexCoder::expandVertices(vector<int> *vertexList)
#endif
{
#ifdef THREADS
   if (numThreads > 1)
   {
      // Synchronize threads.
      if (threadNum == 0)
      {
         expandResult = true;
      }
      int i = pthread_barrier_wait(&expandBarrier);
      if ((i != PTHREAD_BARRIER_SERIAL_THREAD) && (i != 0))
      {
         pthread_exit(NULL);
      }
      if (terminate)
      {
         if (threadNum == 0)
         {
            return(true);
         }
         pthread_exit(NULL);
      }

      // Update vertices.
      while (true)
      {
         pthread_mutex_lock(&expandMutex);
         int v = -1;
         if (expandResult && !vertexList->empty())
         {
            v = vertexList->back();
            vertexList->pop_back();
         }
         pthread_mutex_unlock(&expandMutex);
         if (v == -1)
         {
            break;
         }
         if (!children[v]->expand())
         {
            expandResult = false;
         }
      }

      // Re-group threads.
      pthread_barrier_wait(&expandBarrier);

      return(expandResult);
   }
   else
   {
#endif
   for (int i = 0, j = (int)vertexList->size(); i < j; i++)
   {
      if (!children[(*vertexList)[i]]->expand() ||
          !children[(*vertexList)[i]]->generateCode(false))
      {
         return(false);
      }
   }
   return(true);

#ifdef THREADS
}
#endif


}


// Expand coder one level deeper.
bool GraphHash::VertexCoder::expand()
{
   int         i, i2;
   VertexCoder *child;

   pair<Graph::Vertex *, Graph::Edge *> key;
   map<pair<Graph::Vertex *, Graph::Edge *>,
       VertexCoder *, ltcmpConnection>::iterator itr;

   if (expanded)
   {
      return(true);
   }
   expanded = true;
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
               child = new VertexCoder(vertex->edges[i]->target,
                                       vertex->edges[i], this, generation + 1, vertexMap);
               assert(child != NULL);
               if (child == NULL)
               {
                  return(false);
               }
               children.push_back(child);
               (*vertexMap)[key] = child;
               expanded          = false;
            }
            else
            {
               child = itr->second;

               // Share next generation.
               if (child->generation == generation + 1)
               {
                  children.push_back(child);
                  expanded = false;
               }
            }
         }
         else
         {
            key.first  = vertex->edges[i]->source;
            key.second = vertex->edges[i];
            if ((itr = vertexMap->find(key)) == vertexMap->end())
            {
               child = new VertexCoder(vertex->edges[i]->source,
                                       vertex->edges[i], this, generation + 1, vertexMap);
               assert(child != NULL);
               if (child == NULL)
               {
                  return(false);
               }
               children.push_back(child);
               (*vertexMap)[key] = child;
               expanded          = false;
            }
            else
            {
               child = itr->second;

               // Share next generation.
               if (child->generation == generation + 1)
               {
                  children.push_back(child);
                  expanded = false;
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
            if (!children[i]->expand())
            {
               return(false);
            }
         }
         if (!children[i]->expanded)
         {
            expanded = false;
         }
      }
   }
   return(true);
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
