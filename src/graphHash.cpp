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
bool GraphHash::hash(Graph *graph,
                     int numThreads, bool hashLabels)
#else
bool GraphHash::hash(Graph *graph, bool hashLabels)
#endif
{
   int             i, i2;
   VertexCoder     *child;
   VertexCoderLink *link;

   vector<int> vertexList;

   map<pair<Graph::Vertex *, bool>, VertexCoder *,
       VertexCoder::ltcmpVertices> *vertexMap;

   // Clear code.
   memset(code, 0, MD5_SIZE);

   // Save graph.
   this->graph = graph;

   // Sort vertices by label.
   vertices.clear();
   if ((int)graph->vertices.size() == 0)
   {
      return(false);
   }
   for (i = 0, i2 = (int)graph->vertices.size(); i < i2; i++)
   {
      vertices.push_back(graph->vertices[i]);
   }
   sort(vertices.begin(), vertices.end(), ltcmpVertexLabels);

   // Generate code for graph using vertex coders.
   if (vertexCoder != NULL)
   {
      delete vertexCoder;
   }
   vertexCoder = new VertexCoder();
   if (vertexCoder == NULL)
   {
      return(false);
   }
   for (i = 0, i2 = (int)vertices.size(); i < i2; i++)
   {
      vertexMap = new map<pair<Graph::Vertex *, bool>,
                          GraphHash::VertexCoder *,
                          VertexCoder::ltcmpVertices>;
      if (vertexMap == NULL)
      {
         return(false);
      }
      child = new VertexCoder(vertices[i],
                              false, 1, vertexMap);
      if (child == NULL)
      {
         delete vertexMap;
         return(false);
      }
      link = new VertexCoderLink(NULL, child, true);
      if (link == NULL)
      {
         delete vertexMap;
         delete child;
         return(false);
      }
      vertexCoder->children.push_back(link);
   }
#ifdef THREADS
   bool result = vertexCoder->generateCode(numThreads,
                                           hashLabels, &vertexList);
#else
   bool result =
      vertexCoder->generateCode(hashLabels, &vertexList);
#endif
   if (result)
   {
      memcpy(code, vertexCoder->code, MD5_SIZE);
   }
   return(result);
}


// Get graph hash.
// Valid after graph hash.
// Hash size is MD5_SIZE bytes (see md5.h).
unsigned char *GraphHash::getHash()
{
   if (vertexCoder != NULL)
   {
      return(vertexCoder->code);
   }
   else
   {
      return(NULL);
   }
}



// Less-than comparison of vertices by label.
bool GraphHash::ltcmpVertexLabels(Graph::Vertex *a, Graph::Vertex *b)
{
   return(a->label < b->label);
}


// Vertex coder link constructors.
GraphHash::VertexCoderLink::VertexCoderLink()
{
   edge    = NULL;
   coder   = NULL;
   creator = false;
}


GraphHash::VertexCoderLink::VertexCoderLink(
   Graph::Edge *edge, VertexCoder *coder, bool creator)
{
   this->edge    = edge;
   this->coder   = coder;
   this->creator = creator;
}


// Vertex coder constructors.
GraphHash::VertexCoder::VertexCoder()
{
   vertex     = NULL;
   expanded   = false;
   halo       = false;
   generation = 0;
   codeValid  = false;
   vertexMap  = NULL;
}


GraphHash::VertexCoder::VertexCoder(Graph::Vertex *vertex,
                                    bool halo, int generation,
                                    map<pair<Graph::Vertex *, bool>,
                                        VertexCoder *,
                                        ltcmpVertices> *vertexMap)
{
   pair<Graph::Vertex *, bool> key;

   this->vertex      = vertex;
   this->halo        = halo;
   this->generation  = generation;
   this->vertexMap   = vertexMap;
   key.first         = vertex;
   key.second        = halo;
   (*vertexMap)[key] = this;
   expanded          = false;
   codeValid         = false;
}


// Vertex coder destructor.
GraphHash::VertexCoder::~VertexCoder()
{
   int i, i2;

   // Purge children from non-creating coders to
   // prevent duplicate deletions.
   if (vertex == NULL)
   {
      for (i = 0, i2 = (int)children.size(); i < i2; i++)
      {
         if (children[i]->coder->vertexMap != NULL)
         {
            delete children[i]->coder->vertexMap;
         }
         children[i]->coder->purgeChildren();
      }
   }
   for (i = 0, i2 = (int)children.size(); i < i2; i++)
   {
      if (children[i] != NULL)
      {
         delete children[i]->coder;
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
      if (children[i]->creator)
      {
         children[i]->coder->purgeChildren();
      }
      else
      {
         delete children[i];
         children[i] = NULL;
      }
   }
}


// Generate code.
#ifdef THREADS
bool GraphHash::VertexCoder::generateCode(int         numThreads,
                                          bool        hashLabels,
                                          vector<int> *vertexList)
{
   // Start additional update threads.
   this->numThreads = numThreads;
   terminate        = false;
   int  n      = 0;
   bool result = true;
   if (numThreads > 1)
   {
      if (pthread_barrier_init(&updateBarrier, NULL,
                               numThreads) != 0)
      {
         return(false);
      }
      if (pthread_mutex_init(&updateMutex, NULL) != 0)
      {
         pthread_barrier_destroy(&updateBarrier);
         return(false);
      }
      threads = new pthread_t[numThreads - 1];
      if (threads == NULL)
      {
         pthread_mutex_destroy(&updateMutex);
         pthread_barrier_destroy(&updateBarrier);
         return(false);
      }
      struct ThreadInfo *info;
      for (int i = 0; i < numThreads - 1; i++)
      {
         info = new struct ThreadInfo;
         if (info == NULL)
         {
            result = false;
            break;
         }
         info->coder      = this;
         info->threadNum  = i + 1;
         info->vertexList = vertexList;
         if (pthread_create(&threads[i], NULL, updateThread,
                            (void *)info) != 0)
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
      updateVertices(0, vertexList);
      for (int i = 0; i < n; i++)
      {
         pthread_join(threads[i], NULL);
         pthread_detach(threads[i]);
      }
      pthread_mutex_destroy(&updateMutex);
      pthread_barrier_destroy(&updateBarrier);
      delete threads;
   }
   return(result);
}


// Vertex update thread.
void *GraphHash::VertexCoder::updateThread(void *arg)
{
   struct ThreadInfo      *info  = (struct ThreadInfo *)arg;
   GraphHash::VertexCoder *coder = info->coder;
   int threadNum = info->threadNum;

   vector<int> *vertexList = info->vertexList;
   delete info;
   while (true)
   {
      coder->updateVertices(threadNum, vertexList);
   }
   return(NULL);
}


#endif

bool GraphHash::VertexCoder::generateCode(bool hashLabels, vector<int> *vertexList)
{
   int               i, j, k, n, s;
   struct MD5Context md5c;
   unsigned char     *input, *p;

   vector<Graph::Edge *> edges;

   // Root?
   s = (int)children.size();
   if (vertex == NULL)
   {
      // Iteratively expand children and generate codes without
      // labels to encode structure.
      vertexList->clear();
      for (i = 0; i < s; i++)
      {
         vertexList->push_back(i);
      }
#ifdef THREADS
      if (!updateVertices(0, vertexList))
#else
      if (!updateVertices(vertexList))
#endif
      {
         return(false);
      }

      // Continue until fully expanded or all hashes are distinct.
      while (true)
      {
         // All vertices are fully expanded?
         for (i = 0; i < s; i++)
         {
            if (!children[i]->coder->expanded)
            {
               break;
            }
         }
         if (i == s)
         {
            break;
         }

         // All vertices are distinct?
         sort(children.begin(), children.end(), ltcmpCodes);
         i = 0;
         j = s - 1;
         vertexList->clear();
         while (true)
         {
            for ( ; i < j; i++)
            {
               if ((!children[i]->coder->expanded ||
                    !children[i + 1]->coder->expanded) &&
                   (ltcmpCodes(children[i], children[i + 1]) == 0))
               {
                  break;
               }
            }
            if (i == j)
            {
               break;
            }

            // Expand and generate codes for range of equivalent children.
            for (k = i + 1; k < s &&
                 ltcmpCodes(children[i], children[k]) == 0; k++)
            {
            }
            for ( ; i < k; i++)
            {
               if (!children[i]->coder->expanded)
               {
                  vertexList->push_back(i);
               }
            }
            if (i == s)
            {
               break;
            }
         }
         if ((int)vertexList->empty())
         {
            // All distinct.
            break;
         }
#ifdef THREADS
         if (!updateVertices(0, vertexList))
#else
         if (!updateVertices(vertexList))
#endif
         {
            return(false);
         }
      }
   }

   // Recursively generate code.
   for (i = 0; i < s; i++)
   {
      if (!children[i]->coder->generateCode(hashLabels))
      {
         return(false);
      }
   }
   sort(children.begin(), children.end(), ltcmpCodes);
   MD5Init(&md5c);
   n = 0;
   if (vertex != NULL)
   {
      if (hashLabels)
      {
         n += sizeof(vertex->label);
      }
      if ((int)children.size() > 0)
      {
         if (hashLabels)
         {
            for (i = 0; i < (int)children.size(); i++)
            {
               n += sizeof(children[i]->edge->label);
            }
         }
         n += (int)children.size();
      }
   }
   input = new unsigned char[(s * MD5_SIZE) + n];
   if (input == NULL)
   {
      return(false);
   }
   if (vertex != NULL)
   {
      p = input;
      if (hashLabels)
      {
         memcpy(input, &vertex->label, sizeof(vertex->label));
         p += sizeof(vertex->label);
      }
      if ((int)children.size() > 0)
      {
         for (i = 0; i < (int)children.size(); i++)
         {
            edges.push_back(children[i]->edge);
         }
         sortEdges(edges, vertex, hashLabels);
         for (i = 0; i < (int)edges.size(); i++)
         {
            if (hashLabels)
            {
               memcpy(p, &edges[i]->label, sizeof(edges[i]->label));
               p += sizeof(edges[i]->label);
            }
            if (edges[i]->directed)
            {
               if (edges[i]->source == vertex)
               {
                  *p = 0;
               }
               else
               {
                  *p = 1;
               }
            }
            else
            {
               *p = 2;
            }
            p++;
         }
      }
   }
   for (i = 0; i < s; i++)
   {
      memcpy(&input[n], children[i]->coder->code, MD5_SIZE);
      n += MD5_SIZE;
   }
   MD5Update(&md5c, input, n);
   MD5Final(code, &md5c);
   delete [] input;
   codeValid = true;
   return(true);
}


// Sort edges.
void GraphHash::VertexCoder::sortEdges(vector<Graph::Edge *>& edges,
                                       Graph::Vertex          *vertex,
                                       bool                   hashLabels)
{
   int         n = (int)edges.size();
   int         iPos, iMin;
   int         d1, d2;
   Graph::Edge *tmpEdge;

   for (iPos = 0; iPos < n; iPos++)
   {
      iMin = iPos;
      if (edges[iMin]->directed)
      {
         if (edges[iMin]->source == vertex)
         {
            d1 = 0;
         }
         else
         {
            d1 = 1;
         }
      }
      else
      {
         d1 = 2;
      }
      for (int i = iPos + 1; i < n; i++)
      {
         if (edges[i]->directed)
         {
            if (edges[i]->source == vertex)
            {
               d2 = 0;
            }
            else
            {
               d2 = 1;
            }
         }
         else
         {
            d2 = 2;
         }
         if (d2 < d1)
         {
            iMin = i;
         }
         else if (d2 == d1)
         {
            if (hashLabels && (edges[i]->label < edges[iMin]->label))
            {
               iMin = i;
            }
         }
      }

      if (iMin != iPos)
      {
         tmpEdge     = edges[iPos];
         edges[iPos] = edges[iMin];
         edges[iMin] = tmpEdge;
      }
   }
}


// Update vertices.
#ifdef THREADS
bool GraphHash::VertexCoder::updateVertices(int         threadNum,
                                            vector<int> *vertexList)
#else
bool GraphHash::VertexCoder::updateVertices(vector<int> *vertexList)
#endif
{
#ifdef THREADS
   if (numThreads > 1)
   {
      // Synchronize threads.
      if (threadNum == 0)
      {
         updateResult = true;
      }
      int i = pthread_barrier_wait(&updateBarrier);
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
         pthread_mutex_lock(&updateMutex);
         int v = -1;
         if (updateResult && !vertexList->empty())
         {
            v = vertexList->back();
            vertexList->pop_back();
         }
         pthread_mutex_unlock(&updateMutex);
         if (v == -1)
         {
            break;
         }
         if (!children[v]->coder->expand() ||
             !children[v]->coder->generateCode(false))
         {
            updateResult = false;
         }
      }

      // Re-group threads.
      pthread_barrier_wait(&updateBarrier);

      return(updateResult);
   }
   else
   {
#endif
   for (int i = 0, j = (int)vertexList->size(); i < j; i++)
   {
      if (!children[(*vertexList)[i]]->coder->expand() ||
          !children[(*vertexList)[i]]->coder->generateCode(false))
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
   int i, i2;

   if (expanded)
   {
      return(true);
   }
   expanded = true;
   if (children.empty())
   {
      for (i = 0, i2 = (int)vertex->edges.size(); i < i2; i++)
      {
         if (vertex == vertex->edges[i]->source)
         {
            if (!expand(vertex->edges[i]->target,
                        vertex->edges[i], halo))
            {
               return(false);
            }
            if (expanded)
            {
               if (!halo)
               {
                  if (!expand(vertex->edges[i]->target,
                              vertex->edges[i], true))
                  {
                     return(false);
                  }
               }
            }
         }
         else
         {
            if (!expand(vertex->edges[i]->source,
                        vertex->edges[i], halo))
            {
               return(false);
            }
            if (expanded)
            {
               if (!halo)
               {
                  if (!expand(vertex->edges[i]->source,
                              vertex->edges[i], true))
                  {
                     return(false);
                  }
               }
            }
         }
      }
   }
   else
   {
      for (i = 0, i2 = (int)children.size(); i < i2; i++)
      {
         if (!children[i]->coder->halo)
         {
            if (!children[i]->coder->expanded &&
                (children[i]->creator))
            {
               if (!children[i]->coder->expand())
               {
                  return(false);
               }
            }
            if (!children[i]->coder->expanded)
            {
               expanded = false;
            }
         }
      }
   }
   return(true);
}


// Expand to child.
bool GraphHash::VertexCoder::expand(Graph::Vertex *childVertex,
                                    Graph::Edge *edge, bool halo)
{
   VertexCoder     *child;
   VertexCoderLink *link;

   pair<Graph::Vertex *, bool> key;
   map<pair<Graph::Vertex *, bool>,
       VertexCoder *, ltcmpVertices>::iterator itr;

   key.first  = childVertex;
   key.second = halo;
   if ((itr = vertexMap->find(key)) == vertexMap->end())
   {
      child = new VertexCoder(childVertex, halo,
                              generation + 1, vertexMap);
      if (child == NULL)
      {
         return(false);
      }
      link = new VertexCoderLink(edge, child, true);
      if (link == NULL)
      {
         return(false);
      }
      children.push_back(link);
      if (!halo)
      {
         expanded = false;
      }
      else
      {
         child->expanded = true;
      }
   }
   else
   {
      child = itr->second;
      if (child->generation == generation + 1)
      {
         link = new VertexCoderLink(edge, child, false);
         if (link == NULL)
         {
            return(false);
         }
         children.push_back(link);
      }
   }
   return(true);
}


// Less-than comparison of codes.
bool GraphHash::VertexCoder::ltcmpCodes(VertexCoderLink *a, VertexCoderLink *b)
{
   if (memcmp(a->coder->code, b->coder->code, MD5_SIZE) < 0)
   {
      return(true);
   }
   else
   {
      return(false);
   }
}
