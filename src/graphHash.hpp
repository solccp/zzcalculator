/*
 * This software is provided under the terms of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * Copyright (c) 2007-2011 Tom Portegys, All Rights Reserved.
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

#ifndef __GRAPH_HASH__
#define __GRAPH_HASH__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include "graph.hpp"
#include "md5.h"
#include <cstring>
using namespace std;

class GraphHash
{
public:

   // Graph.
   Graph *graph;

   // Hash code.
   unsigned char code[MD5_SIZE];

   // Constructor.
   GraphHash();

   // Destructor.
   ~GraphHash();

   // Hash graph.
#ifdef THREADS
   bool hash(Graph *graph, int numThreads = 1,
             bool hashLabels = true);

#else
   bool hash(Graph *graph, bool hashLabels = true);
#endif

   // Get graph hash.
   // Valid after graph hash.
   // Hash size is MD5_SIZE bytes (see md5.h).
   unsigned char *getHash();


   // Vertices.
   vector<Graph::Vertex *> vertices;

   // Less-than comparison of vertices by label.
   static bool ltcmpVertexLabels(Graph::Vertex *a, Graph::Vertex *b);

   class VertexCoder;

   // Vertex coder link.
   class VertexCoderLink
   {
public:
      Graph::Edge *edge;
      VertexCoder *coder;
      bool        creator;
      VertexCoderLink();
      VertexCoderLink(Graph::Edge *edge,
                      VertexCoder *coder, bool creator);
   };

   // Vertex coder.
   class VertexCoder
   {
public:
      Graph::Vertex             *vertex;
      unsigned char             code[MD5_SIZE];
      vector<VertexCoderLink *> children;
      struct ltcmpVertices
      {
         bool operator()(pair<Graph::Vertex *, bool> a,
                         pair<Graph::Vertex *, bool> b) const
         {
            if (a.first < b.first) { return(true); }
            else if (a.first == b.first)
            {
               if (a.second && !b.second) { return(true); }
            }
            return(false);
         }
      };
      VertexCoder();
      VertexCoder(Graph::Vertex *vertex,
                  bool halo, int generation,
                  map<pair<Graph::Vertex *, bool>,
                      VertexCoder *,
                      ltcmpVertices> *vertexMap = NULL);
      ~VertexCoder();
#ifdef THREADS
      bool generateCode(int numThreads, bool hashLabels,
                        vector<int> *vertexList = NULL);
#endif
      bool generateCode(bool        hashLabels,
                        vector<int> *vertexList = NULL);

private:
#ifdef THREADS
      bool updateVertices(int threadNum, vector<int> *vertexList);

      bool updateResult;
      static void *updateThread(void *threadInfo);

#else
      bool updateVertices(vector<int> *vertexList);
#endif
      bool expanded;
      bool halo;
      int  generation;
      bool codeValid;
      bool expand();
      bool expand(Graph::Vertex *childVertex,
                  Graph::Edge *, bool halo);
      void purgeChildren();
      void sortEdges(vector<Graph::Edge *>&,
                     Graph::Vertex *, bool hashLabels);

      map<pair<Graph::Vertex *, bool>,
          VertexCoder *, ltcmpVertices> *vertexMap;

      static bool ltcmpCodes(VertexCoderLink *a, VertexCoderLink *b);

#ifdef THREADS
      int               numThreads;
      bool              terminate;
      pthread_barrier_t updateBarrier;
      pthread_mutex_t   updateMutex;
      pthread_t         *threads;
      struct ThreadInfo
      {
         VertexCoder *coder;
         int         threadNum;
         vector<int> *vertexList;
      };
#endif
   };

   // Root vertex coder.
   VertexCoder *vertexCoder;
};
#endif
