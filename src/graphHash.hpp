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
   bool hash(Graph *graph, int numThreads = 1, bool hashLabels = true);

#else
   bool hash(Graph *graph, bool hashLabels = true);
#endif

   // Print graph and its hash.
   void print(FILE *fp = stdout);

   // Dump graph and graph hash in Graphvis "dot" format.
   void dump(FILE *fp = stdout);

   // Vertices.
   vector<Graph::Vertex *> vertices;

   // Less-than comparison of vertices by label.
   static bool ltcmpVertices(Graph::Vertex *a, Graph::Vertex *b);

   // Vertex coder.
   class VertexCoder
   {
public:
      Graph::Vertex         *vertex;
      Graph::Edge           *parentEdge;
      unsigned char         code[MD5_SIZE];
      bool                  codeValid;
      vector<VertexCoder *> children;
      VertexCoder           *creator;
      struct ltcmpConnection
      {
         bool operator()(pair<Graph::Vertex *, Graph::Edge *> a,
                         pair<Graph::Vertex *, Graph::Edge *> b) const
         {
            if (a.first < b.first) { return(true); }
            if ((a.first == b.first) && (a.second < b.second)) { return(true); }
            return(false);
         }
      };
      VertexCoder();
      VertexCoder(Graph::Vertex *vertex, Graph::Edge *parentEdge,
                  VertexCoder *creator, int generation,
                  map<pair<Graph::Vertex *, Graph::Edge *>,
                      VertexCoder *, ltcmpConnection> *vertexMap = NULL);
      ~VertexCoder();
#ifdef THREADS
      bool generateCode(int numThreads, bool hashLabels = true,
                        vector<int> *vertexList = NULL);
#endif
      bool generateCode(bool hashLabels = true, vector<int> *vertexList = NULL);
      void printCode(FILE *fp = stdout);
      void dump(FILE *fp = stdout);

private:
#ifdef THREADS
      bool expandVertices(int threadNum, vector<int> *vertexList);

      bool expandResult;
      static void *expandThread(void *threadInfo);

#else
      bool expandVertices(vector<int> *vertexList);
#endif
      bool expanded;
      int  generation;
      bool expand();

      map<pair<Graph::Vertex *, Graph::Edge *>,
          VertexCoder *, ltcmpConnection> *vertexMap;
      static bool ltcmpCode(VertexCoder *a, VertexCoder *b);
      void purgeChildren();

#ifdef THREADS
      int               numThreads;
      bool              terminate;
      pthread_barrier_t expandBarrier;
      pthread_mutex_t   expandMutex;
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
