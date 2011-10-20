// For conditions of distribution and use, see copyright notice in mona.hpp

/*
 * Graph.
 * The graph can have labeled or non-labeled vertices and edges.
 * An undirected graph contains a pair of bi-directional edges
 * between all connected vertices.
 */

#ifndef __GRAPH__
#define __GRAPH__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
using namespace std;

class Graph
{
public:

   // Null label.
   enum { NULL_LABEL=(unsigned short)(-1) };

   class Edge;

   // Graph vertex.
   class Vertex
   {
public:
      unsigned short label;
      vector<Edge *> edges;
      Vertex(unsigned short label = NULL_LABEL);
      void listEdges(vector<Edge *>& edgeList);
   };

   // Graph edge.
   class Edge
   {
public:
      unsigned short label;
      Vertex         *source;
      Vertex         *target;
      bool           directed;
      Edge(unsigned short label = NULL_LABEL);
   };

   // Vertices.
   vector<Vertex *> vertices;

   // Constructor.
   Graph();

   // Destructor.
   ~Graph();

   // Add vertex.
   Vertex *addVertex(unsigned short label = NULL_LABEL);

   // Connect vertices.
   Edge *connectVertices(Vertex *source, Vertex *target,
                         bool directed, unsigned short label = NULL_LABEL);

   // Get vertex by label (returns first found).
   Vertex *getVertex(unsigned short label);

   // Load and save.
};
#endif
