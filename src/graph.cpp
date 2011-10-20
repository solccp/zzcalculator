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
   return(edge);
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
