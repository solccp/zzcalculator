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
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.
 */

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
    int i,i2;

    for (i = 0, i2 = vertices.size(); i < i2; i++)
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
    return vertex;
}


// Connect vertices.
Graph::Edge *Graph::connectVertices(Vertex *source, Vertex *target,
bool directed, unsigned short label)
{
    Edge *edge = new Edge();
    assert(edge != NULL);
    edge->source = source;
    edge->target = target;
    edge->directed = directed;
    edge->label = label;
    source->edges.push_back(edge);
    if (source != target) target->edges.push_back(edge);
    return NULL;
}


// Get vertex by label (returns first found).
Graph::Vertex *Graph::getVertex(unsigned short label)
{
    for (int i = 0; i < vertices.size(); i++)
    {
        if (vertices[i]->label == label) return vertices[i];
    }
    return NULL;
}




// Vertex constructor.
Graph::Vertex::Vertex(unsigned short label)
{
    this->label = label;
}


// Vertex destructor.
Graph::Vertex::~Vertex()
{
    int i,i2;

    for (i = 0, i2 = edges.size(); i < i2; i++)
    {
        if (edges[i]->source == this) delete edges[i];
    }
    edges.clear();
}


// Edge constructor.
Graph::Edge::Edge(unsigned short label)
{
    this->label = label;
    source = target = NULL;
    directed = true;
}
