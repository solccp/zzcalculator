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
                Vertex(unsigned short label=NULL_LABEL);
                ~Vertex();
        };

        // Graph edge.
        class Edge
        {
            public:
                unsigned short label;
                Vertex *source;
                Vertex *target;
                bool directed;
                Edge(unsigned short label=NULL_LABEL);
        };

        // Vertices.
        vector<Vertex *> vertices;

        // Constructor.
        Graph();

        // Destructor.
        ~Graph();

        // Add vertex.
        Vertex *addVertex(unsigned short label=NULL_LABEL);

        // Connect vertices.
        Edge *connectVertices(Vertex *source, Vertex *target,
            bool directed, unsigned short label=NULL_LABEL);

        // Get vertex by label (returns first found).
        Vertex *getVertex(unsigned short label);

        // Load and save.
        void load(FILE *fp);
        void save(FILE *fp);

        // Print.
        void print(FILE *fp=stdout);
        void print(char *label, FILE *fp=stdout);

        // Dump graph in Graphviz "dot" format.
        void dump(FILE *fp=stdout);
        void dump(char *label, FILE *fp=stdout);
        void dumpSub(FILE *fp=stdout);
};
#endif
