
#include <iostream>
#include "stdint.h"
#include <cstring>
#include <cstdio>
#include "graphHash.hpp"

using namespace std;



extern "C"
{
    void make_hash_(int32_t* nat, int32_t* num_of_edges, int32_t* edges, char* key)
    {
        char tmp_key[3];

        Graph *graph = new Graph();
        for(int i=0; i< *nat; ++i)
        {
            graph->addVertex((unsigned short)i);
        }
        for (int i=0; i< *num_of_edges ; ++i)
        {
            Graph::Vertex *v1,*v2;
            v1 = graph->vertices[edges[i*2]-1];
            v2 = graph->vertices[edges[i*2+1]-1];
            graph->connectVertices(v1, v2, false);
        }
        GraphHash *hash = new GraphHash(graph, false);

        for (int i=0 ; i< MD5_SIZE; ++i)
        {
            sprintf(tmp_key, "%02X", hash->code[i]);
            key[i*2]=tmp_key[0];
            key[i*2+1]=tmp_key[1];
        }

        delete hash;
        delete graph;
    }
};
