

#ifndef __HASH_DATABASE__
#define __HASH_DATABASE__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "polynomial.hpp"
#include "md5.h"
#include "graphHash.hpp"

#include <boost/format.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>


#include <boost/tuple/tuple.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>



typedef std::vector<std::pair<int, int> > EdgeList;
typedef boost::array<double,3> Vector3;
typedef std::vector<Vector3> GeomCoord;
const double CC_Tol = 1.4;

using boost::multi_index_container;
using namespace boost::multi_index;


class hash_entry
{
public:
    friend class boost::serialization::access;
    std::string hash;
    Polynomial poly;
    int natoms;
    EdgeList edges;

    hash_entry()
    {
        natoms = 0;
    }

    hash_entry(int nat, const EdgeList& edgelist, const Polynomial& py)
    {
        natoms = nat;
        edges = edgelist;
        poly = py;
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int file_version)
    {
        ar & hash & natoms & edges & poly ;
    }
};

typedef multi_index_container<
    hash_entry,
    indexed_by<
      hashed_unique<member<hash_entry,std::string,&hash_entry::hash> >
    >
> hash_database_type;

std::string make_hash(int nat, const EdgeList &edges)
{
    ostringstream ostr;
   
    Graph *graph = new Graph();
    for(int i=0; i< nat; ++i)
    {
        graph->addVertex((unsigned short)i);
    }
    for (int i=0; i< edges.size() ; ++i)
    {
        Graph::Vertex *v1,*v2;
        v1 = graph->vertices[edges[i].first];
        v2 = graph->vertices[edges[i].second];
        graph->connectVertices(v1, v2, false);
    }
    GraphHash *hash = new GraphHash ;
    hash->hash(graph, false) ;

    for (int i=0 ; i< MD5_SIZE; ++i)
    {
        ostr << boost::format("%02X") % static_cast<unsigned short>(hash->code[i]);
    }

    delete hash;
    delete graph;
    return ostr.str();
}
void add_hash_entry(hash_database_type& hash_database, int nat, const EdgeList& edges, const Polynomial& poly)
{
    hash_entry he(nat, edges, poly);
    he.hash = make_hash(nat, edges);
    std::pair<hash_database_type::iterator, bool> hdi=hash_database.insert(he);

    if ( hdi.second == false ) 
    {
        if (hdi.first->poly != he.poly)
        {
            cerr << "collision!!" << endl;
            cout << hdi.first->hash;
            cout << hdi.first->poly << endl;
            cout << he.hash << endl;
            cout << he.poly << endl;
        }
    }
}

void save_database_to_file(const hash_database_type& hash_database, const std::string& file_name)
{
    std::ofstream ofs(file_name.c_str());
    boost::archive::text_oarchive oa(ofs);
    oa << hash_database ;
}
void load_database_from_file(hash_database_type& hash_database, const std::string& file_name)
{
    std::ifstream ifs(file_name.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> hash_database ;
}

EdgeList findEdges(const GeomCoord& coords)
{
    EdgeList res;
    double dis;
    Vector3 p1, p2, diff;
    for (int i=0; i < coords.size(); ++i)
    {
        for(int j=0; j<i; ++j)
        {
            p1 = coords[i];
            p2 = coords[j];
            for(int k=0; k<3;++k) diff[k]=p1[k]-p2[k];   
            dis = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);

            if ( dis <= (CC_Tol+0.1))
            {
                res.push_back(make_pair<int,int>(i,j));
            }
        }
    }
    return res;
}

class Structure
{
public:
    GeomCoord coords;
    int nat;
    EdgeList edges;
    Structure(const std::string& xyz_filename, bool loadAll = false)
    {
        load_xyz(xyz_filename, loadAll);
        edges = findEdges(coords);
    }
    void load_xyz(const std::string& xyz_filename, bool loadAll)
    {
        std::string str, s;
        Vector3 xyz;
        stringstream ss (stringstream::in | stringstream::out);
        ifstream fin(xyz_filename.c_str());
        std::getline(fin, str);
        ss.clear();
        ss << str;
        ss >> nat;
        std::getline(fin, str);
        for (int i=0; i< nat; ++i)
        {
            std::getline(fin, str);
            ss.clear();
            ss << str;
            ss >> s >> xyz[0] >> xyz[1] >> xyz[2];
            
            if ( loadAll ) 
            {
                coords.push_back(xyz);
            }
            else if (s == "C" || s == "c")
            {
                coords.push_back(xyz);
            }
        }
        nat = coords.size() ; 
    }
};


#endif
