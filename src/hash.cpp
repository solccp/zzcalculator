
#include <iostream>
#include "stdint.h"
#include <cstring>
#include <cstdio>
#include "graphHash.hpp"
#include <sstream>
#include <cstdlib>

#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

// a non-portable native binary archive
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

using namespace std;



class Term
{
public:
    friend class boost::serialization::access;
    std::vector<int64_t> m_coeffs;
    int64_t m_leadpow;
    int m_block_size;
    int m_order;
    Term()
    {
        m_order = 0;
        m_block_size = 0;
        m_leadpow = 0;
    }
    Term(int order, int block_size, int64_t leadpow, int64_t *coeffs)
    {
        this->m_order = order ;
        this->m_block_size = block_size;
        this->m_leadpow = leadpow;
        for (int i=0; i<leadpow;++i)
        {
            m_coeffs.push_back(coeffs[i]);
        }
    }
    friend std::ostream& operator<<(std::ostream&, Term&);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int file_version)
    {
        ar & m_order & m_block_size & m_leadpow & m_coeffs ; 
    }

};

std::ostream& operator<<(std::ostream& os, Term& term)
{
    os << "[" ; 
    for (int i = term.m_leadpow-1; i>=0; --i)
    {
        os << term.m_coeffs[i]; 
        if (i>0) os << "," ;
    }
    os << "]" ;
    if ( term.m_order >= 2 )
    {
        os << "x^" << term.m_order ;
    }
    else if (term.m_order == 1)
    {
        os << "x";
    }
    return os;
}
class Polynomial
{
public:
    friend class boost::serialization::access;
    std::vector<Term> m_terms;
    Polynomial(const std::vector<Term>& terms): m_terms(terms)
    {}
    Polynomial()
    {
    }
    bool operator!=(Polynomial& rhs)
    {
        stringstream ss1(stringstream::in | stringstream::out);
        stringstream ss2(stringstream::in | stringstream::out);

        ss1 << *this;
        ss2 << rhs;

        return (ss1.str() != ss2.str()) ;

    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int file_version)
    {
        ar & m_terms; 
    }
    friend std::ostream& operator<<(std::ostream&, Polynomial&);
};

std::ostream& operator<<(std::ostream& os, Polynomial& poly)
{
    for(int i=0; i< poly.m_terms.size(); ++i)
    {
        if (i>0)
        {
            os << "+ " ;
        }
        os << poly.m_terms[i] << " " ;
    }
    return os;
}


static std::map<std::string, Polynomial> hash_database;

extern "C"
{
    void make_hash_(int32_t* nat, int32_t* num_of_edges, int32_t* edges, char* key)
    {
        char tmp_key[3];
//        int nthreads = 1;
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
        GraphHash *hash = new GraphHash() ; //(graph, false);
//        char* ncpu = getenv ( "OMP_NUM_THREADS" );
//        if ( ncpu != NULL )
//        {
//            int tmp_ncpus = atoi ( ncpu );
//            if (tmp_ncpus != 0)
//                nthreads = tmp_ncpus;
//        }

        hash->hash(graph, false);

        for (int i=0 ; i< MD5_SIZE; ++i)
        {
            sprintf(tmp_key, "%02X", hash->code[i]);
            key[i*2]=tmp_key[0];
            key[i*2+1]=tmp_key[1];
        }

        delete hash;
        delete graph;
    }

    void get_polynomial_order_kernel_(char key[32], int32_t* order)
    {
        std::string hash_key;
        hash_key.assign(key,32);
        std::map<std::string, Polynomial>::iterator p = hash_database.find(hash_key);
        if ( p!= hash_database.end() )
        {
            *order = p->second.m_terms.size()-1 ;
        }
        else
        {
            *order = -1 ;
        }
    }

    void get_polynomial_kernel_(char key[32], int32_t* block_size, int64_t* poly)
    {
        std::string hash_key;
        hash_key.assign(key,32);
        std::map<std::string, Polynomial>::iterator p = hash_database.find(hash_key);
        if ( p!= hash_database.end() )
        {
            for(int i=0; i< p->second.m_terms.size(); ++i)
            {
                int offset = (*block_size)+1 ;
                poly[i*offset] = p->second.m_terms[i].m_leadpow ;
                for(int j=0; j<p->second.m_terms[i].m_leadpow; ++j)
                {
                    poly[i*offset+j+1] = p->second.m_terms[i].m_coeffs[j] ;
                }
            }
        }
        else
        {
            throw "some logic error";
        }
    }

    void add_polynomial_kernel_(char key[32], int64_t *nterms, int32_t *block_size, int64_t *poly)
    {
        std::string hash_key;
        hash_key.assign(key,32);
//        cout << hash_key.size() << endl;
//        cout << "add poly" << endl;
//        cout << "hash: " << hash_key << endl;
//        cout << "terms: " << *nterms << endl;

        std::vector<Term> terms;
        for(int i=0; i< *nterms; ++i)
        {
            int offset = (*block_size)+1 ;
            Term t(i, *block_size, poly[i*offset], &poly[i*offset+1]);
            if ( t.m_leadpow == 1 && t.m_coeffs[0] ==0 ) 
                continue;
            terms.push_back(t);
        }
        std::pair<std::map<std::string, Polynomial>::iterator,bool> ret;
        std::pair<std::string, Polynomial> new_entry(hash_key, Polynomial(terms));
        ret = hash_database.insert(new_entry);
/*        if ( ret.second == false ) 
        {
            if ( new_entry.second != ret.first->second )
            {
                cout << "same hash" << endl;
                cout << "two str have the same HASH!!!" << endl;
                cout << hash_key << endl;
                cout << new_entry.first << endl;
                cout << new_entry.second << endl;
                cout << ret.first->second << endl;
            }
        }
*/
    }

    void print_all_database_entry_()
    {
        std::map<std::string, Polynomial>::iterator p = hash_database.begin();
        for ( ;p!= hash_database.end(); ++p)
        {
            cout << (*p).first << "  " << (*p).second << endl;
        }
    }

    void save_database_to_file_() //(int len, char* filename)
    {
        cout << "database size: " << hash_database.size() << endl;
        std::ofstream ofs("database.new");
        {
            boost::archive::binary_oarchive oa(ofs);
//            boost::archive::text_oarchive oa(ofs);
            oa << hash_database ;
        }
    }

    void load_database_from_file_() // (int len, char* filename)
    {
        std::ifstream ifs("database");
        {
            boost::archive::binary_iarchive ia(ifs);
//            boost::archive::text_iarchive ia(ifs);
            ia >> hash_database ;
        }
        cout << "Loaded database size: " << hash_database.size() << endl;
    }


};
