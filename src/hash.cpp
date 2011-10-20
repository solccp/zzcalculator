
#include "polynomial.hpp"
#include "hash_database.hpp"
#include <cstring>

using namespace std;

static hash_database_type hash_database;


extern "C" 
{
    void make_hash_(int32_t* nat, int32_t* num_of_edges, int32_t* edges, char* key)
    {

        EdgeList edgesL;
        for (int i=0; i< *num_of_edges ; ++i)
        {
            edgesL.push_back(make_pair(edges[i*2]-1, edges[i*2+1]-1));
        }
        std::string strkey = make_hash(*nat, edgesL); 
        memcpy(key, &strkey[0], 32);
    }

    void add_polynomial_kernel_new_(int32_t* nat, int32_t* num_of_edges, int32_t* edges, int64_t *nterms, int32_t *block_size, int64_t *poly)
    {
        EdgeList edgelist;

        for (int i=0; i< *num_of_edges ; ++i)
        {
            edgelist.push_back(make_pair(edges[i*2]-1, edges[i*2+1]-1));
        }

        std::vector<Term> terms;
        for(int i=0; i< *nterms; ++i)
        {
            int offset = (*block_size)+1 ;
            Term t(i, *block_size, poly[i*offset], &poly[i*offset+1]);
            if ( t.m_leadpow == 1 && t.m_coeffs[0] ==0 ) 
                continue;
            terms.push_back(t);
        }
        add_hash_entry(hash_database, *nat, edgelist, Polynomial(terms));
    }


    void save_database_to_file_() //(int len, char* filename)
    {
        save_database_to_file(hash_database, "database.new");
    }

    void load_database_from_file_() // (int len, char* filename)
    {
        load_database_from_file(hash_database, "database");
        cerr << "Loaded database size: " << hash_database.size() << endl;
    }

    void print_database_()
    {
        cout << "database size: " << hash_database.size() << endl;
        for(hash_database_type::iterator wit=hash_database.begin(),wit_end=hash_database.end(); wit!=wit_end;++wit)
        {
            std::cout<< wit->hash << " " << wit->natoms << " " << wit->poly << std::endl;
        }


    }
    void get_polynomial_order_kernel_(char key[32], int32_t* order)
    {
        std::string hash_key;
        hash_key.assign(key,32);
        hash_database_type::iterator p = hash_database.find(hash_key);
        if ( p!= hash_database.end())
        {
            *order = p->poly.m_terms.size()-1 ;
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
        hash_database_type::iterator p = hash_database.find(hash_key);
        if ( p!= hash_database.end() )
        {
            for(int i=0; i< p->poly.m_terms.size(); ++i)
            {
                int offset = (*block_size)+1 ;
                poly[i*offset] = p->poly.m_terms[i].m_leadpow ;
                for(int j=0; j<p->poly.m_terms[i].m_leadpow; ++j)
                {
                    poly[i*offset+j+1] = p->poly.m_terms[i].m_coeffs[j] ;
                }
            }
        }
        else
        {
            throw "some logic error";
        }
    }


};
