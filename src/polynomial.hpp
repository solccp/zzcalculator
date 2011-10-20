
#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <stdint.h>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <cstdlib>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
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
    {}
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
    friend std::ostream& operator<<(std::ostream&, const Term&);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int file_version)
    {
        ar & m_order & m_block_size & m_leadpow & m_coeffs ; 
    }

};

std::ostream& operator<<(std::ostream& os, const Term& term)
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
    bool operator!=(const Polynomial& rhs) const
    {
        return !(*this == rhs) ;
    }
    bool operator==( const Polynomial& rhs) const
    {
        stringstream ss1(stringstream::in | stringstream::out);
        stringstream ss2(stringstream::in | stringstream::out);

        ss1 << *this;
        ss2 << rhs;

        return (ss1.str() == ss2.str()) ;

    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int file_version)
    {
        ar & m_terms; 
    }
    friend std::ostream& operator<<(std::ostream&, const Polynomial&);
};

    std::ostream& operator<<(std::ostream& os, const Polynomial& poly)
    {
        if ( poly.m_terms.empty() )
        {
            os << "0" ;
        }   
        else
        {
            for(int i=0; i< poly.m_terms.size(); ++i)
            {
                if (i>0)
                {
                    os << "+ " ;
                }
                os << poly.m_terms[i] << " " ;
            }
        }
        return os;
    }



#endif
