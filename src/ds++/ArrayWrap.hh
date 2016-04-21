//----------------------------------*-C++-*----------------------------------//
/*!
  \file    ds++/ArrayWrap.hh
  \author  Paul Henning
  \brief   Declaration of class ArrayWrap
  \note    Copyright 2016 Los Alamos National Security, LLC.
  \version $Id$
*/
//---------------------------------------------------------------------------//
#ifndef ArrayWrap_hh
#define ArrayWrap_hh

#include "Assert.hh"
#include "ds++/config.h"

#ifdef HAS_CXX11_ARRAY

// This file is deprecated by std::array.  However, as of Oct, 2013, PGI 13.7
// and XLC 12.1 still did not provide std::array I am keeping this class for
// those two compilers.

// This class should not be removed unless
// jayenne/clubimc/src/imc/*_Particle.* routines can be compiled with
// std::array.  

#else

namespace rtt_dsxx
{

//==============================================================================
/*!
 * \class ArrayWray
 * This is sort of like the TR1 array<> class, but it uses our range checking
 * mechanism and doesn't introduce other iterator types. The reason for using
 * this class is to allow zero-length arrays, and to provide normal container
 * iterator stuff.
 */
template<class T, unsigned N>
class ArrayWrap
{
  public:
    typedef T& reference;
    typedef T const & const_reference;
    typedef T* iterator;
    typedef T const * const_iterator;
    typedef unsigned size_type;
    typedef T value_type;

  public:
    
    bool empty() const { return false; }

    size_type size() const { return N; }

    reference operator[](unsigned const i)
    {
        Require(i < N);
        return d_data[i];
    }

    const_reference operator[](unsigned const i) const
    {
        Require(i < N);
        return d_data[i];
    }

    iterator begin() { return d_data; }
    iterator end() { return d_data + N; }
    const_iterator begin() const { return d_data; }
    const_iterator end() const { return d_data + N; }

    reference front() { return d_data[0]; }
    reference back() { return d_data[N-1]; }
    
    T* c_array() { return d_data; }
    T const * c_array() const { return d_data; }

    T const * operator+(const unsigned n) const
    {
        Require(n < N);
        return d_data + n;
    }

    T * operator+(const unsigned n) 
    {
        Require(n < N);
        return d_data + n;
    }

  private:
    T d_data[N];
};


// ---------------------------------------------------------------------------


/* Specialization of the class for zero-length arrays */
template<class T> class ArrayWrap<T, 0>
{
  public:
    typedef T&  reference;
    typedef T const & const_reference;
    typedef T* iterator;
    typedef T const * const_iterator;
    typedef unsigned size_type;

  public:

    bool empty() const { return true; }
    size_type size() const { return 0; }

    iterator begin() { return iterator(0); }
    iterator end() { return iterator(0); }
    const_iterator begin() const { return const_iterator(0); }
    const_iterator end() const { return const_iterator(0); }

    reference operator[](const unsigned) 
    {
        Require(0);
        return d_data;
    }

    const_reference operator[](const unsigned) const
    {
        Require(0);
        return d_data;
    }


    T* c_array() { return 0; }
    T const * c_array() const { return 0; }

  private:
    T d_data;

};

}

#endif // HAS_CXX11_ARRAY

#endif
