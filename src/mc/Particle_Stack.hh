//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Particle_Stack.hh
 * \author Thomas M. Evans
 * \date   Fri Dec 21 10:12:18 2001
 * \brief  Particle_Stack.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Particle_Stack_hh__
#define __mc_Particle_Stack_hh__

#include "ds++/SP.hh"
#include <vector>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Particle_Stack
 *
 * \brief Stack class for holding particle types (PT).
 *
 * The Particle_Stack class is an implementation of a STL-like stack class.
 * That is, storage is last-in-first-out.  The only difference is that it is
 * defined on std::vector, and it provides and overloaded [] operator for
 * subscripting access.
 *
 * This was originally implemented because of deficiencies in the KCC
 * implementation of std::stack; however, we found that we needed the extra
 * subscripting functionality.
 *
 */
//===========================================================================//

template<class PT>
class Particle_Stack
{
  public:
    // Typedefs.
    typedef typename std::vector<PT>::value_type value_type;
    typedef typename std::vector<PT>::size_type  size_type;

    // Particle banks and containers.
    typedef Particle_Stack<rtt_dsxx::SP<PT> > Census;
    typedef Particle_Stack<rtt_dsxx::SP<PT> > Bank;

  private:
    // Container holding data in the stack.
    std::vector<PT> c;

  public:
    //! Constructor.
    explicit Particle_Stack(const std::vector<PT> &ct = std::vector<PT>()) 
	: c(ct) {}
    
    //! Query if the stack is empty.
    bool empty() const { return c.empty(); }
    
    //! Return the size of the stack.
    size_type size() const { return c.size(); }
    
    //! Get a reference to the object on the top of the stack.
    value_type& top() { return c.back(); } 

    //! Get a const reference to the object on the top of the stack.
    const value_type& top() const { return c.back(); }
    
    //! Push an object onto the top stack.
    void push(const value_type &x) { c.push_back(x); }

    //! Remove an object from the top of the stack.
    void pop() { c.pop_back(); }

    //! Overloaded operator [] for viewing elements sequentially.
    const value_type& operator[](int i) const { return c[i]; }
};

} // end namespace rtt_mc

#endif                          // __mc_Particle_Stack_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Particle_Stack.hh
//---------------------------------------------------------------------------//
