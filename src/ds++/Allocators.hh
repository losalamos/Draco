//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    Allocators.hh
 * \author  Geoffrey Furnish
 * \date    Mon Jan 27 12:11:50 1997
 * \brief   Memory management classes for ds++.
 * \note    Copyright (c) 1997-2010 Los Alamos National Security, LLC
 * \version $Id$
 */
//---------------------------------------------------------------------------//
#ifndef __ds_Allocators_hh__
#define __ds_Allocators_hh__

#include "Assert.hh"

//---------------------------------------------------------------------------//
/*! \note The allocators in this file are provided for use by the DS++
 * containers, or even by DS++ clients.  What is most important to understand
 * about these allocators is that they follow the STL model for allocators.
 * Specifically, they /only/ do /allocation/.  They do not concern themselves
 * with object creation per se.  When you get back memory from one of these
 * allocators, you have to initialize that memory with C++ object constructor
 * calls yourself.  Likewise, before you release memory, you must destroy the
 * objects contained in the memory pool before calling the release method on
 * the allocator.
 */
//---------------------------------------------------------------------------//

#include <memory>

/*! \note \c MSIPL_THROW is defined in KCC's memory header file. In order to
 *  allow compilation with other compilers, we define \c MSIPL_THROW to state
 *  that functions declared with \c MSIPL_THROW will not throw exceptions.
 */
#ifndef MSIPL_THROW
#define MSIPL_THROW throw()
#endif

namespace rtt_dsxx
{

// These functions were in the April '96 draft, but disappeared by December
// '96.  Too bad, b/c they are extremely helpful for building allocators.

template<class T>
inline T *ds_allocate( int size, const T * /* hint */ )
{
    return (T *) (::operator new( static_cast<size_t>(size) * sizeof(T) ) );
}

template<class T>
inline void ds_deallocate( T *buf ) { ::operator delete( buf ); }

//===========================================================================//
/*!
 * \class Simple_Allocator
 * \brief Minimalistic storage allocator
 *
 * This allocator class does the least sophisticated job possible.
 * Specifically, it allocates exactly and only the amount of memory required
 * to hold the requested number of objects of type \c T.  This is very similar
 * to new \c T[n], except that object constructors are \b not called.
 */
//===========================================================================//

template<class T>
class Simple_Allocator 
{
  public:
    typedef size_t                        size_type;
    typedef ptrdiff_t                     difference_type;
    typedef T*                            pointer;
    typedef const T*                      const_pointer;
    typedef T&                            reference;
    typedef const T&                      const_reference;
    typedef T                             value_type;
    typedef T *iterator;
    typedef const T *const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    static T *fetch( int n, const T *hint = 0 )
    {
        return ds_allocate( n, hint );
    }

    static T *validate( T *v, int /* n */ ) { return v; }

    static void release( T *v, int n =0 )
    {
        ds_deallocate( validate( v, n ) );
    }

    //#ifdef _CRAYT3E
    //size_type max_size () const MSIPL_THROW
    // Cray T3E backend sometimes incorrectly uses signed comparisons
    // instead of unsigned comparisons.  Therefore for T3E, we use
    // a max_size that is positive even if misinterpreted by backend.
    //    { return sizeof(T)==1 ? size_type( size_type (-1)/2u ) :
    //      ( size_type(1) > size_type (size_type (-1)/sizeof (T)) ) ?
    //            size_type(1) :
    //               size_type (size_type (-1)/sizeof (T)); }
    //#else
    // { return // max (size_type (1), size_type (size_type (-1)/sizeof (T))); 
    size_type max_size () const MSIPL_THROW
    { return
            ( size_type(1) > size_type (size_type (-1)/sizeof (T)) )
            ? size_type(1):size_type(size_type(-1)/sizeof (T));
    }
    //#endif /*_CRAYT3E*/

};

//===========================================================================//
/*!
 * \class Guarded_Allocator
 * \brief Range checking allocator
 *
 * This allocator class allocates enough memory to hold the requested number
 * of objects of the specified type \c T, but also allocates a little extra to
 * pad each end of the memory region.  This extra padding area is seeded with
 * a specific byte pattern on allocation, and is checked for this same pattern
 * on deallocation.  If the memory has been corrupted at the time of
 * deallocation, you have a clear indication of an erroneous program, and an
 * exception is consequently thrown.
 */
//===========================================================================//

template<class T>
class Guarded_Allocator {
  public:
    typedef size_t                        size_type;
    typedef ptrdiff_t                     difference_type;
    typedef T*                            pointer;
    typedef const T*                      const_pointer;
    typedef T&                            reference;
    typedef const T&                      const_reference;
    typedef T                             value_type;
    typedef T *iterator;
    typedef const T *const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    static T *fetch( int n, const T* hint = 0 )
    {
        Assert( n >= 0 );
        T *v = ds_allocate( n+2, hint );

        // Write magic info into T[0] and T[n+1]
        char *pb = reinterpret_cast<char *>(v);
        char *pe = reinterpret_cast<char *>( v+n+1 );
        for( size_t i=0; i < sizeof(T); i++ ) 
        {
            /*! \bug char and hex are different length data.  The following
             *  casts will issue warning C4309 'truncation of constant value'
             */
            pb[i] = static_cast<char>( 0xE8 );
            pe[i] = static_cast<char>( 0xE9 );
        }
	
        return v+1;
    }

    static T *validate( T *v, int Remember(n) )
    {
        v--;
        Check( guard_elements_ok(v,n) );
        return v;
    }

    //! \brief Check magic data in T[0] and T[n+1]
    static bool guard_elements_ok( T *v, int n )
    {
        char *pb = reinterpret_cast<char *>(v);
        char *pe = reinterpret_cast<char *>( v+n+1 );

        for( size_t i=0; i < sizeof(T); i++ )
        {
            /*! \bug char and hex are different length data.  The following
             *  casts will issue warning C4309 'truncation of constant value'
             */
            if (pb[i] != static_cast<char>( 0xE8 )) return false;
            if (pe[i] != static_cast<char>( 0xE9 )) return false;
        }

        return true;
    }

    static void release( T *v, int n )
    {
        if (!v) return;

        ds_deallocate( validate( v, n ) );
    }

    size_type max_size () const MSIPL_THROW
#ifdef _CRAYT3E
    // Cray T3E backend sometimes incorrectly uses signed comparisons
    // instead of unsigned comparisons.  Therefore for T3E, we use
    // a max_size that is positive even if misinterpreted by backend.
    { return sizeof(T)==1 ? size_type( size_type (-1)/2u ) :
            ( size_type(1) > size_type (size_type (-1)/sizeof (T)) ) ? size_type(1) :
            size_type (size_type (-1)/sizeof (T)); }
#else
    { return // max (size_type (1), size_type (size_type (-1)/sizeof (T))); 
            ( size_type(1) > size_type (size_type (-1)/sizeof (T)) ) ? size_type(1):size_type(size_type(-1)/sizeof (T));}
#endif /*_CRAYT3E*/

};

//===========================================================================//
/*!
 * \class alloc_traits
 * \brief Traits class for specifying client allocators
 *
 * This allocator traits class is provided to allow traits specialization of
 * the default allocator to be used by ds++ container classes for various user
 * defined types.  The default ds++ allocator will be the Guarded_Allocator
 * (eventually, once I am sure it is working right), but individual user
 * defined types can override this default on a case by case basis using
 * template specialization of this class.  For example, to change the default
 * allocator for wombats, you could put the following at the top of wombat.hh:
 *
 * \code
 * #include "Allocators.hh"
 * template<> class alloc_traits<wombat> {
 *   public:
 *     typedef Simple_Allocator<wombat> Default_Allocator;
 * }
 * \endcode
 *
 * This way, when a client code does this:
 *
 * \code 
 * #include "wombat.hh"
 * #include "Mat.hh"
 *
 * Mat1<wombat> wombatmat;
 * \endcode
 *
 * The \c Simple_Allocator will be used rather than the \c Guarded_Allocator. 
 */
//===========================================================================//

//! \bug alloc_traits is not being tested.
//lint -e758 Ignore warning about this function not being referenced until
//           the unit test has been beefed up.
//lint -e526 Not sure about this warning "Simple_Allocator<<1>> not defined.
//lint -e768 Default_Allocator not referenced

template<class T>
class alloc_traits
{
  public:
    typedef Simple_Allocator<T> Default_Allocator;
};

template<> class alloc_traits<int>
{
  public:
    typedef Guarded_Allocator<int> Default_Allocator;
};

} // end of rtt_dsxx

#endif // __ds_Allocators_hh__

//---------------------------------------------------------------------------//
// end of ds++/Allocators.hh
//---------------------------------------------------------------------------//
