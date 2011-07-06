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
#include  <limits>

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
 *
 * [2011/07/06 KT] I am updating this class to make it more standardized to
 * modern design for Allocator classes.  See
 * http://www.codeproject.comm/KB/cpp/allocator.aspx. This article points to
 * additional resources:
 * - Nicolai M. Josuttis, "The C++ Standard Library: A Tutorial and Reference"
 * - ISO C++ Standard, 1998.
 * - Stanley B. Lippman, "Inside the C++ Object Model"
 * - Andrei Alexandrescu, "Modern C++ Design: Generic Programming and Design
 *   Patterns Applied"
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

// template<class T>
// inline T *ds_allocate( int size, const T * /* hint */ )
// {
//     return (T *) (::operator new( static_cast<size_t>(size) * sizeof(T) ) );
// }

// template<class T>
// inline void ds_deallocate( T *buf ) { ::operator delete( buf ); }

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

template<typename T>
class Simple_Allocator 
{
  public:

    // TYPEDEFS
    
    typedef       size_t                  size_type;
    typedef       ptrdiff_t               difference_type;
    typedef       T*                      pointer;
    typedef const T*                      const_pointer;
    typedef       T&                      reference;
    typedef const T&                      const_reference;
    typedef       T                       value_type;
    typedef       T*                      iterator;
    typedef const T*                      const_iterator;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    // CONSTRUCTORS

    inline Simple_Allocator( void) {/*empty*/}
    inline ~Simple_Allocator(void) {/*empty*/}
    inline explicit Simple_Allocator(Simple_Allocator const & ) {/*empty*/}

    // ADDRESS

    inline pointer       address(reference       r) { return &r; }
    inline const_pointer address(const_reference r) { return &r; }

    // MEMORY ALLOCATION
    
    static pointer fetch(
        size_type n, typename std::allocator<void>::const_pointer hint = 0 )
    {
        return allocate(n,hint);//ds_allocate( n, hint );
    }
    static pointer allocate(
        size_type n, typename std::allocator<void>::const_pointer = 0 )
    {
        return reinterpret_cast<pointer>(::operator new(n * sizeof(T) ));
    }
    static pointer validate( T *v, int /* n */ ) { return v; }
    static void release( T *v, int n =0 )
    {
        deallocate(v,n);
        return;
        // ds_deallocate( validate( v, n ) );
    }
    static void deallocate(pointer p, size_type)
    {
        ::operator delete(p);
        return;
    }

    //  SIZE

    inline size_type max_size() const
    {
        return std::numeric_limits<size_type>::max() / sizeof(T);
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
    // inline size_type max_size () const MSIPL_THROW
    // { return
    //         ( size_type(1) > size_type (size_type (-1)/sizeof (T)) )
    //         ? size_type(1):size_type(size_type(-1)/sizeof (T));
    // }
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

template<typename T>
class Guarded_Allocator
{
  public:

    // TYPEDEFS
    
    typedef       size_t                  size_type;
    typedef       ptrdiff_t               difference_type;
    typedef       T*                      pointer;
    typedef const T*                      const_pointer;
    typedef       T&                      reference;
    typedef const T&                      const_reference;
    typedef       T                       value_type;
    typedef       T*                      iterator;
    typedef const T*                      const_iterator;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    // CONSTRUCTORS

    inline Guarded_Allocator( void) {/*empty*/}
    inline ~Guarded_Allocator(void) {/*empty*/}
    inline explicit Guarded_Allocator(Guarded_Allocator const & ) {/*empty*/}

    // ADDRESS

    inline pointer       address(reference       r) { return &r; }
    inline const_pointer address(const_reference r) { return &r; }

    // MEMORY ALLOCATION
    
    //static pointer fetch( int n, const T* hint = 0 )
    static pointer fetch(
        size_type n, typename std::allocator<void>::const_pointer hint = 0 )
    {
        return allocate(n,hint);//ds_allocate( n, hint );
    }
    static pointer allocate(
        size_type n, typename std::allocator<void>::const_pointer = 0 )
    {
        // Assert( n >= 0 );
        // pointer v = ds_allocate( n+2, hint );
        pointer v = reinterpret_cast<pointer>(::operator new((n+2) * sizeof(T) ));

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

    static pointer validate( pointer v, size_type Remember(n) )
    {
        v--;
        Check( guard_elements_ok(v,n) );
        return v;
    }

    //! \brief Check magic data in T[0] and T[n+1]
    static bool guard_elements_ok( pointer v, size_type n )
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

    static void release( pointer v, size_type n )
    {
        deallocate(v,n);
        return;
    }

    static void deallocate( pointer v, size_type n ) 
    {
        if (!v) return;
        validate( v, n );
        ::operator delete(v-1);
        return;
    }
    inline size_type max_size () const // MSIPL_THROW
    {
        // subtract 2 from max() to treat guarded values at begining and end.
        return (std::numeric_limits<size_type>::max() / sizeof(T) ) - 2;
    }
// #ifdef _CRAYT3E
//     // Cray T3E backend sometimes incorrectly uses signed comparisons
//     // instead of unsigned comparisons.  Therefore for T3E, we use
//     // a max_size that is positive even if misinterpreted by backend.
//     { return sizeof(T)==1 ? size_type( size_type (-1)/2u ) :
//             ( size_type(1) > size_type (size_type (-1)/sizeof (T)) ) ? size_type(1) :
//             size_type (size_type (-1)/sizeof (T)); }
// #else
//     { return // max (size_type (1), size_type (size_type (-1)/sizeof (T))); 
//             ( size_type(1) > size_type (size_type (-1)/sizeof (T)) ) ? size_type(1):size_type(size_type(-1)/sizeof (T));}
// #endif /*_CRAYT3E*/

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

template<typename T>
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
