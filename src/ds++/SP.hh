//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/SP.hh
 * \author Geoffrey Furnish, Thomas Evans
 * \date   Tue Feb  4 11:27:53 2003
 * \brief  Smart Point class file.
 * \note   Copyright (C) 2003-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_ds_SP_HH
#define RTT_ds_SP_HH

#include <typeinfo>
#include "ds++/config.h"
#include "Assert.hh"

namespace rtt_dsxx
{
 
//===========================================================================//
/*!
 * \struct SPref
 * 
 * \brief Reference holder struct for SP class.
 */
//===========================================================================//

struct SPref 
{
    //! Number of references.
    int refs;

    //! Constructor
    SPref(int r = 1) : refs(r) { }

#if DRACO_DIAGNOSTICS & 2
    // To prevent refcounts from being included in memory debugging.
    // These are problematic because we allocate reference counters even
    // for null smart pointers.
    void *operator new(size_t const n)
    {
        return malloc(n);
    }
    void  operator delete(void *const ptr)
    {
        free(ptr);
    }
#endif
};

//===========================================================================//
/*!
 * \class SP
 * 
 * \brief Smart pointer implementation that does reference counting.
 *
 * The smart pointer provides a "safe" encapsulation for a standard C++
 * pointer.  Consider: A function new's an object and return the pointer to
 * that object as its return value.  Now it is the caller's responsibility to
 * free the object.  What if the caller passes the pointer to other objects
 * or functions?  What if it is not known which will be deleted first or
 * last?
 *
 * Instead the function can return a "smart pointer".  This SP class uses
 * reference counting to determine the number of current users of a pointer.
 * Each time an SP goes out of scope, the reference count is decremented.
 * When the last user of a pointer is done, the pointer is freed.
 *
 * Note: I am calling this a "smart pointer", not a "safe pointer".  There
 * are clearly ways you can hose this.  In particular, when you bind an SP<T>
 * to a T*, you yield all rights to the T*.  You'd better not squirrel the
 * bare pointer away somewhere and expect to clandestinely use it in other
 * ways or places--death will be sure to follow.  Consequently then, the
 * safest way to use this smart pointer, is to bind it to the contained
 * pointer and then always use the smart pointer.  Immediately returning the
 * smart pointer as a return value, allowing the original bare pointer to go
 * out of scope never to be seen again, is one good example of how to use
 * this.
 * 
 * One good example of bad usage is assigning the same dumb pointer to
 * multiple SPs.  Consider:
 * \code
 *     SP<Foo> f1;
 *     SP<Foo> f2;
 *     // ...
 *     Foo *f = new Foo;
 *     f1 = f;
 *     // ...
 *     f2 = f; // bad, now f1 and f2 assume they "own" f!
 * \endcode
 * Unfortunately, there is no way to check if another SP owns the dumb
 * pointer that you give to a SP.  This is simply something that needs to be
 * watched by the programmer.
 * 
 * \note
 * Having an std::vector or other array-based container of SPs can have
 * non-obvious implications for object lifetime. Since operations like
 * pop_back() or clear() do not call the destructor of the SP, SP's that are
 * in the "slop" between the vector size and the capacity are still holding
 * references.  It is only once those slots are re-assigned, or the whole
 * vector is deleted, or a resize operation leads to a reallocation, that
 * those SPs are killed.
 */
/*!
 * \example ds++/test/tstSP.cc
 *
 * rtt_dsxx::SP (smart pointer) usage example.
 */
// revision history:
// -----------------
// 0) original
// 1) 020403 : updated with doxygen comments; minor refactoring
// 
//===========================================================================//

template<typename T>
class SP 
{
  private: 
    // >>> DATA

    //! Raw pointer held by smart pointer.
    T *p;

    //! Pointer to reference counter.
    SPref *r;

  private:
    // >>> IMPLEMENTATION

    // Free the pointer.
    inline void free();

    //! All derivatives of SP are friends. 
    template<class X> friend class SP;

  public:
    //! Default constructor.
    SP() : p(NULL), r(new SPref) { Ensure (r); Ensure (r->refs == 1); }

    // Explicit constructor for type T *.
    inline explicit SP(T *p_in);

    // Explicit constructor for type X *.
    template<class X>
    inline explicit SP(X *px_in);

    // Copy constructor for SP<T>.
    inline SP(const SP<T> &sp_in);

    // Copy constructor for SP<X>.
    template<class X>
    inline SP(const SP<X> &spx_in);

    //! Destructor, memory is released when count goes to zero.
    ~SP(void) { free(); }

    // Assignment operator for type T *.
    inline SP<T>& operator=(T *p_in);

    // Assignment operator for type X *.
    template<class X>
    inline SP<T>& operator=(X *px_in);

    // Assignment operator for type SP<T>.
    inline SP<T>& operator=(const SP<T> sp_in);

    // Assignment operator for type SP<X>.
    template<class X>
    inline SP<T>& operator=(const SP<X> spx_in);

    //! Access operator.
    T* operator->() const { Require(p); return p; }

    //! Dereference operator.
    T& operator*() const { Require(p); return *p; }

    //! Get the base-class pointer; better know what you are doing.
    T* bp() const { return p; }

    //! Boolean conversion operator.
    operator bool() const { return p != NULL; }

    //! Operator not.
    bool operator!() const { return p == NULL; }

    //! Equality operator for T*.
    bool operator==(const T *p_in) const { return p == p_in; }

    //! Inequality operator for T*.
    bool operator!=(const T *p_in) const { return p != p_in; }

    //! Equality operator for SP<T>.
    bool operator==(const SP<T> &sp_in) const { return p == sp_in.p; }

    //! Inequality operator for SP<T>.
    bool operator!=(const SP<T> &sp_in) const { return p != sp_in.p; }
};

DLL_PUBLIC void incompatible(std::type_info const &X, std::type_info const &T);

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Do equality check with a free pointer.
 */
template<typename T>
bool operator==(const T *pt, const SP<T> &sp)
{
    return sp == pt;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do inequality check with a free pointer.
 */
template<typename T>
bool operator!=(const T *pt, const SP<T> &sp)
{
    return sp != pt;
}

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Explicit constructor for type T *.
 *
 * This constructor is used to initialize a SP with a pointer, ie.
 * \code
 *     Foo *f = new Foo;
 *     SP<Foo> spf(f);   // f now owned by spf
 *     SP<Foo> spf2 = f; // error! does not do implicit conversion
 * \endcode
 * Once a pointer is "given" to a SP, the SP takes control. This means that
 * spf will delete f when the last SP to f is destroyed.
 *
 * \param p_in pointer to type T
 */
template<typename T>
SP<T>::SP(T *p_in)
    : p(p_in),
      r(new SPref)
{
    Ensure (r);
    Ensure (r->refs == 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Explicit constructor for type X *.
 *
 * This constructor is used to initialize a base class smart pointer of type
 * T with a derived class pointer of type X or, equivalently, any types in
 * which X * is convertible to T * through a dynamic_cast.  Consider,
 * \code
 *     class Base {//...}; 
 *     class Derived : public Base {//...};
 *     
 *     SP<Base> spb(new Derived); // spb base class SP to Derived type
 *     Derived *d = new Derived;
 *     SP<Base> spb2(d);          // different syntax
 *     SP<Base> spb3 = d;         // error! no implicit conversion
 *     
 *     Derived *d2;
 *     SP<Base> spb4(d2);         // error! cannot initialize with NULL
 *                                // pointer of different type than T
 * \endcode
 * The pointer to X must not be equal to NULL.  The SP owns the pointer when
 * it is constructed.
 *
 * \param px_in pointer to type X that is convertible to T *
 */
template<typename T>
template<class X>
SP<T>::SP(X *px_in)
    : p(NULL), r(NULL)
{
    Require (px_in);

    // make a dynamic cast to check that we can cast between X and T
    T *np = dynamic_cast<T *>(px_in);

    // check that we have made a successfull cast if px exists
    if(!np)
        incompatible(typeid(X), typeid(T));

    // assign the pointer and reference
    p = np;
    r = new SPref;

    Ensure (p);
    Ensure (r);
    Ensure (r->refs == 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for SP<T>.
 *
 * \param sp_in smart pointer of type SP<T>
 */
template<typename T>
SP<T>::SP(const SP<T> &sp_in)
    : p(sp_in.p),
      r(sp_in.r)
{
    Require (r);

    // advance the reference to T
    r->refs++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for SP<X>.
 *
 * This copy constructor requires that X * is convertible to T * through a
 * dynamic_cast.  The pointer in spx_in can point to NULL; however, it must
 * still be convertible to T * through a dynamic_cast.
 *
 * \param spx_in smart pointer of type SP<X>
 */
template<typename T>
template<class X>
SP<T>::SP(const SP<X> &spx_in)
    :p(NULL), r(NULL)
{
    Require (spx_in.r);

    // make a pointer to T *
    T *np = dynamic_cast<T *>(spx_in.p);
    if (spx_in.p ? np == 0 : false)
        incompatible(typeid(X), typeid(T));

    // assign the pointer and reference
    p = np;
    r = spx_in.r;
    
    // advance the reference to T
    r->refs++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type T *.
 *
 * The assignment operator checks the existing count of the smart pointer and
 * then assigns the smart pointer to the raw pointer.  As in copy
 * construction, the smart pointer owns the pointer after assignment. Here is
 * an example of usage:
 * \code
 *     SP<Foo> f;         // has reference to NULL pointer
 *     f      = new Foo;  // now has 1 count of Foo
 *     Foo *g = new Foo; 
 *     f      = g;        // f's original pointer to Foo is deleted
 *                        // because count goes to zero; f now has
 *                        // 1 reference to g
 * \endcode
 * 
 * \param p_in pointer to T
 */
template<typename T>
SP<T>& SP<T>::operator=(T *p_in)
{
    // check if we already own this pointer
    if (p == p_in)
	return *this;

    // next free the existing pointer
    free();
    
    // now make add p_in to this pointer and make a new reference to it
    p = p_in;
    r = new SPref;

    Ensure (r->refs == 1);
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type X *.
 *
 * This assignment requires that X * is convertible to T * through a
 * dynamic_cast.  It follows the same principle as SP(X*); however,
 * this is assignment:
 * \code
 *     SP<Base> b;
 *     b = new Derived;
 * \endcode
 * The pointer to X must not be equal to NULL.
 *
 * \param px_in pointer to type X * that is convertible to type T * through a
 * dynamic cast
 */
template<typename T>
template<class X>
SP<T>& SP<T>::operator=(X *px_in)
{
    Require (px_in);

    // do a dynamic cast to ensure convertiblility between T* and X*
    T *np = dynamic_cast<T *>(px_in);
    if (!np)
        incompatible(typeid(X), typeid(T));

    // now assign this to np (using previously defined assignment operator)
    *this = np;
    
    Ensure (r->refs == 1);
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type SP<T>.
 *
 * \param sp_in smart pointer of type SP<T>
 */
template<typename T>
SP<T>& SP<T>::operator=(const SP<T> sp_in)
{
    Require (sp_in.r);

    // see if they are equal
    if (this == &sp_in || p == sp_in.p)
	return *this;

    // free the existing pointer
    free();

    // assign p and r to sp_in
    p = sp_in.p;
    r = sp_in.r;

    // add the reference count and return
    r->refs++;
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type SP<X>.
 *
 * This assignment requires that X * is convertible to T * through a
 * dynamic_cast.  The pointer in spx_in can point to NULL; however, it must
 * still be convertible to T * through a dynamic_cast.
 *
 * \param spx_in smart pointer of type SP<X>
 */
template<typename T>
template<class X>
SP<T>& SP<T>::operator=(const SP<X> spx_in)
{
    Require (spx_in.r);

    // make a pointer to T *
    T *np = dynamic_cast<T *>(spx_in.p);
    if (spx_in.p ? np == 0 : false)
        incompatible(typeid(X), typeid(T));

    // check to see if we are holding the same pointer (and np is not NULL);
    // to NULL pointers to the same type are defined to be equal by the
    // standard 
    if (p == np && p)
    {
	// if the pointers are the same the reference count better be the
	// same 
	Check (r == spx_in.r);
	return *this;
    }

    // we don't need to worry about the case where p == np and np == NULL
    // because this combination is impossible; if np is NULL then it belongs
    // to a different smart pointer; in other words, if p == np and np ==
    // NULL then r != spx_in.r

    // free the existing pointer and reference
    free();

    // assign new values
    p = np;
    r = spx_in.r;

    // advance the counter and return
    r->refs++;
    return *this;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Decrement the count and free the pointer if count is zero.
 *
 * Note that it is perfectly acceptable to call delete on a NULL pointer.
 */
template<typename T>
void SP<T>::free()
{
    Require (r);
    
    // if the count goes to zero then we free the data
    if (--r->refs == 0)
    {
	delete p;
	delete r;
    }
}

} // end namespace rtt_dsxx

#endif                          // RTT_ds_SP_HH

//---------------------------------------------------------------------------//
// end of ds++/SP.hh
//---------------------------------------------------------------------------//
