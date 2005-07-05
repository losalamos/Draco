//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Safe_Ptr.hh
 * \author Paul Henning
 * \brief  Safe pointer-like class for scalars.
 * \note   Copyright &copy; 2005 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_Safe_Ptr_hh
#define rtt_dsxx_Safe_Ptr_hh

#include "Assert.hh"
#include <typeinfo>
#include <iostream>

namespace rtt_dsxx
{

//! A reference counting object for Safe_Ptr
struct Safe_Ptr_Ref
{
    //! Constructor
    explicit Safe_Ptr_Ref(const int r = 1) 
	: refs(r)
	, deleted(false) { }

    //! Number of references.  
    /*!  
      This is signed so that we can check for over-released references. You
      could replace this with a vector of pointers back to the Safe_Ptrs if
      you wanted even more information.
    */
    int refs;

    //! If delete has been called since the last assignment
    bool deleted;
};


//! Specialized Assert.hh-like test
/*! 
  If a Safe_Ptr member function throws an exception, and the stack unwind
  to handle that exception causes ~Safe_Ptr to be called, the dtor might
  throw an exception as well, because the internal state is hosed.  In this
  replacement for Require/Check/Ensure, we mark the internal state of
  Safe_Ptr as bad.  The dtor becomes trivial in this case, thus preventing an
  exception cascade.

  We don't provide an Insist replacement, since this stuff only exists when
  any DBC level is active.

  Of course, this macro only makes sense in the definition of Safe_Ptr, since
  it assumes that there is a variable called d_bad_state available to set.
*/
#define DBC_Check(COND) \
    if (!(COND)) { \
	d_bad_state = true; \
	rtt_dsxx::toss_cookies_ptr( #COND, __FILE__, __LINE__ ); \
    } 

#define DBC_Check_Eq(VAR, VAL)			\
    if (VAR != VAL) { \
	std::cout << "Expecting " #VAR " == " << VAL << ", got " \
		  << VAR << " at File " << __FILE__ << ":" << __LINE__ << std::endl; \
	d_bad_state = true; \
	rtt_dsxx::toss_cookies_ptr( #VAR "==" #VAL, __FILE__, __LINE__ );	\
    } 


//! A safe pointer-like class for scalars
/*! 

  Except for strange situations, you should use the \c DBC_Ptr macro in
  DBC_Ptr.hh instead of using this class explicitly.

  void pointers are not reference counted!

  \sa DBC_Ptr.hh
 */
template<class T>
class Safe_Ptr 
{
  public:
    //! Default constructor.
    /*! The pointer is initialized to \c NULL. */
    Safe_Ptr() 
	: d_p(0)
	, d_r(0)
	, d_bad_state(false) 
    { 
    }

    // Explicit constructor for type T *.
    inline explicit Safe_Ptr(T *p_in);

    // Copy constructor for Safe_Ptr<T>.
    inline Safe_Ptr(const Safe_Ptr<T> &sp_in);

    //! Destructor, memory is released when count goes to zero.
    ~Safe_Ptr() { if(!d_bad_state) decrement_rc(); }

    // Assignment operator for type T *.
    inline Safe_Ptr<T>& operator=(T *p_in);

    // Assignment operator for type Safe_Ptr<T>.
    inline Safe_Ptr<T>& operator=(const Safe_Ptr<T> sp_in);

    // Explicit constructor for type X *.
    template<class X>
    inline explicit Safe_Ptr(X *px_in);

    // Copy constructor for Safe_Ptr<X>.
    template<class X>
    inline Safe_Ptr(const Safe_Ptr<X> &spx_in);

    // Assignment operator for type X *.
    template<class X>
    inline Safe_Ptr<T>& operator=(X *px_in);

    // Assignment operator for type Safe_Ptr<X>.
    template<class X>
    inline Safe_Ptr<T>& operator=(const Safe_Ptr<X> spx_in);

    //! Set the pointer to \c NULL
    void release_data() { decrement_rc(); d_r = 0; d_p = 0; }

    //! Call \c delete on the pointer
    void delete_data();

    //! Access operator.
    T* operator->() const 
    { DBC_Check(d_p); DBC_Check(!d_r->deleted); return d_p; }

    //! Dereference operator.
    T& operator*() const 
    { DBC_Check(d_p); DBC_Check(!d_r->deleted); return *d_p; }

    size_t ref_count() const { return (d_r?d_r->refs:0); }

    //! Boolean conversion operator.
    operator bool() const 
    { return d_p != 0; }

    //! Operator not.
    bool operator!() const 
    { return d_p == 0; }

    //! Equality operator for T*.
    bool operator==(const T *p_in) const 
    { return d_p == p_in; }

    //! Inequality operator for T*.
    bool operator!=(const T *p_in) const 
    { return d_p != p_in; }

    //! Equality operator for Safe_Ptr<T>.
    bool operator==(const Safe_Ptr<T> &sp_in) const 
    { return d_p == sp_in.d_p; }

    //! Inequality operator for Safe_Ptr<T>.
    bool operator!=(const Safe_Ptr<T> &sp_in) const 
    { return d_p != sp_in.d_p; }


  private: 
    // >>> DATA

    //! Pointer to reference counter.
    Safe_Ptr_Ref* d_r;

    //! Pointer to the data
    T* d_p;

    //! This is true IFF we have thrown an exception.  
    /*! This is used to "turn off" the destructor */
    mutable bool d_bad_state;

  private:
    // >>> IMPLEMENTATION

    // Decrement_Rc the pointer.
    inline void decrement_rc();

    //! All derivatives of SP are friends. 
    template<class X> friend class Safe_Ptr;


};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Do equality check with a decrement_rc pointer.
 */
template<class T> inline bool
operator==(const T *pt, const Safe_Ptr<T> &sp)
{
    return sp == pt;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do inequality check with a decrement_rc pointer.
 */
template<class T> inline bool
operator!=(const T *pt, const Safe_Ptr<T> &sp)
{
    return sp != pt;
}


//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Explicit constructor for type T *.
 *
 * This constructor is used to initialize a Safe_Ptr with a pointer, ie.
 * \code
 *     Foo *f = new Foo;
 *     Safe_Ptr<Foo> spf(f);   // f now owned by spf
 *     Safe_Ptr<Foo> spf2 = f; // error! does not do implicit conversion
 * \endcode
 *
 * \param p_in pointer to type T
 */
template<class T>
Safe_Ptr<T>::Safe_Ptr(T *p_in)
    : d_p(p_in)
    , d_r(new Safe_Ptr_Ref)
    , d_bad_state(false)
    
{
    if(d_p) 
    {
	d_r = new Safe_Ptr_Ref;
	DBC_Check(d_r);
	DBC_Check_Eq(d_r->refs, 1);
    }
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
 *     Safe_Ptr<Base> spb(new Derived); // spb base class Safe_Ptr to Derived type
 *     Derived *d = new Derived;
 *     Safe_Ptr<Base> spb2(d);          // different syntax
 *     Safe_Ptr<Base> spb3 = d;         // error! no implicit conversion
 *     
 *     Derived *d2;
 *     Safe_Ptr<Base> spb4(d2);         // error! cannot initialize with NULL
 *                                // pointer of different type than T
 * \endcode
 * The pointer to X must not be equal to NULL.  The Safe_Ptr owns the pointer when
 * it is constructed.
 *
 * \param px_in pointer to type X that is convertible to T *
 */
template<class T>
template<class X>
Safe_Ptr<T>::Safe_Ptr(X *px_in)
    : d_p(0)
    , d_r(0)
    , d_bad_state(false)
{
    DBC_Check(px_in);

    // make a dynamic cast to check that we can cast between X and T
    T *np = dynamic_cast<T *>(px_in);

    // check that we have made a successfull cast if px exists
    DBC_Check(np);

    // assign the pointer and reference
    d_p = np;
    d_r = new Safe_Ptr_Ref;

    DBC_Check(d_p);
    DBC_Check(d_r);
    DBC_Check_Eq(d_r->refs, 1);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for Safe_Ptr<T>.
 *
 * \param sp_in smart pointer of type Safe_Ptr<T>
 */
template<class T>
Safe_Ptr<T>::Safe_Ptr(const Safe_Ptr<T> &sp_in)
    : d_p(sp_in.d_p)
    , d_r(sp_in.d_r)
    , d_bad_state(false)
{
    // advance the reference to T
    if(d_r)
	++d_r->refs;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor for Safe_Ptr<X>.
 *
 * This copy constructor requires that X * is convertible to T * through a
 * dynamic_cast.  The pointer in spx_in can point to NULL; however, it must
 * still be convertible to T * through a dynamic_cast.
 *
 * \param spx_in smart pointer of type Safe_Ptr<X>
 */
template<class T>
template<class X>
Safe_Ptr<T>::Safe_Ptr(const Safe_Ptr<X> &spx_in)
    : d_p(0)
    , d_r(0)
    , d_bad_state(false)
{
    DBC_Check(spx_in.d_r);

    // make a pointer to T *
    T *np = dynamic_cast<T *>(spx_in.p);
    DBC_Check(spx_in.d_p ? np != 0 : true);

    // assign the pointer and reference
    d_p = np;
    d_r = spx_in.d_r;
    
    // advance the reference to T
    if(d_r)
	++d_r->refs;
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
 *     Safe_Ptr<Foo> f;         // has reference to NULL pointer
 *     f      = new Foo;  // now has 1 count of Foo
 *     Foo *g = new Foo; 
 *     f      = g;        // f's original pointer to Foo is deleted
 *                        // because count goes to zero; f now has
 *                        // 1 reference to g
 * \endcode
 * 
 * \param p_in pointer to T
 */
template<class T>
Safe_Ptr<T>& Safe_Ptr<T>::operator=(T *p_in)
{
    // check if we already own this pointer
    if (d_p == p_in)
	return *this;

    // next decrement_rc the existing pointer
    decrement_rc();
    
    // now make add p_in to this pointer and make a new reference to it
    d_p = p_in;

    if(d_p)
    {
	d_r = new Safe_Ptr_Ref;
	DBC_Check_Eq(d_r->refs,1);
    }
    else
    {
	d_r = 0;
    }
    return *this;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type X *.
 *
 * This assignment requires that X * is convertible to T * through a
 * dynamic_cast.  It follows the same principle as Safe_Ptr(X*); however,
 * this is assignment:
 * \code
 *     Safe_Ptr<Base> b;
 *     b = new Derived;
 * \endcode
 * The pointer to X must not be equal to NULL.
 *
 * \param px_in pointer to type X * that is convertible to type T * through a
 * dynamic cast
 */
template<class T>
template<class X>
Safe_Ptr<T>& Safe_Ptr<T>::operator=(X *px_in)
{
    DBC_Check(px_in);

    // do a dynamic cast to ensure convertiblility between T* and X*
    T *np = dynamic_cast<T *>(px_in);
    DBC_Check(np);

    // now assign this to np (using previously defined assignment operator)
    *this = np;
    
    DBC_Check_Eq(d_r->refs,1);
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type Safe_Ptr<T>.
 *
 * \param sp_in smart pointer of type Safe_Ptr<T>
 */
template<class T>
Safe_Ptr<T>& Safe_Ptr<T>::operator=(const Safe_Ptr<T> sp_in)
{
    DBC_Check(sp_in.d_r);

    // see if they are equal
    if (this == &sp_in || d_p == sp_in.d_p)
    {
	DBC_Check(d_r == sp_in.d_r);
	return *this;
    }

    // decrement_rc the existing pointer
    decrement_rc();

    // assign p and r to sp_in
    d_p = sp_in.d_p;
    d_r = sp_in.d_r;

    // add the reference count and return
    if(d_r)
	++d_r->refs;
    return *this;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type Safe_Ptr<X>.
 *
 * This assignment requires that X * is convertible to T * through a
 * dynamic_cast.  The pointer in spx_in can point to NULL; however, it must
 * still be convertible to T * through a dynamic_cast.
 *
 * \param spx_in smart pointer of type Safe_Ptr<X>
 */
template<class T>
template<class X>
Safe_Ptr<T>& Safe_Ptr<T>::operator=(const Safe_Ptr<X> spx_in)
{
    DBC_Check(spx_in.d_r);

    // make a pointer to T *
    T *np = dynamic_cast<T *>(spx_in.d_p);
    DBC_Check(spx_in.d_p ? np != 0 : true);

    // check to see if we are holding the same pointer (and np is not NULL);
    // to NULL pointers to the same type are defined to be equal by the
    // standard 
    if (d_p == np && d_p)
    {
	// if the pointers are the same the reference count better be the
	// same 
	DBC_Check(d_r == spx_in.d_r);
	return *this;
    }

    // we don't need to worry about the case where p == np and np == NULL
    // because this combination is impossible; if np is NULL then it belongs
    // to a different smart pointer; in other words, if p == np and np ==
    // NULL then r != spx_in.r

    // free the existing pointer and reference
    decrement_rc();

    // assign new values
    d_p = np;
    d_r = spx_in.d_r;

    // advance the counter and return
    if(d_r)
	++d_r->refs;
    return *this;
}


template<class T> void
Safe_Ptr<T>::delete_data()
{
    if(d_r)
    {
	// Check for dangling pointers.  If delete_data() gets called during
	// the unwind from a memory leak exception thrown in decrement_rc,
	// d_r->refs could be zero.  d_bad_state should be set if this is
	// true, so just skip the delete.

	if(d_bad_state) return; 

	if(d_r->refs > 1)
	{
	    d_bad_state = true;
	    std::cout << "***DBC_Ptr error: dangling pointer\n\tSafe_Ptr<"
		      << typeid(T).name() << "> at addr " << this 
		      << "\n\tdeleted data at addr "
		      << d_p << " still referenced by " <<  d_r->refs - 1 
		      << " other Safe_Ptrs" << std::endl;
	    rtt_dsxx::toss_cookies_ptr("dangling pointer", __FILE__, __LINE__);
	}

	DBC_Check_Eq(d_r->refs,1);


	// Check for double deletes. 
	if(d_r->deleted)
	{
	    d_bad_state = true;
	    std::cout << "***DBC_Ptr error: double delete\n\tSafe_Ptr<"
		      << typeid(T).name() << "> at addr " << this 
		      << "\n\tcalled delete_data on previous deleted addr "
		      << d_p << std::endl;
	    rtt_dsxx::toss_cookies_ptr("double delete", __FILE__, __LINE__);
	}
	d_r->deleted = true;
	delete d_p;
    } else {
	// If there isn't a reference counter, there shouldn't be data.
	DBC_Check(!d_p);       
    }

}




//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Decrement the count and decrement_rc the pointer if count is zero.
 *
 * Note that it is perfectly acceptable to call delete on a NULL pointer.
 */
template<class T>
void Safe_Ptr<T>::decrement_rc()
{
    if(d_r)
    {
	DBC_Check(d_r->refs > 0);

	// if the count goes to zero then we free the reference.
	if (--d_r->refs == 0)
	{
	    if(!d_r->deleted)
	    {
		d_bad_state = true;
		std::cout << "***DBC_Ptr error: memory leak\n\tSafe_Ptr<"
			  << typeid(T).name() << "> at addr " << this 
			  << "\n\treleased last handle "
		    "to undeleted data at addr "
			  << d_p << std::endl;
		rtt_dsxx::toss_cookies_ptr("memory leak", __FILE__, __LINE__);
	    }
	    delete d_r; d_r = 0;
	}
    }
}

} // end namespace rtt_dsxx

#endif
