//----------------------------------*-C++-*----------------------------------//
// SP.hh
// Geoffrey Furnish
// 2 December 1994
//---------------------------------------------------------------------------//
// @> A "smart pointer" facility.
//---------------------------------------------------------------------------//

#ifndef __ds_SP_hh__
#define __ds_SP_hh__

#include "Assert.hh"

NAMESPACE_DS_BEG

// autodoc: noprint SPref

//===========================================================================//
// class SP<T> - A templated "smart pointer" class.

// This class provides an encapsulation of a pointer of a specified type,
// enabling an accurate determination of when the pointer is in use and when
// it can be freed.  Consider: A function new's an object and return the
// pointer to that object as its return value.  Now it is the caller's
// responsibility to free the object.  What if the caller passes the pointer
// to other objects or functions?  What if it is not known which will be
// deleted first or last?
//
// Instead the function can return a "smart pointer".  This SP class uses
// reference counting to determine the number of current users of a pointer.
// Each time an SP goes out of scope, the reference count is decremented.
// When the last user of a pointer is done, the pointer is freed.
//
// Note:  I am calling this a "smart pointer", not a "safe pointer".  There
// are clearly ways you can hose this.  In particular, when you bind an SP<T>
// to a T*, you yield all rights to the T*.  You'd better not squirrel the
// bare pointer away somewhere and expect to clandestinely use it in other
// ways or places--death will be sure to follow.  Consequently then, the
// safest way to use this smart pointer, is to bind it to the contained
// pointer and then always use the smart pointer.  Immediately returning the
// smart pointer as a return value, allowing the orriginal bare pointer to go
// out of scope never to be seen again, is one good example of how to use
// this. 
//
// In other words, as with other parts of DS++, the SP<T> class is intended to
// facilitate legitimate use.  I am not bending over backwards to make certain
// it is inviolable.  The user is supposed to abstain from heinous misuse.
//===========================================================================//

class SPref {
    int refs;

    SPref( int _refs =1 ) : refs(_refs) {}

    template<class T> friend class SP;
};

template<class T>
class SP {
    T *p;
    SPref *r;

    void validate() const  { Insist( p, "No dumb pointer bound "
				     "to this smart pointer." ); }

    void detach()
    {
	if ( --r->refs == 0) {
	    delete p;
	    delete r;
	}
    }

  public:
    SP() : p(0) { r = new SPref; }

    explicit SP( T *_p ) : p(_p) { r = new SPref; }

    template<class X>
    SP( X *px )
    {
	T *np = dynamic_cast<T *>( px );
	if (!np) throw "Incompatible dumb pointer type.";
	p = np;
	r = new SPref;
    }

    SP( const SP<T>& sp ) : p(sp.p), r(sp.r) { r->refs++; }

    template<class X>
    SP( const SP<X>& spx )
    {
	X *px = spx.p;
	T *np = dynamic_cast<T *>( px );
	if (!np) throw "Incompatible smart pointer types.";
	p = np;
	r = spx.r;
	r->refs++;
    }

    ~SP() { detach(); }

    SP& operator=( T *np )
    {
	if (p == np) return *this; // It could happen.
	detach();
	p = np;
	r = new SPref(1);
	return *this;
    }

    template<class X>
    SP& operator=( X *px )
    {
	T *np = dynamic_cast<T *>( px );
	if (!np) throw "Incompatible smart pointer types.";
	return *this = np;
    }

    SP& operator=( const SP<T>& sp )
    {
	if ( &sp == this ) return *this;
	if ( p == sp.p ) return *this;
	detach();
	p = sp.p;
	r = sp.r;
	r->refs++;
	return *this;
    }

    template<class X>
    SP& operator=( const SP<X>& spx )
    {
	X *px = spx.p;
	T *np = dynamic_cast<T *>( px );
	if (!np) throw "Incompatible smart pointer types.";
	if (p == np) {
	// If we have the same pointer, we'd better be sharing the reference
	// holder too!
	    Assert( r == spx.r );
	    return *this;	// It could happen.
	}
	detach();
	p = np;
	r = spx.r;
	r->refs++;
	return *this;
    }
    
    T *operator->() const { validate(); return p; }
    T& operator*() const  { validate(); return *p; }

    T *bp () const { return p; } // Better know what you're doing!

    operator bool() const { return p != 0; }
    bool operator!() const { return p == 0; }

    bool operator==( const T *pt ) const { return p == pt; }
    bool operator!=( const T *pt ) const { return p != pt; }

    bool operator==( const SP<T>& sp ) const { return p == sp.p; }
    bool operator!=( const SP<T>& sp ) const { return p != sp.p; }

    template<class X> friend class SP;
};

template<class T>
bool operator==( const T *pt, const SP<T>& sp )
{
    return pt == sp.bp();
}

template<class T>
bool operator!=( const T *pt, const SP<T>& sp )
{
    return pt != sp.bp();
}

NAMESPACE_DS_END

#endif                          // __ds_SP_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/SP.hh
//---------------------------------------------------------------------------//
