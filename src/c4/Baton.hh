//----------------------------------*-C++-*----------------------------------//
// Baton.hh
// Geoffrey Furnish
// Thu Apr 13 15:54:13 1995
//---------------------------------------------------------------------------//
// @> A spinlock which transmits a value between nodes.
//---------------------------------------------------------------------------//

#ifndef __c4_Baton_hh__
#define __c4_Baton_hh__

#include "c4/NodeInfo.hh"

C4_NAMESPACE_BEG

//===========================================================================//
// class Baton<T> - Spinlock which passes a value between nodes.

// The Baton<T> class is used to handle a case where a node needs to do work
// which depends on an answer from the preceeding node.  <More blurb...>
//
// IMPLEMENTATION NOTE:
// We really ought to be able to derive from the template parameter.  That
// would best represent the incremental addition of the baton "pass along"
// capability, while leaving all other aspects of the class available in the
// usual way.  Unfortunately, since pathetic Cfront doesn't support
// derivation from a template parameter, we use a sorry substitue instead.
//===========================================================================//

template<class T>
class Baton : public NodeInfo {

    T v;

// Disable object copying, since we don't want the sychronization point to
// float around in the call tree.

    Baton( const Baton<T>& ) {}
    Baton<T>& operator=( const Baton<T>& b ) { return *this; }

  public:
    Baton();
    Baton( T _v );
    ~Baton();

    T& operator=( const T& u )  { v = u;  return v; }
    T& operator+=( const T& u ) { v += u; return v; }
    T& operator-=( const T& u ) { v -= u; return v; }
    T& operator*=( const T& u ) { v *= u; return v; }
    T& operator/=( const T& u ) { v /= u; return v; }

    operator T () const { return v; }
};

C4_NAMESPACE_END

#endif                          // __c4_Baton_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Baton.hh
//---------------------------------------------------------------------------//
