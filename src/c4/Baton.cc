//----------------------------------*-C++-*----------------------------------//
// Baton.cc
// Geoffrey Furnish
// Thu Apr 13 15:54:14 1995
//---------------------------------------------------------------------------//
// @> A spinlock which transmitts a value between nodes.
//---------------------------------------------------------------------------//

#include "c4/Baton.hh"
#include "c4/BSwap.hh"

//---------------------------------------------------------------------------//
// Default ctor.  The type T had better have a default ctor too, or you need
// to avoid using this method.  Anyway, this is for defered assignment.
//---------------------------------------------------------------------------//

template<class T>
Baton<T>::Baton()
{
    if (node > 0) {
	BSwap<T> bs(node-1);

	bs.recv(v);
    }
}

//---------------------------------------------------------------------------//
// Constructor.  If not node zero, wait for the previous node to tell us his
// result. 
//---------------------------------------------------------------------------//

template<class T>
Baton<T>::Baton( T _v )
    : v(_v)
{
    if (node > 0) {
	BSwap<T> bs(node-1);

	bs.recv(v);
    }
}

//---------------------------------------------------------------------------//
// Destructor.  If not lastnode, send our result to the next node, enabling
// him to procede.
//---------------------------------------------------------------------------//

template<class T>
Baton<T>::~Baton()
{
    if (node < lastnode) {
	BSwap<T> bs(node+1);

	bs.send(v);
    }
}

//---------------------------------------------------------------------------//
//                              end of Baton.cc
//---------------------------------------------------------------------------//
