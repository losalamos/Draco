//----------------------------------*-C++-*----------------------------------//
// TstPCG_MatVec.cc
// Dave Nystrom
// Fri May  9 13:39:24 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "linalg/test/TstPCG_MatVec.hh"

//---------------------------------------------------------------------------//
// Destructor.
//---------------------------------------------------------------------------//

template<class T>
TstPCG_MatVec<T>::~TstPCG_MatVec() {}

//---------------------------------------------------------------------------//
// Constructor.
//---------------------------------------------------------------------------//

template<class T>
TstPCG_MatVec<T>::TstPCG_MatVec( int _nxs, int _nys )
    : PCG_MatVec<T>(), nxs(_nxs), nys(_nys) {}

//---------------------------------------------------------------------------//
// Evaluate matrix-vector product.
//---------------------------------------------------------------------------//

template<class T>
void TstPCG_MatVec<T>::MatVec( Mat1<T>& b, Mat1<T>&x )
{
    for( int ix = 0; ix < nxs; ix++ ) {
	for( int iy = 0; iy < nys; iy++ ) {
	    int indva   = ix     + nxs*(iy);
	    int indvqr0 = ix     + nxs*(iy);
	    int indvqr1 = ix - 1 + nxs*(iy);
	    int indvqr2 = ix + 1 + nxs*(iy);
	    int indvqr3 = ix     + nxs*(iy-1);
	    int indvqr4 = ix     + nxs*(iy+1);

	    b(indva) = 4.0*x(indvqr0);

	    if( ix !=     0 ) b(indva) -= x(indvqr1);
	    if( ix != nxs-1 ) b(indva) -= x(indvqr2);
	    if( iy !=     0 ) b(indva) -= x(indvqr3);
	    if( iy != nys-1 ) b(indva) -= x(indvqr4);
	}
    }
}

//---------------------------------------------------------------------------//
//                              end of TstPCG_MatVec.cc
//---------------------------------------------------------------------------//

