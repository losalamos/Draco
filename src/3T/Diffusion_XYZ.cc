//----------------------------------*-C++-*----------------------------------//
// Diffusion_XYZ.cc
// Geoffrey Furnish
// Thu Dec  4 09:30:51 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/Diffusion_XYZ.hh"

//---------------------------------------------------------------------------//
// Initialize the diffusion solver.
//---------------------------------------------------------------------------//

template<class MT>
Diffusion_XYZ<MT>::Diffusion_XYZ( const SP<MT>& spm_ )
    : MT::Coord_Mapper( spm_->get_Mesh_DB() ),
      spm(spm_),
      A( ncp, nct )
{}

//---------------------------------------------------------------------------//
// Solve the diffusion equation.  
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_XYZ<MT>::solve( const typename MT::fcdsf& D,
			       const typename  MT::cell_array& rhs,
			       double dt, typename MT::cell_array& x )
{
// Clear the matrix.

    A = 0.;

// Now initialize the matrix.

    for( int i=0; i < ncp; i++ )
	A(i,goff+i) = 1.;

// Now solve the matrix equation A.x = rhs.

    x = rhs;
}

//---------------------------------------------------------------------------//
//                              end of Diffusion_XYZ.cc
//---------------------------------------------------------------------------//
