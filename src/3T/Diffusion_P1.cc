//----------------------------------*-C++-*----------------------------------//
// Diffusion_P1.cc
// Geoffrey M. Furnish
// Thu May 28 13:16:55 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/Diffusion_P1.hh"

template<class MT>
Diffusion_P1<MT>::Diffusion_P1( const Diffusion_DB& diffdb,
                                const SP<MT>& spm_, const pcg_DB& pcg_db )
    : MT::Coord_Mapper( spm_->get_Mesh_DB() ),
      Diffusion_DB( diffdb ),
      spm(spm_),
      A( static_cast<MT::Coord_Mapper>(*this), spm->get_diag_offsets() ),
      pcg_ctrl( pcg_db, ncp ),

//       vc( spm->get_vc() ),

      dx( ncx ), dy( ncy ), dz( ncz )
{
// Fetch the face locations from the mesh.
    const Mat1<double>& xf = spm->get_xf();
    const Mat1<double>& yf = spm->get_yf();
    const Mat1<double>& zf = spm->get_zf();

// Calculate the cell widths from the locations of the faces.

    for( int i=0; i < ncx+1; i++ )
        dx(i) = xf(i+1) - xf(i);

    for( int j=0; j < ncy+1; j++ )
        dy(j) = yf(j+1) - yf(j);

    for( int k=0; k < ncz+1; k++ )
        dz(k) = zf(k+1) - zf(k);
}

template<class MT>
void Diffusion_P1<MT>::solve( const fcdsf& D,
                              const ccsf& sigmaabar,
                              const ccsf& Qbar_r,
                              const fcdsf& Fprime,
                              const bssf& f_b,
                              ccsf& phi,
                              fcdsf& F )
{
}

template<class MT>
void Diffusion_P1<MT>::solve( const fcdsf& D,
                              const ccsf& sigmaabar,
                              const ccsf& Qbar_r,
                              ccsf& phi )
{
}

//---------------------------------------------------------------------------//
//                              end of Diffusion_P1.cc
//---------------------------------------------------------------------------//
