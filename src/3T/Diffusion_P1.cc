//----------------------------------*-C++-*----------------------------------//
// Diffusion_P1.cc
// Geoffrey M. Furnish
// Thu May 28 13:16:55 1998
//---------------------------------------------------------------------------//
// @> Class for solving P1 Diffusion equations to support P1-3T package.
//---------------------------------------------------------------------------//

#include "3T/Diffusion_P1.hh"

#include "c4/C4_Req.hh"

template<class MT>
Diffusion_P1<MT>::Diffusion_P1( const Diffusion_DB& diffdb,
                                const SP<MT>& spm_, const pcg_DB& pcg_db )
    : MT::Coord_Mapper( spm_->get_Mesh_DB() ),
      Diffusion_DB( diffdb ),
      spm(spm_),
      A( static_cast<MT::Coord_Mapper>(*this), spm->get_diag_offsets() ),
      pcg_ctrl( pcg_db, ncp ),

      dx( ncx ), dy( ncy ), dz( ncz ),

      Dprime( spm ), Dtwidle( spm ), Dhat( spm )
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
// Calculate Dprime.  This is a fcdsf which represents the D for the same
// face, as seen by the cell on the other side of the face.  For cells
// interior to a processor subdomain, this is easy.  For cells on the
// boundaries of processor subdomains, message communication is required.

    calculate_Dprime( D );

// Calculate Dtwidle.  This is a fcdsf computed from D and Dprime.  Note that 
// Dtwidle can only be calculated on the interior faces, since it involves
// geometric data  which is unavailable for the exterior faces.

    calculate_Dtwidle( D, Dprime );

// Calculate Dhat.  This is an fcdsf.  For interior faces, it is Dtwidle.
// For faces on the exterior boundary of the system, it is a function of
// Dtwidle and the boundary conditions.

    Dhat = Dtwidle;
    calculate_Dhat_on_boundaries();

// Calculate A from sigmaabar, Dhat, and geometric data.

// Calculate Ftwidle.  This is an fcdsf computed from D, Dprime, and Fprime.
// Like Dtwidle, it can only be calculated on interior faces because it
// involves gemetric data which does not exist for exterior faces.

// Calculate Fhat.  This is an fcdsf.  Like Dhat, Fhat is Ftwidle on interior 
// faces, and is a function of D, Fprime, f_b, and various geometric data on
// exterior faces.

// Calculate b from Qbar_r, Fhat, and various geometric data.

// Assert that A is symmetric.

// Solve A.phi = b using PCG.

// Calculate new values of F from phi, Dtwidle, Ftwidle, D and boundary
// data. 

}

template<class MT>
void Diffusion_P1<MT>::solve( const fcdsf& D,
                              const ccsf& sigmaabar,
                              const ccsf& Qbar_r,
                              ccsf& phi )
{
// Calculate Dprime.  This is a fcdsf which represents the D for the same
// face, as seen by the cell on the other side of the face.  For cells
// interior to a processor subdomain, this is easy.  For cells on the
// boundaries of processor subdomains, message communication is required.

    calculate_Dprime( D );

// Calculate Dtwidle.  This is a fcdsf computed from D and Dprime.  Note that 
// Dtwidle can only be calculated on the interior faces, since it involves
// geometric data  which is unavailable for the exterior faces.

    calculate_Dtwidle( D, Dprime );

// Calculate Dhat.  This is an fcdsf.  For interior faces, it is Dtwidle.
// For faces on the exterior boundary of the system, it is 0, since this
// solve method is used when solving the conduction equations, which always
// user reflective boundary conditions, for which alpha=0, beta=1.

// Calculate A from sigmaabar, Dhat, and geometric data.

// Calculate b from Qbar_r.  Do not include Fhat terms in this case.

// Solve A.phi = b using PCG.
}

//---------------------------------------------------------------------------//
// This method calculcates Dprime, which in Randy's notation, is
// $D^{c'}_{f'}$, the value of D on a face as seen from the cell on the other 
// side of that face.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_Dprime( const fcdsf& D )
{
    using namespace C4;

// These are the receive buffers.
    Mat2<double> Dbot( ncx, ncy ), Dtop( ncx, ncy );

// These are the receive request handles.
    C4_Req brreq, trreq;

// Post receives for inbound data.
    if (node > 0)
        RecvAsync( brreq, &Dbot(0,0), ncx*ncy, node-1 );
    if (node < lastnode)
        RecvAsync( trreq, &Dtop(0,0), ncx*ncy, node+1 );

// Allocate send buffers.
    Mat2<double> sbot( ncx, ncy ), stop( ncx, ncy );

// Fill up the send buffers
    for( int j=0; j < ncy; j++ )
        for( int i=0; i < ncx; i++ ) {
            sbot( i, j ) = D( local_cell_index(i,j,zoff), 4 );
            stop( i, j ) = D( local_cell_index(i,j,zoff+nczp-1), 5 );
        }

// Send the send buffers.
    if (node > 0)
        Send( sbot.begin(), sbot.size(), node-1 );
    if (node < lastnode)
        Send( stop.begin(), stop.size(), node+1 );

// Now wait on the inbound data, because we need it to continue.
    brreq.wait();
    trreq.wait();

// Clear out Dprime (should be unnecessary...).
    Dprime = 0.;

// Loop through cells and set value of each face.
    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index(i,j,k);

            // Our left face = right face of our neighbor on left.
                if (i == 0)
                    Dprime( c, 0 ) = D( c, 0 );
                else
                    Dprime( c, 0 ) = D( local_cell_index( i-1, j, k ), 1 );

            // Our rightt face = left face of our neighbor on right.
                if (i == ncx-1)
                    Dprime( c, 1 ) = D( c, 1 );
                else
                    Dprime( c, 1 ) = D( local_cell_index( i+1, j, k ), 0 );

            // Our front face = back face of our neighbor to the front.
                if (j == 0)
                    Dprime( c, 2 ) = D( c, 2 );
                else
                    Dprime( c, 2 ) = D( local_cell_index( i, j-1, k ), 3 );

            // Our back face = front face of our neighbor to the rear.
                if (j == ncy-1)
                    Dprime( c, 3 ) = D( c, 3 );
                else
                    Dprime( c, 3 ) = D( local_cell_index( i, j+1, k ), 2 );

            // Our bottom face = top face of our neighbor below.
                if (k == 0)
                    Dprime( c, 4 ) = D( c, 4 );
                else if (k == zoff)
                    Dprime( c, 4 ) = Dbot( i, j );
                else
                    Dprime( c, 4 ) = D( local_cell_index( i, j, k-1 ), 5 );

            // Our top face = bottom face of neighbor above.
                if (k == ncz-1)
                    Dprime( c, 5 ) = D( c, 5 );
                else if (k == zoff+nczp-1)
                    Dprime( c, 5 ) = Dtop( i, j );
                else
                    Dprime( c, 5 ) = D( local_cell_index( i, j, k+1 ), 4 );
            }
}

//---------------------------------------------------------------------------//
// Calculates the value of Dtwidle on interior faces, sets the others to 0.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_Dtwidle( const fcdsf& D, const fcdsf& Dp )
{
    Dtwidle = 0.;

// Loop through cells and set value of each face.
    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index(i,j,k);

            // Handle x faces.
                if (i > 0)
                    Dtwidle( c, 0 ) = 2.*D(c,0)*Dp(c,0)*dx(i) /
                        ( D(c,0)*dx(i-1) + Dp(c,0)*dx(i) );
                if (i < ncx-1)
                    Dtwidle( c, 1 ) = 2.*D(c,1)*Dp(c,1)*dx(i) /
                        ( D(c,1)*dx(i+1) + Dp(c,1)*dx(i) );

            // Handle y faces.
                if (j > 0)
                    Dtwidle( c, 2 ) = 2.*D(c,2)*Dp(c,2)*dy(j) /
                        ( D(c,2)*dy(j-1) + Dp(c,0)*dy(j) );
                if (j < ncy-1)
                    Dtwidle( c, 3 ) = 2.*D(c,3)*Dp(c,3)*dy(j) /
                        ( D(c,3)*dy(j+1) + Dp(c,3)*dy(j) );

            // Handle z faces.
                if (k > 0)
                    Dtwidle( c, 4 ) = 2.*D(c,4)*Dp(c,4)*dz(k) /
                        ( D(c,4)*dz(k-1) + Dp(c,4)*dz(k) );
                if (k < ncz-1)
                    Dtwidle( c, 5 ) = 2.*D(c,5)*Dp(c,4)*dz(k) /
                        ( D(c,5)*dz(k+1) + Dp(c,5)*dz(k) );
            }
}

//---------------------------------------------------------------------------//
// Loop over boundary faces and apply boundary conditions to Dhat.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_Dhat_on_boundaries()
{
}

//---------------------------------------------------------------------------//
//                              end of Diffusion_P1.cc
//---------------------------------------------------------------------------//
