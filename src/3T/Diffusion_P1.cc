//----------------------------------*-C++-*----------------------------------//
// Diffusion_P1.cc
// Geoffrey M. Furnish
// Thu May 28 13:16:55 1998
//---------------------------------------------------------------------------//
// @> Class for solving P1 Diffusion equations to support P1-3T package.
//---------------------------------------------------------------------------//

#include "3T/Diffusion_P1.hh"

#include "c4/C4_Req.hh"

using dsxx::Mat1;
using dsxx::Mat2;

template<class MT>
Diffusion_P1<MT>::Diffusion_P1( const Diffusion_DB& diffdb,
                                const dsxx::SP<MT>& spm_, const pcg_DB& pcg_db )
    : MT::Coord_Mapper( spm_->get_Mesh_DB() ),
      Diffusion_DB( diffdb ),
      spm(spm_),
      A( static_cast<MT::Coord_Mapper>(*this), spm->get_diag_offsets() ),
      b( spm ),
      
      pcg_ctrl( pcg_db, ncp ),

      dx( ncx ), dy( ncy ), dz( ncz ),
      fdeltal( spm ),

      Dprime( spm ),      Dtwidle( spm ), Dhat( spm ),
      Fprimeprime( spm ), Ftwidle( spm ), Fhat( spm )
{
// Fetch the face locations from the mesh.
    const Mat1<double>& xf = spm->get_xf();
    const Mat1<double>& yf = spm->get_yf();
    const Mat1<double>& zf = spm->get_zf();

// Calculate the cell widths from the locations of the faces.

    for( int i=0; i < ncx; i++ )
        dx(i) = xf(i+1) - xf(i);

    for( int j=0; j < ncy; j++ )
        dy(j) = yf(j+1) - yf(j);

    for( int k=0; k < ncz; k++ )
        dz(k) = zf(k+1) - zf(k);

// Calculate the delta-l that goes with each face.
    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index(i,j,k);
                fdeltal(c,0) = fdeltal(c,1) = dx(i);
                fdeltal(c,2) = fdeltal(c,3) = dy(j);
                fdeltal(c,4) = fdeltal(c,5) = dz(k);
            }

// Construct the objects we'll be neededing later to do the pcg solve.

    spmv = new MatVec_3T< Diffusion_P1<MT> >( this );
    precond = new PreCond< Diffusion_P1<MT> >( ncp, pc_meth, this );
}

template<class MT>
void Diffusion_P1<MT>::solve( ccsf& phi,
                              fcdsf& F,
			      const fcdsf& D,
                              const ccsf& sigmaabar,
                              const ccsf& Qbar_r,
                              const fcdsf& Fprime,
                              const bssf& f_b )
{
// Calculate Dprime.  This is a fcdsf which represents the D for the same
// face, as seen by the cell on the other side of the face.  For cells
// interior to a processor subdomain, this is easy.  For cells on the
// boundaries of processor subdomains, message communication is required.

    calculate_opposing_face( D, Dprime );

// Calculate Dtwidle.  This is a fcdsf computed from D and Dprime.  Note that 
// Dtwidle can only be calculated on the interior faces, since it involves
// geometric data  which is unavailable for the exterior faces.

    calculate_Dtwidle( D, Dprime );

// Calculate Dhat.  This is an fcdsf.  For interior faces, it is Dtwidle.
// For faces on the exterior boundary of the system, it is a function of
// Dtwidle and the boundary conditions.

    Dhat = Dtwidle;
    calculate_Dhat_on_boundaries( D );

// Calculate A from sigmaabar, Dhat, and geometric data.

    calculate_A( sigmaabar );

// Calculate Fprimeprime.  Like Dprime compared to D, Fprimeprime is the
// value of Fprime as seen by the cell on the other side of the face.

    calculate_opposing_face( Fprime, Fprimeprime );

// Calculate Ftwidle.  This is an fcdsf computed from D, Dprime, and Fprime.
// Like Dtwidle, it can only be calculated on interior faces because it
// involves gemetric data which does not exist for exterior faces.

    calculate_Ftwidle( D, Fprime );

// Calculate Fhat.  This is an fcdsf.  Like Dhat, Fhat is Ftwidle on interior 
// faces, and is a function of D, Fprime, f_b, and various geometric data on
// exterior faces.

    Fhat = Ftwidle;
    calculate_Fhat_on_boundaries( D, Fprime, f_b );

// Calculate b from Qbar_r, Fhat, and various geometric data.

    calculate_b( Qbar_r, Fhat );

// Solve A.phi = b using PCG.

    solve_A_phi_equals_b( phi );

// Calculate new values of F from phi, Dtwidle, Ftwidle, D and boundary
// data. 

    calculate_new_F( D, Fprime, f_b, phi, F );
}

template<class MT>
void Diffusion_P1<MT>::solve( ccsf& phi,
			      const fcdsf& D,
                              const ccsf& sigmaabar,
                              const ccsf& Qbar_r )
{
// Calculate Dprime.  This is a fcdsf which represents the D for the same
// face, as seen by the cell on the other side of the face.  For cells
// interior to a processor subdomain, this is easy.  For cells on the
// boundaries of processor subdomains, message communication is required.

    calculate_opposing_face( D, Dprime );

// Calculate Dtwidle.  This is a fcdsf computed from D and Dprime.  Note that 
// Dtwidle can only be calculated on the interior faces, since it involves
// geometric data  which is unavailable for the exterior faces.

    calculate_Dtwidle( D, Dprime );

// Calculate Dhat.  This is an fcdsf.  For interior faces, it is Dtwidle.
// For faces on the exterior boundary of the system, it is 0, since this
// solve method is used when solving the conduction equations, which always
// user reflective boundary conditions, for which alpha=0, beta=1.

    Dhat = Dtwidle;

// Calculate A from sigmaabar, Dhat, and geometric data.

    calculate_A( sigmaabar );

// Calculate b from Qbar_r.  Do not include Fhat terms in this case.

    calculate_b( Qbar_r );

// Solve A.phi = b using PCG.
}

//---------------------------------------------------------------------------//
// This method calculcates the value of a discontinuous face centered
// quantity on the opposite side of the face.  That is, the value of the face 
// as seen by the cell on the other side of the face.  For cells interior to
// a processor's subdomain, this info can be determined trivially.  However,
// for cells on the exterior boundary of a processor's subdomain, this data
// is only available on the other processor, so messaging has to be done to
// make it available over here.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_opposing_face( const fcdsf& X, fcdsf& Xprime )
{
    using namespace C4;

// These are the receive buffers.
    Mat2<double> Xbot( ncx, ncy ), Xtop( ncx, ncy );

// These are the receive request handles.
    C4_Req brreq, trreq;

// Post receives for inbound data.
    if (node > 0)
        RecvAsync( brreq, Xbot.begin(), Xbot.size(), node-1 );
    if (node < lastnode)
        RecvAsync( trreq, Xtop.begin(), Xtop.size(), node+1 );

// Allocate send buffers.
    Mat2<double> sbot( ncx, ncy ), stop( ncx, ncy );

// Fill up the send buffers
    for( int j=0; j < ncy; j++ )
        for( int i=0; i < ncx; i++ ) {
            sbot( i, j ) = X( local_cell_index(i,j,zoff), 4 );
            stop( i, j ) = X( local_cell_index(i,j,zoff+nczp-1), 5 );
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
    Xprime = 0.;

// Loop through cells and set value of each face.
    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index(i,j,k);

            // Our left face = right face of our neighbor on left.
                if (i == 0)
                    Xprime( c, 0 ) = X( c, 0 );
                else
                    Xprime( c, 0 ) = X( local_cell_index( i-1, j, k ), 1 );

            // Our rightt face = left face of our neighbor on right.
                if (i == ncx-1)
                    Xprime( c, 1 ) = X( c, 1 );
                else
                    Xprime( c, 1 ) = X( local_cell_index( i+1, j, k ), 0 );

            // Our front face = back face of our neighbor to the front.
                if (j == 0)
                    Xprime( c, 2 ) = X( c, 2 );
                else
                    Xprime( c, 2 ) = X( local_cell_index( i, j-1, k ), 3 );

            // Our back face = front face of our neighbor to the rear.
                if (j == ncy-1)
                    Xprime( c, 3 ) = X( c, 3 );
                else
                    Xprime( c, 3 ) = X( local_cell_index( i, j+1, k ), 2 );

            // Our bottom face = top face of our neighbor below.
                if (k == 0)
                    Xprime( c, 4 ) = X( c, 4 );
                else if (k == zoff)
                    Xprime( c, 4 ) = Xbot( i, j );
                else
                    Xprime( c, 4 ) = X( local_cell_index( i, j, k-1 ), 5 );

            // Our top face = bottom face of neighbor above.
                if (k == ncz-1)
                    Xprime( c, 5 ) = X( c, 5 );
                else if (k == zoff+nczp-1)
                    Xprime( c, 5 ) = Xtop( i, j );
                else
                    Xprime( c, 5 ) = X( local_cell_index( i, j, k+1 ), 4 );
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
void Diffusion_P1<MT>::calculate_Dhat_on_boundaries( const fcdsf& D )
{
// Calculate faces perpendicular to x (left and right).

    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
        {
            int cl = local_cell_index( 0, j, k );
            int cr = local_cell_index( ncx-1, j, k );

            Assert( Dhat(cl,0) == 0. );
            Dhat( cl, 0 ) = 2.*alpha_left*D(cl,0) /
                ( alpha_left - 2.*beta_left*D(cl,0)/dx(0) );

            Assert( Dhat(cr,1) == 0. );
            Dhat( cr, 1 ) = 2.*alpha_right*D(cr,1) /
                ( alpha_right - 2.*beta_right*D(cr,1)/dx(ncx-1) );
        }

// Calculate faces perpendicular to y (front and back).

    for( int k=zoff; k < zoff+nczp; k++ )
        for( int i=0; i < ncx; i++ )
        {
            int cf = local_cell_index( i, 0, k );
            int cb = local_cell_index( i, ncy-1, k );

            Assert( Dhat(cf,2) == 0. );
            Dhat( cf, 2 ) = 2.*alpha_front*D(cf,2) /
                ( alpha_front - 2.*beta_front*D(cf,2)/dy(0) );

            Assert( Dhat(cb,3) == 0. );
            Dhat( cb, 3 ) = 2.*alpha_back*D(cb,3) /
                ( alpha_back - 2.*beta_back*D(cb,3)/dy(ncy-1) );
        }

// Calculate faces perpendicular to z (bottom and top).

    if (zoff == 0)
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index( i, j, 0 );

                Assert( Dhat(c,4) == 0. );
                Dhat( c, 4 ) = 2.*alpha_bottom*D(c,4) /
                    ( alpha_bottom - 2.*beta_bottom*D(c,4)/dz(0) );
            }

    if (zoff+nczp == ncz)
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index( i, j, ncz-1 );

                Assert( Dhat(c,5) == 0. );
                Dhat( c, 5 ) = 2.*alpha_top*D(c,5) /
                    ( alpha_top - 2.*beta_top*D(c,5)/dz(ncz-1) );
            }
}

//---------------------------------------------------------------------------//
// Calculates the value of Ftwidle on interior faces, sets the others to 0.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_Ftwidle( const fcdsf& D, const fcdsf& Fprime )
{
    Ftwidle = 0.;               // should be unnecessary.

// Loop through cells and set value of each face.
    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index(i,j,k);

            // Handle x faces.
                if (i > 0)
                    Ftwidle( c, 0 ) =
                        ( Dprime(c,0)*dx(i)*Fprime(c,0) -
                          D(c,0)*dx(i-1)*Fprimeprime(c,0) ) /
                        ( D(c,0)*dx(i-1) + Dprime(c,0)*dx(i) );
                if (i < ncx-1)
                    Ftwidle( c, 1 ) =
                        ( Dprime(c,1)*dx(i)*Fprime(c,1) -
                          D(c,1)*dx(i+1)*Fprimeprime(c,1) ) /
                        ( D(c,1)*dx(i+1) + Dprime(c,1)*dx(i) );

            // Handle y faces.
                if (j > 0)
                    Ftwidle( c, 2 ) =
                        ( Dprime(c,2)*dy(j)*Fprime(c,2) -
                          D(c,2)*dy(j-1)*Fprimeprime(c,2) ) /
                        ( D(c,2)*dy(j-1) + Dprime(c,2)*dy(j) );
                if (j < ncy-1)
                    Ftwidle( c, 3 ) =
                        ( Dprime(c,3)*dy(j)*Fprime(c,3) -
                          D(c,3)*dy(j+1)*Fprimeprime(c,3) ) /
                        ( D(c,3)*dy(j+1) + Dprime(c,3)*dy(j) );

            // Handle z faces.
                if (k > 0)
                    Ftwidle( c, 4 ) =
                        ( Dprime(c,4)*dz(k)*Fprime(c,4) -
                          D(c,4)*dz(k-1)*Fprimeprime(c,4) ) /
                        ( D(c,4)*dz(k-1) + Dprime(c,4)*dz(k) );
                if (k < ncz-1)
                    Ftwidle( c, 5 ) =
                        ( Dprime(c,5)*dz(k)*Fprime(c,5) -
                          D(c,5)*dz(k+1)*Fprimeprime(c,5) ) /
                        ( D(c,5)*dz(k+1) + Dprime(c,5)*dz(k) );
            }
}

//---------------------------------------------------------------------------//
// Calculate exterior boundary face contributions to Fhat.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_Fhat_on_boundaries( const fcdsf& D,
                                                     const fcdsf& Fprime,
                                                     const bssf& f_b )
{
#ifdef RMR
    const Mat2<double>& f_b_left   = f_b.face(0);
    const Mat2<double>& f_b_right  = f_b.face(1);
    const Mat2<double>& f_b_front  = f_b.face(2);
    const Mat2<double>& f_b_back   = f_b.face(3);
    const Mat2<double>& f_b_bottom = f_b.face(4);
    const Mat2<double>& f_b_top    = f_b.face(5);

// Calculate faces perpendicular to x (left and right).

    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
        {
            int cl = local_cell_index( 0, j, k );
            int cr = local_cell_index( ncx-1, j, k );

            Assert( Fhat(cl,0) == 0. );
            Fhat( cl, 0 ) =
                ( alpha_left*Fprime(cl,0) -
                  2.*D(cl,0)*f_b_left(j,k)/dx(0) ) /
                ( alpha_left - 2.*beta_left*D(cl,0)/dx(0) );

            Assert( Fhat(cr,1) == 0. );
            Fhat( cr, 1 ) =
                ( alpha_right*Fprime(cr,1) -
                  2.*D(cr,1)*f_b_right(j,k)/dx(ncx-1) ) /
                ( alpha_right - 2.*beta_right*D(cr,1)/dx(ncx-1) );
        }

// Calculate faces perpendicular to y (front and back).

    for( int k=zoff; k < zoff+nczp; k++ )
        for( int i=0; i < ncx; i++ )
        {
            int cf = local_cell_index( i, 0, k );
            int cb = local_cell_index( i, ncy-1, k );

            Assert( Fhat(cf,2) == 0. );
            Fhat( cf, 2 ) = 
                ( alpha_front*Fprime(cf,2) -
                  2.*D(cf,2)*f_b_front(i,k)/dy(0) ) /
                ( alpha_front - 2.*beta_front*D(cf,2)/dy(0) );

            Assert( Fhat(cb,3) == 0. );
            Fhat( cb, 3 ) = 
                ( alpha_back*Fprime(cb,3) -
                  2.*D(cb,3)*f_b_back(i,k)/dy(ncy-1) ) /
                ( alpha_back - 2.*beta_back*D(cb,3)/dy(ncy-1) );
        }

// Calculate faces perpendicular to z (bottom and top).

    if (zoff == 0)
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index( i, j, 0 );

                Assert( Fhat(c,4) == 0. );
                Fhat( c, 4 ) = 
                    ( alpha_bottom*Fprime(c,4) -
                      2.*D(c,4)*f_b_bottom(i,j)/dz(0) ) /
                    ( alpha_bottom - 2.*beta_bottom*D(c,4)/dz(0) );
            }

    if (zoff+nczp == ncz)
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index( i, j, ncz-1 );

                Assert( Fhat(c,5) == 0. );
                Fhat( c, 5 ) = 
                    ( alpha_top*Fprime(c,5) -
                      2.*D(c,5)*f_b_top(i,j)/dz(ncz-1) ) /
                    ( alpha_top - 2.*beta_top*D(c,5)/dz(ncz-1) );
            }
#endif
}

//---------------------------------------------------------------------------//
// Calculate the coefficient matrix.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_A( const ccsf& sigmaabar )
{
    A = 0.;

// Note that A(c,3) is the diagonal.
//           A(c,2) is for interaction with cell (i-1,j,k)
//           A(c,4) is for interaction with cell (i+1,j,k)
//           A(c,1) is for interaction with cell (i,j-1,k)
//           A(c,5) is for interaction with cell (i,j+1,k)
//           A(c,0) is for interaction with cell (i,j,k-1)
//           A(c,6) is for interaction with cell (i,j,k+1)

    for( int c=0; c < ncp; c++ )
    {
        int i = I(c), j = J(c), k = K(c);

        A(c,3) = sigmaabar(c)*dx(i)*dy(j)*dz(k)
        // left face
            + dy(j)*dz(k) * Dhat(c,0) / dx(i)
        // right face
            + dy(j)*dz(k) * Dhat(c,1) / dx(i)
        // front face
            + dx(i)*dz(k) * Dhat(c,2) / dy(j)
        // back face
            + dx(i)*dz(k) * Dhat(c,3) / dy(j)
        // bottom face
            + dx(i)*dy(j) * Dhat(c,4) / dz(k)
        // top face
            + dx(i)*dy(j) * Dhat(c,5) / dz(k);

        if (i > 0)
            A(c,2) = -dy(j)*dz(k) * Dhat(c,0) / dx(i);
        if (i < ncx-1)
            A(c,4) = -dy(j)*dz(k) * Dhat(c,1) / dx(i);

        if (j > 0)
            A(c,1) = -dx(i)*dz(k) * Dhat(c,2) / dy(j);
        if (j < ncy-1)
            A(c,5) = -dx(i)*dz(k) * Dhat(c,3) / dy(j);

        if (k > 0)
            A(c,0) = -dx(i)*dy(j) * Dhat(c,4) / dz(k);
        if (k < ncz-1)
            A(c,6) = -dx(i)*dy(j) * Dhat(c,5) / dz(k);
    }
}

//---------------------------------------------------------------------------//
// Calculate the first term of b, involving only Qbar_r and V^c.  This is the 
// version used when solving the conduction equations, for which f_b, and
// therefore Fhat, are 0.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_b( const ccsf& Qbar_r )
{
    b = 0.;

    for( int c=0; c < ncp; c++ )
    {
        int i = I(c), j = J(c), k = K(c);
        b(c) = Qbar_r(c) * dx(i) * dy(j) * dz(k);
    }
}

//---------------------------------------------------------------------------//
// Calculate the full form of b, including the sum over faces.  This is the
// version used when solving full radiation diffusion.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_b( const ccsf& Qbar_r, const fcdsf& Fh )
{
    calculate_b( Qbar_r );

    for( int c=0; c < ncp; c++ )
    {
        int i = I(c), j = J(c), k = K(c);

        b(c) -= ( (Fhat(c,0) + Fhat(c,1)) * dy(j) * dz(k) +
                  (Fhat(c,2) + Fhat(c,3)) * dx(i) * dz(k) +
                  (Fhat(c,4) + Fhat(c,5)) * dx(i) * dy(j) );
    }
}

//---------------------------------------------------------------------------//
// This method actually solves the equations which have been built up in all
// the other methods.  Since this is a 7 banded SPD system, it is ideal for
// congugate gradient.  We invoke PCG for the actual solve.
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::solve_A_phi_equals_b( ccsf& phi )
{
    Assert( phi.size() == ncp );
    Assert( b.size() == ncp );

// Grr, have to make aliases since Mesh_XYZ::cell_array is no longer
// publicly derived from dsxx::Mat1.

    dsxx::Mat1<double> x( phi.begin(), ncp );
    dsxx::Mat1<double> rhs( b.begin(), ncp );

// Now solve the matrix equation A.x = rhs.

    pcg_ctrl.pcg_fe( x, rhs, spmv, precond );
}

//---------------------------------------------------------------------------//
// Calculate the updated values of F based on the results of our
// computation. 
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_P1<MT>::calculate_new_F( const fcdsf& D, const fcdsf& Fprime,
                                        const bssf& f_b, const ccsf& phi,
					fcdsf& F )
{
#ifdef RMR
    F = 0.;                     // should be unnecessary...

    {
    // Do the interior part.  Can do this with expression templates.

        fcdsf fphic( spm ), fphicprime( spm );
        fphic = phi;
        calculate_opposing_face( fphic, fphicprime );

        F = Ftwidle - Dtwidle*( fphicprime - fphic )/fdeltal;
    }

    const Mat2<double>& f_b_left   = f_b.face(0);
    const Mat2<double>& f_b_right  = f_b.face(1);
    const Mat2<double>& f_b_front  = f_b.face(2);
    const Mat2<double>& f_b_back   = f_b.face(3);
    const Mat2<double>& f_b_bottom = f_b.face(4);
    const Mat2<double>& f_b_top    = f_b.face(5);

// Calculate faces perpendicular to x (left and right).

    for( int k=zoff; k < zoff+nczp; k++ )
        for( int j=0; j < ncy; j++ )
        {
            int cl = local_cell_index( 0, j, k );
            int cr = local_cell_index( ncx-1, j, k );

            F(cl,0) =
                ( alpha_left*Fprime(cl,0) -
                  2.*D(cl,0)*( f_b_left(j,k) -
                               alpha_left * phi(cl) ) / dx(0) ) /
                ( alpha_left - 2.* beta_left * D(cl,0)/dx(0) );

            F(cr,1) =
                ( alpha_right*Fprime(cr,1) -
                  2.*D(cr,1)*( f_b_right(j,k) -
                               alpha_right * phi(cr) ) / dx(ncx-1) ) /
                ( alpha_right - 2.* beta_right * D(cr,1)/dx(ncx-1) );
        }

// Calculate faces perpendicular to y (front and back).

    for( int k=zoff; k < zoff+nczp; k++ )
        for( int i=0; i < ncx; i++ )
        {
            int cf = local_cell_index( i, 0, k );
            int cb = local_cell_index( i, ncy-1, k );

            F(cf,2) =
                ( alpha_front*Fprime(cf,2) -
                  2.*D(cf,2)*( f_b_front(i,k) -
                               alpha_front * phi(cf) ) / dy(0) ) /
                ( alpha_front - 2.* beta_front * D(cf,2)/dy(0) );

            F(cb,3) =
                ( alpha_back*Fprime(cb,3) -
                  2.*D(cb,3)*( f_b_back(i,k) -
                               alpha_back * phi(cb) ) / dy(ncy-1) ) /
                ( alpha_back - 2.* beta_back * D(cb,3)/dy(ncy-1) );
        }

// Calculate faces perpendicular to z (bottom and top).

    if (zoff == 0)
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index( i, j, 0 );

                F(c,4) =
                    ( alpha_bottom*Fprime(c,4) -
                      2.*D(c,4)*( f_b_bottom(i,j) -
                                  alpha_bottom * phi(c) ) / dz(0) ) /
                    ( alpha_bottom - 2.* beta_bottom * D(c,4)/dz(0) );
            }

    if (zoff+nczp == ncz)
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
            {
                int c = local_cell_index( i, j, ncz-1 );

                F(c,5) =
                    ( alpha_top*Fprime(c,5) -
                      2.*D(c,5)*( f_b_top(i,j) -
                                  alpha_top * phi(c) ) / dz(ncz-1) ) /
                    ( alpha_top - 2.* beta_top * D(c,5)/dz(ncz-1) );
            }
#endif
}

//---------------------------------------------------------------------------//
//                              end of Diffusion_P1.cc
//---------------------------------------------------------------------------//
