//----------------------------------*-C++-*----------------------------------//
// Diffusion_XYZ.cc
// Geoffrey Furnish
// Thu Dec  4 09:30:51 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/Diffusion_XYZ.hh"

#include "util/dump.hh"
#include "util/ADFile.hh"

//---------------------------------------------------------------------------//
// Initialize the diffusion solver.
//---------------------------------------------------------------------------//

template<class MT>
Diffusion_XYZ<MT>::Diffusion_XYZ( const SP<MT>& spm_, const pcg_DB& pcg_db )
    : MT::Coord_Mapper( spm_->get_Mesh_DB() ),
      spm(spm_),
      A( ncp, nct ),
      AA( static_cast<MT::Coord_Mapper>(*this), spm->get_diag_offsets() ),
      pcg_ctrl( pcg_db, ncp ),

      vc( spm->get_vc() ),

      dx( spm->get_dx() ), dy( spm->get_dy() ), dz( spm->get_dz() )
{
    spmv = new MatVec_3T< MT, Diffusion_XYZ<MT> >( spm, this );
    precond = new PreCond<double>();
}

//---------------------------------------------------------------------------//
// Solve the diffusion equation.  
//---------------------------------------------------------------------------//

template<class MT>
void Diffusion_XYZ<MT>::solve( const typename MT::fcdsf& D,
			       typename  MT::cell_array<double>& r,
			       double dt, typename MT::cell_array<double>& x,
			       const typename MT::fcdsf& Eb )
{
    const Mat1<double>& xA = spm->get_xA();
    const Mat1<double>& yA = spm->get_yA();
    const Mat1<double>& zA = spm->get_zA();

// Clear the matrix.

    A = 0.;
    AA = 0.;

// Now initialize the matrix.

    for( int n=0; n < ncp; n++ ) {
	int i = I(n), j = J(n), k = K(n);
	double d, fac;

	A(n,goff+n) = 1.;
        AA(n,3) = 1.;

    // left face.

    //	d = D( xf(i), yc(j), zc(k) );
	d = D( n, 0 );
	fac = (d * dt) / (vc(n) * dx);
	A(n,goff+n) += (fac * xA(i));
	AA(n,3) += (fac * xA(i));
	if (i > 0) {
	    A(n, goffset(i-1,j,k)) -= fac * xA(i);
	    AA(n,2) -= fac * xA(i);
	} else
// 	    r(n) += fac * xA(i) * E( xc(i) - dx, yc(j), zc(k), t );
	    r(n) += fac * xA(i) * Eb( n, 0 );

    // right face.

// 	d = D( xf(i+1), yc(j), zc(k) );
	d = D( n, 1 );
	fac = (d * dt) / (vc(n) * dx);
	A(n,goff+n) += (fac * xA(i+1));
	AA(n,3) += (fac * xA(i+1));
	if (i < ncx-1) {
	    A(n, goffset(i+1,j,k)) -= fac * xA(i+1);
	    AA(n,4) -= fac * xA(i+1);
	} else
// 	    r(n) += fac * xA(i+1) * E( xc(i)+dx, yc(j), zc(k), t );
	    r(n) += fac * xA(i+1) * Eb( n, 1 );

    // front face.

// 	d = D( xc(i), yf(j), zc(k) );
	d = D( n, 2 );
	fac = (d * dt) / (vc(n) * dy);
	A(n,goff+n) += (fac * yA(j));
	AA(n,3) += (fac * yA(j));
	if (j > 0) {
	    A(n, goffset(i,j-1,k)) -= fac * yA(j);
	    AA(n,1) -= fac * yA(j);
	} else
// 	    r(n) += fac * yA(j) * E( xc(i), yc(j) - dy, zc(k), t );
	    r(n) += fac * yA(j) * Eb( n, 2 );

    // back face.

// 	d = D( xc(i), yf(j+1), zc(k) );
	d = D( n, 3 );
	fac = (d * dt) / (vc(n) * dy);
	A(n,goff+n) += (fac * yA(j+1));
	AA(n,3) += (fac * yA(j+1));
	if (j < ncy-1) {
	    A(n, goffset(i,j+1,k)) -= fac * yA(j+1);
	    AA(n,5) -= fac * yA(j+1);
	} else
// 	    r(n) += fac * yA(j+1) * E( xc(i), yc(j)+dy, zc(k), t );
	    r(n) += fac * yA(j+1) * Eb( n, 3 );

    // bottom face.

// 	d = D( xc(i), yc(j), zf(k) );
	d = D( n, 4 );
	fac = (d * dt) / (vc(n) * dz);
	A(n,goff+n) += (fac * zA(k));
	AA(n,3) += (fac * zA(k));
	if (k > 0) {
	    A(n, goffset(i,j,k-1)) -= fac * zA(k);
	    AA(n,0) -= fac * zA(k);
	} else
// 	    r(n) += fac * zA(k) * E( xc(i), yc(j), zc(k) - dz, t );
	    r(n) += fac * zA(k) * Eb( n, 4 );

    // top face.

// 	d = D( xc(i), yc(j), zf(k+1) );
	d = D( n, 5 );
	fac = (d * dt) / (vc(n) * dz);
	A(n,goff+n) += (fac * zA(k+1));
	AA(n,3) += (fac * zA(k+1));
	if (k < ncz-1) {
	    A(n, goffset(i,j,k+1)) -= fac * zA(k+1);
	    AA(n,6) -= fac * zA(k+1);
	} else
// 	    r(n) += fac * zA(k+1) * E( xc(i), yc(j), zc(k)+dz, t );
	    r(n) += fac * zA(k+1) * Eb( n, 5 );

    // Here comes a pig face.    (With appologies to Dr. Seus :-).
    }

//     print2_Mat( A, "A" );

// Grr, have to make aliases since Mesh_XYZ::cell_array is not longer
// publicly derived from dsxx::Mat1.

    Mat1<double> xx( &x[0], ncp );
    Mat1<double> rr( &r[0], ncp );


// Now solve the matrix equation A.x = rhs.

//     pcg_ctrl.pcg_fe( x, r, spmv, precond );
    pcg_ctrl.pcg_fe( xx, rr, spmv, precond );

    int pcgits = spmv->get_iterations();
}

//---------------------------------------------------------------------------//
//                              end of Diffusion_XYZ.cc
//---------------------------------------------------------------------------//
