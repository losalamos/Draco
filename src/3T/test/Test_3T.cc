//----------------------------------*-C++-*----------------------------------//
// Test_3T.cc
// Geoffrey M. Furnish
// Wed Nov 19 17:05:08 1997
//---------------------------------------------------------------------------//
// @> Test problem template.
//---------------------------------------------------------------------------//

#include "3T/test/Test_3T.hh"

#include <iostream.h>
#include <algorithm>

#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "c4/Baton.hh"

using namespace C4;

#include "util/dump.hh"
#include "util/ADFile.hh"

#include "linalg/pcg_DB.hh"
#include "linalg/PCG_MatVec.hh"
#include "linalg/PCG_PreCond.hh"
#include "linalg/PCG_Ctrl.hh"

#include "3T/MatVec_3T.hh"
#include "3T/PreCond.hh"

using dsxx::SP;

//---------------------------------------------------------------------------//
// Constructor.  Perform basic initialization of the test problem.
//---------------------------------------------------------------------------//

template<class MT, class Problem>
Test_3T<MT, Problem>::Test_3T( const SP<MT>& spm_,
			       const Run_DB& rdb,
                               const Diffusion_DB& diffdb,
			       const typename Problem::params& p,
			       const pcg_DB& pcg_db_,
                               int verbose_ )
    : Run_DB( rdb ), Problem(p),
      MT::Coord_Mapper( spm_->get_Mesh_DB() ),
      spm(spm_),
      Df( spm ),
      pcg_db(pcg_db_),
      verbose(verbose_)
{
    spd = new Diffusion_XYZ<MT>( diffdb, spm, pcg_db_ );

    if (node == 0)
	adf = new ADFile( "test.dat", IDX::WRITE, 500 );
    else
	adf = (ADFile *) 1;	// Don't /eeeeven/ ask.

// Now build up D.  D does not depend on time, so we can do this here.

    MT::fcdsf xf(spm);
    MT::fcdsf yf(spm);
    MT::fcdsf zf(spm);
    spm->get_xloc(xf);
    spm->get_yloc(yf);
    spm->get_zloc(zf);

    for( int f=0; f < 6; f++ )
	for( int c=0; c < ncp; c++ )
	{
	    Df(c,f) = D( xf(c,f), yf(c,f), zf(c,f) );
	}
}

//---------------------------------------------------------------------------//
// Clean up after the run is done.
//---------------------------------------------------------------------------//

template<class MT, class Problem>
Test_3T<MT, Problem>::~Test_3T()
{
    if (node == 0)
	delete adf;
}

//---------------------------------------------------------------------------//
// Run the simulation.  Contains t=0 setup, timestep cycling, and error
// analysis. 
//---------------------------------------------------------------------------//

template<class MT, class Problem>
void Test_3T<MT, Problem>::run()
{
    char buf[ 100 ];

    cout << "Running the test problem," << endl;

// Set up initial conditions of the problem.

    cell_array_double Eo(spm), En(spm);

    Mat1<double> xc( spm->get_xc() );
    Mat1<double> yc( spm->get_yc() );
    Mat1<double> zc( spm->get_zc() );

// Still need to bust out the solver, but for now, just suck in some things
// that were exported to the mesh class.

    Mat1<double> xf( spm->get_xf() );
    Mat1<double> yf( spm->get_yf() );
    Mat1<double> zf( spm->get_zf() );

    double dx = spm->get_dx(), dy = spm->get_dy(), dz = spm->get_dz();

    Mat1<double> xA( spm->get_xA() );
    Mat1<double> yA( spm->get_yA() );
    Mat1<double> zA( spm->get_zA() );

//    Mat1<double> vc( spm->get_vc() );
    Mat1<double> vc( ncp );
    std::copy( spm->get_vc().begin(), spm->get_vc().end(), vc.begin() );

    cout << "verbose = " << verbose << endl;

// Setup the "old" values for the first timestep.  This is just the analytic
// solution at t=0.

    for( int i=0; i < ncp; i++ ) {
        double x = xc( I(i) );
        double y = yc( J(i) );
        double z = zc( K(i) );

        Eo(i) = E(x,y,z,0.);

        if (verbose > 1) {
            sprintf( buf,
                     "node %d, I(%d)=%d xc(I(%d))=%lf J(%d)=%d yc(J(%d))=%lf",
                     node, i, I(i), i, xc(I(i)), i, J(i), i, yc(J(i)) );
            cout << buf << endl;
        }
    }

    cell_array_double E_analytic(spm), rhs(spm), r(spm);
    MT::fcdsf Eb(spm);

    bool test_passed = true;

// main loop.
    for( int ns=1; ns <= nsteps && test_passed; ns++ ) {

    // Calculate current time level.
	double t = dt * ns;

	if (node == 0) {
	    cout << "Begining timestep " << ns << " t=" << t << endl;
	    sprintf( buf, "Timestep %d", ns );
	    dump( adf, buf, "Timestep Notification" );
	}

    // Calculate analytic solution.
	for( int i=0; i < ncp; i++ )
	    E_analytic(i) = E( xc(I(i)), yc(J(i)), zc(K(i)), t );

// 	diag_out( E_analytic, "E_analytic" );

    // Calculate rhs(t).
	for( int i=0; i < ncp; i++ ) {
	    double x = xc( I(i) );
	    double y = yc( J(i) );
	    double z = zc( K(i) );
	    
	    rhs(i) = Exyz(x,y,z)*dEtdt()

		- (Dy(y)*Dz(z)*Et(t)*Ey(y)*Ez(z)) *
		( Dx(x)*d2Exdxx(x) + dExdx(x)*dDxdx() )

		- (Dx(x)*Dz(z)*Et(t)*Ex(x)*Ez(z)) *
		( Dy(y)*d2Eydyy(y) + dEydy(y)*dDydy() )

		- (Dx(x)*Dy(y)*Et(t)*Ex(x)*Ey(y)) *
		( Dz(z)*d2Ezdzz(z) + dEzdz(z)*dDzdz() );
	}

    // Calculate boundary values.

	Eb = 0.;
	for( int n=0; n < ncp; n++ ) {
	    int i = I(n), j = J(n), k = K(n);

	    if (i == 0)     Eb(n,0) = E( xc(i) - dx, yc(j), zc(k), t );
	    if (i == ncx-1) Eb(n,1) = E( xc(i) + dx, yc(j), zc(k), t );

	    if (j == 0)     Eb(n,2) = E( xc(i), yc(j) - dy, zc(k), t );
	    if (j == ncy-1) Eb(n,3) = E( xc(i), yc(j) + dy, zc(k), t );

	    if (k == 0)     Eb(n,4) = E( xc(i), yc(j), zc(k) - dz, t );
	    if (k == ncz-1) Eb(n,5) = E( xc(i), yc(j), zc(k) + dz, t );
	}

    // Begin setting up to solve the problem.  First, set r.
	for( int i=0; i < ncp; i++ )
	    r(i) = rhs(i)*dt + Eo(i);

    // Ask the solver to solve the system.

	spd->solve( Df, r, dt, En, Eb );

    // Grrr.  Make some more aliases

        Mat1<double> EE_analytic( &E_analytic(0), ncp );
        Mat1<double> EEn( &En(0), ncp );

    // Compute A . E_analytic = b_analytic

        Mat1<double> b_anal( ncp );
        b_anal = 0.;
//         spd->get_matvec()->MatVec( b_anal, E_analytic );
        spd->get_matvec()->MatVec( b_anal, EE_analytic );

    // Compute the product of A and En to see how well it compares to r.
	double r1 = 0.;
	int pcgits = spd->get_matvec()->get_iterations();

	Mat1<double> b( ncp );
// 	spd->get_matvec()->MatVec( b, En );
	spd->get_matvec()->MatVec( b, EEn );

	for( int i=0; i < ncp; i++ ) {
	    double d = b[i] - r[i];
	    r1 += d*d;
	}
        gsum(r1);
        r1 = sqrt( r1 / nct );

    // Compute solution quality metrics
	double s1=0.;
	for( int i=0; i < ncp; i++ ) {
	    double d = E_analytic[i] - En[i];
	    s1 += d*d;
	}
	gsum(s1);
	s1 = sqrt( s1 / nct );

	double s2=0.;
	for( int i=0; i < ncp; i++ ) {
	    double d = E_analytic[i] * E_analytic[i];
	    s2 += d;
	}
	gsum(s2);
	s2 = sqrt( s2 / nct );
	if (node == 0) {
	    cout << "Difference norm: s1=" << s1
		 << " s2=" << s2 << endl;
	    if (s2 > 0.) {
		cout << "s1/s2 = " << s1/s2 << endl;
		cout << "Also, r1=" << r1 << " r1/s2 = " << r1/s2 << endl;
	    }
	    cout << "PCG iterations = " << pcgits << endl;

            if (s2 == 0. && s1 > 0. ||
                s2 > 0. && s1/s2 > 1.e-9 )
                test_passed=false;            
	}

    // Prepare for next iteration
	Eo = En;
    }

    cout << "\n\ntest " << (test_passed ? "passed" : "failed") << "\n\n";
}

//---------------------------------------------------------------------------//
// This routine dumps a cell array out to the tty for visual inspection.
//---------------------------------------------------------------------------//

template<class MT, class Problem>
void Test_3T<MT, Problem>::diag_out( const cell_array_double& data,
				     const char *name ) const
{
    char buf[80];

    Mat1<double> xc( spm->get_xc() );
    Mat1<double> yc( spm->get_yc() );
    Mat1<double> zc( spm->get_zc() );

    if (node == 0) {
    // Dump own stuff.
	for( int i=0; i < ncp; i++ ) {
	    sprintf( buf, "node %d i=%d x=%lf y=%lf z=%lf %s=%lf",
		     node, i, xc(I(i)), yc(J(i)), zc(K(i)),
		     name, data(i) );
	    cout << buf << endl;
	}

    // Now collect from the other nodes.

	for( int n=1; n < nodes; n++ ) {
	    Send( 1, n ); // Tell them to transmit.
	    for( int i=0; i < ncp; i++ ) {
		Recv( buf, 80, n );
		cout << buf << endl;
	    }
	}
    }
    else {
	int dummy;
	Recv( dummy, 0 );

	for( int i=0; i < ncp; i++ ) {
	    sprintf( buf, "node %d i=%d x=%lf y=%lf z=%lf %s=%lf",
		     node, i, xc(I(i)), yc(J(i)), zc(K(i)),
		     name, data(i) );
	    Send( buf, 80, 0 );
	}
    }
}

//---------------------------------------------------------------------------//
//                              end of Test_3T.cc
//---------------------------------------------------------------------------//
