//----------------------------------*-C++-*----------------------------------//
// main.cc
// Geoffrey Furnish
// Thu Sep 11 11:31:43 1997
//---------------------------------------------------------------------------//
// @> Baby solver test case.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Okay, I'm a little frustrated.  The purpose of this test program is to
// just make sure I can set up and solve the equation, without having to use
// Pooma and all the associated machinery.  This will just be a bare bones
// test case.
//---------------------------------------------------------------------------//

#include <iostream.h>
#include <math.h>

#include "Run_DB.hh"
#include "Test_Prob.hh"

#include "ds++/Mat.hh"
#include "ds++/Algorithms.hh"

#include "c4/global.hh"

#include "util/ADFile.hh"
#include "util/dump.hh"
#include "util/g.hh"

void solve( const Mat2<double>& A, Mat1<double>& x, const Mat1<double>& b );

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    int node = C4::node();
    int nodes = C4::nodes();

    char buf[80];
    sprintf( buf, "node %d, nodes=%d\n", node, nodes );
    cout << buf << flush;

    if (node == 0)
	cout << "Initiating baby solver test." << endl;

    Test_Prob prob;
    prob.run();
#if 0
    ADFile *adf = new ADFile( "t1.dat", IDX::WRITE, 500 );

    NML_Group g( "t1" );

    rdb.setup_namelist( g );

    g.readgroup( "t1.in" );
    g.writegroup( "params.out" );

// Uhh

    int ncx = rdb.ncx;
    int coordsys = rdb.coordsys;
    int nsteps = rdb.nsteps;

    double xmin = rdb.xmin;
    double xmax = rdb.xmax;

    double dt = rdb.dt;

    cout << "ncx = " << ncx << endl;
    cout << "coordsys = " << coordsys << endl;
    cout << "nsteps = " << nsteps << endl;

    cout << "xmin = " << xmin << endl;
    cout << "xmax = " << xmax << endl;

    cout << "dt = " << dt << endl;

    double a = rdb.a;
    double b = rdb.b;
    double c = rdb.c;
    double d = rdb.d;
    double e = rdb.e;

// Calculate face locations, face areas, cell volumes, etc.

    double dx = (xmax - xmin) / ncx;

    Mat1<double> xf( ncx+1 );
    for( int i=0; i < ncx+1; i++ )
	xf[i] = i * dx;

// Assume coordsys = cartesian 1d for now.

    Mat1<double> af( ncx+1 );
    for( int i=0; i < ncx+1; i++ )
	af[i] = 1.;

    Mat1<double> E_analytic( ncx );
    Mat1<double> vc( ncx ), xc( ncx );
    for( int i=0; i < ncx; i++ ) {
	vc[i] = 1.;
	xc[i] = .5 * (xf[i] + xf[i+1]);
    }

// Set up initial conditions of the problem.

    Mat1<double> Eo( ncx ), En( ncx ), rhs( ncx ), r( ncx );;

// Set up E(x,0)
    for( int i=0; i < ncx; i++ ) {
	double x = xc[i];
	Eo(i) = a * (c*x*x + d*x + e);
    }

    Mat2<double> A( ncx, ncx );

// main loop

    for( int ns=1; ns <= nsteps; ns++ ) {

    // Calculate current time level.
	double t = dt * ns;

    // Calculate analytic solution.
	for( int i=0; i < ncx; i++ ) {
	    double x = xc[i];
	    E_analytic[i] = (a + b*t)*(c*x*x + d*x + e);
	}

    // Calculate rhs(t).
	rhs = b - 2*c;		// Note, not time dependent for the assumed
				// form of E(x,t).

    // Calculate left and right boundary values.

	double xl = xmin - .5*dx;
	double xr = xmax + .5*dx;
	double E_l = (a + b*t)*(c*xl*xl + d*xl + e);
	double E_r = (a + b*t)*(c*xr*xr + d*xr + e);

    // Begin setting up to solve the problem.  First, set r.
	for( int i=0; i < ncx; i++ )
	    r[i] = rhs[i]*(dt/vc[i]) + Eo[i];

    // Set up coefficient matrix.
	A = 0.;

    // Loop over interior rows.
	for( int row=1; row < ncx-1; row++ ) {
	    int i=row;		// uhh, why not make row => i?

	    A(i,i) = 1.;

	// Now add contributions from gradient on left face.
	    A(i,i-1) -= af(i)*(dt)/(dx*vc(i));// Recal A dot grad E, A is vector
	    A(i,i) += af(i)*(dt)/(dx*vc(i));

	// Now add contributions from gradient on right face.
	    A(i,i+1) -= af(i+1)*(dt)/(dx*vc(i));
	    A(i,i) += af(i+1)*(dt)/(dx*vc(i));
	}

    // Now handle the ends with boundary conditions.

    // Handle cell 0.
	int i=0;
	A(i,i) = 1.
	    + af(i)*(dt)/(dx*vc(i))
	    + af(i+1)*(dt)/(dx*vc(i));
	A(i,i+1) -= af(i+1)*(dt)/(dx*vc(i));
    // Now move A(i,i-1) over to rhs.
	r(i) += E_l*af(i)*(dt)/(dx*vc(i));

    // Handle cell ncx-1.
	i = ncx-1;
	A(i,i) = 1.
	    + af(i)*(dt)/(dx*vc(i))
	    + af(i+1)*(dt)/(dx*vc(i));
	A(i,i-1) -= af(i)*(dt)/(dx*vc(i));
    // Now move A(i,i+1) over to rhs.
	r(i) += E_r*af(i+1)*(dt)/(dx*vc(i));

	print_Mat( A, "A" );

    // Solve Ax = b
	solve( A, En, r );

    // Dump any data.

	cout << "Got En, here it is:" << endl;
	for( int i=0; i < ncx; i++ )
	    cout << En[i] << ' ';
	cout << endl;
	cout << "E_analytic is:" << endl;
	for( int i=0; i < ncx; i++ )
	    cout << E_analytic[i] << ' ';
	cout << endl;

	double s=0.;
	for( int i=0; i < ncx; i++ ) {
	    double d = E_analytic[i] - En[i];
	    s += d*d;
	}
	s = sqrt( s / ncx );
	if (s > 0.)
	    s /= norm( E_analytic.begin(), E_analytic.end() );
	cout << "Difference norm: " << s << endl;

    // Prepare for next iteration.
	Eo = En;
    }
#endif
    if (node == 0)
	cout << "Done with baby solver test." << endl;

    C4::Finalize();
}

void solve( const Mat2<double>& A, Mat1<double>& x, const Mat1<double>& b )
{
    int n = A.nx();

    cout << "Solving Ax = b   (NOT)" << endl;

    Mat2<double> Ainv( A );
    Mat1<double> binv( b );

    gaussj( Ainv, binv );

    x = 0.;
    for( int i=0; i < n; i++ )
	for( int j=0; j < n; j++ )
	    x[i] += Ainv(i,j)*b(j);

    cout << "x=" << endl;
    for( int i=0; i < n; i++ )
	cout << x[i] << ' ';
    cout << endl;

// Let's check that it solved right.
    Mat1<double> Ax( n );
    for( int i=0; i < n; i++ )
	for( int j=0; j < n; j++ )
	    Ax(i) += A(i,j)*x(j);

    cout << "Comparing Ax with b.\n";
    for( int i=0; i < n; i++ )
	cout << "Ax=" << Ax[i] << " b=" << b[i] << endl;

}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
