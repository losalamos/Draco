//----------------------------------*-C++-*----------------------------------//
// Test_Prob.cc
// Geoffrey Furnish
// Tue Sep 30 16:52:42 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "t3/Test_Prob.hh"
#include "t3/Run_DB.hh"
#include "t3/T3_MatVec.hh"
#include "t3/T3_PreCond.hh"

#include "nml/Group.hh"

#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "c4/Baton.hh"

using namespace C4;

#include "linalg/pcg_DB.hh"
#include "linalg/PCG_MatVec.hh"
#include "linalg/PCG_PreCond.hh"
#include "linalg/PCG_Ctrl.hh"

#include "ds++/SP.hh"

#include <math.h>

//---------------------------------------------------------------------------//
// Set up the supporting data structures.
//---------------------------------------------------------------------------//

Test_Prob::Test_Prob()
    : pcg_db("pcg")
{
// Set up namelist input.
    NML_Group g( "t3" );

    Run_DB::setup_namelist( g );
    pcg_db.setup_namelist( g );

    g. readgroup( "t3.in" );
    g.writegroup( "t3.out" );

// Now that we have the user input, we should be able to resize our data
// arrays, etc. 

    nct = ncx * ncy * ncz;
    ncp = nct / nodes + ( (nct % nodes) > node );

    SPINLOCK( cout << "node " << node << " has " << ncp
	      << " cells." << endl << flush );

    nxy = ncx * ncy;

    goff = 0;
    {
	Baton<int> s(goff);
	goff = s;
	s += ncp;
    }

    SPINLOCK( cout << "node " << node
	      << " has global offset " << goff << endl );
    C4::gsync();

// Allocate the arrays.

    xc = Mat1<double>( ncx );
    yc = Mat1<double>( ncy );
    zc = Mat1<double>( ncz );

    xf = Mat1<double>( ncx+1 );
    yf = Mat1<double>( ncy+1 );
    zf = Mat1<double>( ncz+1 );

// Initialize the deltas.

    dx = (xmax - xmin) / ncx;
    dy = (ymax - ymin) / ncy;
    dz = (zmax - zmin) / ncz;

// Initialize the cell center locations.

    for( int i=0; i < ncx; i++ )
	xc(i) = (i+.5)*dx;
    for( int i=0; i < ncy; i++ )
	yc(i) = (i+.5)*dy;
    for( int i=0; i < ncz; i++ )
	zc(i) = (i+.5)*dz;

// Initialize the face locations.

    for( int i=0; i < ncx+1; i++ )
	xf(i) = i * dx;
    for( int i=0; i < ncy+1; i++ )
	yf(i) = i * dy;
    for( int i=0; i < ncz+1; i++ )
	zf(i) = i * dz;

// Initialize the cell volumes we will be concerned with.

    vc = Mat1<double>( ncp );
    for( int i=0; i < ncp; i++ )
	vc(i) = dx * dy * dz;
}

//---------------------------------------------------------------------------//
// Perform the mock simulation.
//---------------------------------------------------------------------------//

void Test_Prob::run()
{
    C4::gsync();

    if (node == 0)
	cout << "Running the test problem, NOT!" << endl;

// Now we need a mapping between local cell id, and global i,j,k.  This will
// be needed in order to calculate the values of E and D.

// Set up initial conditions of the problem.

    Mat1<double> Eo(ncp), En(ncp);

    for( int i=0; i < ncp; i++ ) {
	double x = xc( I(i) );
	double y = yc( J(i) );
	double z = zc( K(i) );

	Eo(i) = E(x,y,z,0.);
    }

    Mat1<double> E_analytic(ncp), rhs(ncp), r(ncp);

    Mat2<double> A( ncp, nct );

// main loop.
    for( int ns=1; ns <= nsteps; ns++ ) {

    // Calculate current time level.
	double t = dt * ns;

	if (node == 0)
	    cout << "Begining timestep " << ns << " t=" << t << endl;

    // Calculate analytic solution.
	for( int i=0; i < ncp; i++ )
	    E_analytic(i) = E( xc(I(i)), yc(J(i)), zc(K(i)), t );

    // Calculate rhs(t).
	for( int i=0; i < ncp; i++ ) {
	    double x = xc( I(i) );
	    double y = yc( J(i) );
	    double z = zc( K(i) );
	    
	    rhs(i) = Exyz(x,y,z)*dEtdt()

		+ D(x,y,z)*Ey(y)*Ez(z)*d2Exdxx(x)
		+ Ey(y)*Ez(z)*dExdx(x)*Dy(y)*Dz(z)*dDxdx()

		+ D(x,y,z)*Ex(x)*Ez(z)*d2Eydyy(y)
		+ Ex(x)*Ez(z)*dEydy(y)*Dx(x)*Dz(z)*dDydy()

		+ D(x,y,z)*Ex(x)*Ey(y)*d2Ezdzz(z)
		+ Ex(x)*Ey(y)*dEzdz(z)*Dx(x)*Dy(y)*dDzdz();
	}

    // Calculate boundary values.

    // Begin setting up to solve the problem.  First, set r.
	for( int i=0; i < ncp; i++ )
	    r(i) = rhs(i)*dt/vc(i) + Eo(i);

    // Set up coefficient matrix.
	A = 0.;

	for( int i=0; i < ncp; i++ ) {
	    A(i,goff+i) = 1.;
	    
	}

    // Solve Ax = b.
	En = 0;

    // The best way to solve Ax = b.
	En = E_analytic;

    // A somewhat better approach for solving Ax = b.
	SP< PCG_MatVec<double> >  matvec  = new T3_MatVec<double>();
	SP< PCG_PreCond<double> > precond = new T3_PreCond<double>();

	PCG_Ctrl<double> pcg_ctrl( pcg_db, ncp );

	pcg_ctrl.pcg_fe( En, r, matvec, precond );

    // Compute solution quality metrics
	double s1=0.;
	for( int i=0; i < ncp; i++ ) {
	    double d = E_analytic[i] - En[i];
	    s1 += d*d;
	}
	gsum(s1);
	s1 = sqrt( s1 / nct );
	double s2=0.;
	for( int i=0; i < ncp; i++ )
	    s2 += E_analytic[i] * E_analytic[i];
	gsum(s2);
	s2 = sqrt( s2 / nct );
	if (s1 > 0.)
	    s1 /= s2;
	if (node == 0)
	    cout << "Difference norm: " << s1 << endl;

    // Prepare for next iteration
	Eo = En;
    }
}

//---------------------------------------------------------------------------//
//                              end of Test_Prob.cc
//---------------------------------------------------------------------------//
