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

#include "util/dump.hh"
#include "util/ADFile.hh"

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
    char buf[80];

// Set up namelist input.
    NML_Group g( "t3" );

    Run_DB::setup_namelist( g );
    pcg_db.setup_namelist( g );

    g. readgroup( "t3.in" );
    g.writegroup( "t3.out" );

    if (node == 0)
	adf = new ADFile( "t3.dat", IDX::WRITE, 500 );
    else
	adf = (ADFile *) 1;	// Don't /eeeeven/ ask.


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

    for( int i=0; i < ncp; i++ ) {
	sprintf( buf, "node %d, i=%d, I(%d)=%d J(%d)=%d K(%d)=%d",
		 node, i, i, I(i), i, J(i), i, K(i) );
	cout << buf << endl;
    }

// Allocate the arrays.

    A = Mat2<double>( ncp, nct );

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
	xc(i) = xmin + (i+.5)*dx;
    for( int i=0; i < ncy; i++ )
	yc(i) = ymin + (i+.5)*dy;
    for( int i=0; i < ncz; i++ )
	zc(i) = zmin + (i+.5)*dz;

    cout << "node " << node << " ncy=" << ncy << endl;
    for( int i=0; i < ncy; i++ )
	cout << "node " << node << " yc(" << i << ")=" << yc(i) << endl;

// Initialize the face locations.

    for( int i=0; i < ncx+1; i++ )
	xf(i) = xmin + i * dx;
    for( int i=0; i < ncy+1; i++ )
	yf(i) = ymin + i * dy;
    for( int i=0; i < ncz+1; i++ )
	zf(i) = zmin + i * dz;

// Initialize the cell volumes we will be concerned with.

    vc = Mat1<double>( ncp );
    for( int i=0; i < ncp; i++ )
	vc(i) = dx * dy * dz;

// Initialize the face areas.

    xA = Mat1<double>( ncx+1 );
    yA = Mat1<double>( ncy+1 );
    zA = Mat1<double>( ncz+1 );

    for( int i=0; i < ncx+1; i++ )
	xA(i) = dy * dz;
    for( int i=0; i < ncy+1; i++ )
	yA(i) = dx * dz;
    for( int i=0; i < ncz+1; i++ )
	zA(i) = dx * dy;
}

Test_Prob::~Test_Prob()
{
    if (node == 0)
	delete adf;
}

//---------------------------------------------------------------------------//
// Perform the mock simulation.
//---------------------------------------------------------------------------//

void Test_Prob::run()
{
    char buf[80];

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

	sprintf( buf,
		 "node %d, I(%d)=%d xc(I(%d))=%lf J(%d)=%d yc(J(%d))=%lf",
		 node, i, I(i), i, xc(I(i)), i, J(i), i, yc(J(i)) );
	cout << buf << endl;
    }

    Mat1<double> E_analytic(ncp), rhs(ncp), r(ncp);

// main loop.
    for( int ns=1; ns <= nsteps; ns++ ) {

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

	{
	    char buf[80];

	    if (node == 0) {
	    // Dump own stuff.
		for( int i=0; i < ncp; i++ ) {
		    sprintf( buf, "node %d i=%d x=%lf y=%lf z=%lf E_analytic=%lf",
			     node, i, xc(I(i)), yc(J(i)), zc(K(i)),
			     E_analytic(i) );
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
		    sprintf( buf, "node %d i=%d x=%lf y=%lf z=%lf E_analytic=%lf",
			     node, i, xc(I(i)), yc(J(i)), zc(K(i)),
			     E_analytic(i) );
		    Send( buf, 80, 0 );
		}
	    }
	}

    // Calculate rhs(t).
	for( int i=0; i < ncp; i++ ) {
	    double x = xc( I(i) );
	    double y = yc( J(i) );
	    double z = zc( K(i) );
	    
	    rhs(i) = Exyz(x,y,z)*dEtdt()

// 		+ D(x,y,z)*Ey(y)*Ez(z)*d2Exdxx(x)
// 		+ Ey(y)*Ez(z)*dExdx(x)*Dy(y)*Dz(z)*dDxdx()

		- (Dy(y)*Dz(z)*Et(t)*Ey(y)*Ez(z)) *
		( Dx(x)*d2Exdxx(x) + dExdx(x)*dDxdx() )

// 		+ D(x,y,z)*Ex(x)*Ez(z)*d2Eydyy(y)
// 		+ Ex(x)*Ez(z)*dEydy(y)*Dx(x)*Dz(z)*dDydy()

		- (Dx(x)*Dz(z)*Et(t)*Ex(x)*Ez(z)) *
		( Dy(y)*d2Eydyy(y) + dEydy(y)*dDydy() )

// 		+ D(x,y,z)*Ex(x)*Ey(y)*d2Ezdzz(z)
// 		+ Ex(x)*Ey(y)*dEzdz(z)*Dx(x)*Dy(y)*dDzdz();

		- (Dx(x)*Dy(y)*Et(t)*Ex(x)*Ey(y)) *
		( Dz(z)*d2Ezdzz(z) + dEzdz(z)*dDzdz() );
	}

    // Calculate boundary values.

    // Begin setting up to solve the problem.  First, set r.
	for( int i=0; i < ncp; i++ )
	//	    r(i) = rhs(i)*dt/vc(i) + Eo(i);
	    r(i) = rhs(i)*dt + Eo(i);

    // Set up coefficient matrix.
	A = 0.;

	for( int n=0; n < ncp; n++ ) {
	    int i = I(n), j = J(n), k = K(n);
	    double d, fac;

	    A(n,goff+n) = 1.;

	// left face.

	    d = D( xf(i), yc(j), zc(k) );
	    fac = (d * dt) / (vc(n) * dx);
	    A(n,goff+n) += (fac * xA(i));
	    if (i > 0)
		A(n, goffset(i-1,j,k)) -= fac * xA(i);
	    else
		r(n) += fac * xA(i) * E( xc(i) - dx, yc(j), zc(k), t );

	// right face.

	    d = D( xf(i+1), yc(j), zc(k) );
	    fac = (d * dt) / (vc(n) * dx);
	    A(n,goff+n) += (fac * xA(i+1));
	    if (i < ncx-1)
		A(n, goffset(i+1,j,k)) -= fac * xA(i+1);
	    else
		r(n) += fac * xA(i+1) * E( xc(i)+dx, yc(j), zc(k), t );

	// front face.

	    d = D( xc(i), yf(j), zc(k) );
	    fac = (d * dt) / (vc(n) * dy);
	    A(n,goff+n) += (fac * yA(j));
	    if (j > 0)
		A(n, goffset(i,j-1,k)) -= fac * yA(j);
	    else
		r(n) += fac * yA(j) * E( xc(i), yc(j) - dy, zc(k), t );

	// back face.

	    d = D( xc(i), yf(j+1), zc(k) );
	    fac = (d * dt) / (vc(n) * dy);
	    A(n,goff+n) += (fac * yA(j+1));
	    if (j < ncy-1)
		A(n, goffset(i,j+1,k)) -= fac * yA(j+1);
	    else
		r(n) += fac * yA(j+1) * E( xc(i), yc(j)+dy, zc(k), t );

	// bottom face.

	    d = D( xc(i), yc(j), zf(k) );
	    fac = (d * dt) / (vc(n) * dz);
	    A(n,goff+n) += (fac * zA(k));
	    if (k > 0)
		A(n, goffset(i,j,k-1)) -= fac * zA(k);
	    else
		r(n) += fac * zA(k) * E( xc(i), yc(j), zc(k) - dz, t );

	// top face.

	    d = D( xc(i), yc(j), zf(k+1) );
	    fac = (d * dt) / (vc(n) * dz);
	    A(n,goff+n) += (fac * zA(k+1));
	    if (k < ncz-1)
		A(n, goffset(i,j,k+1)) -= fac * zA(k+1);
	    else
		r(n) += fac * zA(k+1) * E( xc(i), yc(j), zc(k)+dz, t );

	// Here comes a pig face.    (With appologies to Dr. Seus :-).
	}

	print2_Mat( A, "A" );

    // Solve Ax = b
// 	En = 0;

    // The best way to solve Ax = b.
// 	En = E_analytic;

    // A somewhat better approach for solving Ax = b.
	SP< PCG_MatVec<double> >  matvec  = new T3_MatVec<double>( this );
	SP< T3_MatVec<double> >  t3matvec  = matvec;
	SP< PCG_PreCond<double> > precond = new T3_PreCond<double>();

	PCG_Ctrl<double> pcg_ctrl( pcg_db, ncp );

	pcg_ctrl.pcg_fe( En, r, matvec, precond );

	int pcgits = t3matvec->get_iterations();

    // Calculate A.En, compare to r.

	Mat1<double> b( ncp );
	matvec->MatVec( b, En );

	double r1=0.;
	for( int i=0; i < ncp; i++ ) {
	    double d = b[i] - r[i];
	    r1 += d*d;
	}

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
// 	if (s2 > 0.)
// 	    s1 /= s2;
	if (node == 0) {
	    cout << "Difference norm: s1=" << s1
		 << " s2=" << s2 << endl;
	    if (s2 > 0.) {
		cout << "s1/s2 = " << s1/s2 << endl;
		cout << "Also, r1=" << r1 << " r1/s2 = " << r1/s2 << endl;
	    }
	    cout << "PCG iterations = " << pcgits << endl;
	}

    // Prepare for next iteration
	Eo = En;
    }
}

//---------------------------------------------------------------------------//
//                              end of Test_Prob.cc
//---------------------------------------------------------------------------//
