//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.cc
// RHS Linux User
// Thu Nov 20 22:08:34 1997
//---------------------------------------------------------------------------//
// @> A 3-d cartesian structured mesh facility.
//---------------------------------------------------------------------------//

#include "mesh/Mesh_XYZ.hh"
#include "c4/SpinLock.hh"
#include "c4/Baton.hh"
using namespace C4;

#include <iostream.h>

XYZ_Mapper::XYZ_Mapper( const Mesh_DB& mdb )
    : Mesh_DB( mdb )
{
    nct = ncx * ncy * ncz;
    nxy = ncx * ncy;

    Insist( ncz > nodes, "Current decomposition requires ncz > nodes." );

// Calculate # of cells in z on this processor.
    nczp = ncz / nodes + ( (ncz % nodes) > node );

    ncp = nxy * nczp;

// Calculate the starting z index for this processor's data.

    zoff = 0;
    {
        Baton<int> s(zoff);
        zoff = s;
        s += ncz;
    }

// Calculate the starting global cell index for this processor's data.

    goff = 0;
    {
	Baton<int> s(goff);
	goff = s;
	s += ncp;
    }
}


Mesh_XYZ::Mesh_XYZ( const Mesh_DB& mdb )
    : XYZ_Mapper( mdb ),

      vc( this ),

      xF( this ), yF( this ), zF( this )
{
    char buf[80];

    if (dump_indicies)
        for( int i=0; i < ncp; i++ ) {
            sprintf( buf, "node %d, i=%d, I(%d)=%d J(%d)=%d K(%d)=%d",
                     node, i, i, I(i), i, J(i), i, K(i) );
            cout << buf << endl;
        }

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
	xc(i) = xmin + (i+.5)*dx;
    for( int i=0; i < ncy; i++ )
	yc(i) = ymin + (i+.5)*dy;
    for( int i=0; i < ncz; i++ )
	zc(i) = zmin + (i+.5)*dz;

    if (dump_coords) {
        cout << "node " << node << " ncy=" << ncy << endl;
        for( int i=0; i < ncy; i++ )
            cout << "node " << node << " yc(" << i << ")=" << yc(i) << endl;
    }

// Initialize the face locations.

    for( int i=0; i < ncx+1; i++ )
	xf(i) = xmin + i * dx;
    for( int i=0; i < ncy+1; i++ )
	yf(i) = ymin + i * dy;
    for( int i=0; i < ncz+1; i++ )
	zf(i) = zmin + i * dz;

// Initialize the cell volumes we will be concerned with.

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

// Uhh, initialize the "new" face locations...

    for( int c=0; c < ncp; c++ )
    {
	int i = I(c), j = J(c), k = K(c);

    // Calculate cell center.

	double x = dx * (i + .5), dx2 = dx / 2.;
	double y = dy * (j + .5), dy2 = dy / 2.;
	double z = dz * (k + .5), dz2 = dz / 2.;

    // Now loop through all the faces, and initiailize face locations using
    // deltas from cell center.

    // 0 == left,  1 == right       x variation
    // 2 == front  3 == back        y variation
    // 4 == bottom 5 == top         z variation

	xF(c,0) = x - dx2; xF(c,1) = x + dx2;
	xF(c,2) = x; xF(c,3) = x;
	xF(c,4) = x; xF(c,5) = x;

	yF(c,0) = y; yF(c,1) = y;
	yF(c,2) = y - dy2; yF(c,3) = y + dy2;
	yF(c,4) = y; yF(c,5) = y;

	zF(c,0) = z; zF(c,1) = z;
	zF(c,2) = z; zF(c,3) = z;
	zF(c,4) = z - dz2; zF(c,5) = z + dz2;
    }

// Now initialize the diag offsets which are needed by clients employing
// sparse matrices.

    diags[3] = 0;
    diags[4] = 1;
    diags[5] = ncx;
    diags[6] = ncx*ncy;
    diags[2] = -diags[4];
    diags[1] = -diags[5];
    diags[0] = -diags[6];
}

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.cc
//---------------------------------------------------------------------------//
