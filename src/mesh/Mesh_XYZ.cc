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
    ncp = nct / nodes + ( (nct % nodes) > node );

    nxy = ncx * ncy;

    goff = 0;
    {
	Baton<int> s(goff);
	goff = s;
	s += ncp;
    }

}


Mesh_XYZ::Mesh_XYZ( const Mesh_DB& mdb )
    : XYZ_Mapper( mdb ),
      xF( this ), yF( this ), zF( this )
{
    char buf[80];

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


//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.cc
//---------------------------------------------------------------------------//
