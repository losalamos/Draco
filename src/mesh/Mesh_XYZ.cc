//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.cc
// RHS Linux User
// Thu Nov 20 22:08:34 1997
//---------------------------------------------------------------------------//
// @> A 3-d cartesian structured mesh facility.
//---------------------------------------------------------------------------//

#include "Mesh_XYZ.hh"
#include "c4/SpinLock.hh"
#include "c4/Baton.hh"
#include "Mesh_DB.hh"
using namespace C4;

#include <iostream.h>

XYZ_Mapper::XYZ_Mapper( const Mesh_DB& mdb_in )
    : ncx(mdb_in.ncx), ncy(mdb_in.ncy), ncz(mdb_in.ncz), mdb(mdb_in)
{
    Insist(mdb.isValid(), "Invalid mdb used to initialize XYZ_Mapper.");
    
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
        s += nczp;
    }

// Calculate the starting global cell index for this processor's data.

    goff = 0;
    {
	Baton<int> s(goff);
	goff = s;
	s += ncp;
    }
}


Mesh_XYZ::Mesh_XYZ( const Mesh_DB& mdb_in )
    : XYZ_Mapper( mdb_in ),
      vc( this ),
      dX( this ), dY( this ), dZ( this ),
      xC( this ), yC( this ), zC( this ),
      xF( this ), yF( this ), zF( this ),
      face_norms( this )
{
// Allocate the arrays.

    xc = rtt_dsxx::Mat1<double>( ncx );
    yc = rtt_dsxx::Mat1<double>( ncy );
    zc = rtt_dsxx::Mat1<double>( ncz );

    xf = rtt_dsxx::Mat1<double>( ncx+1 );
    yf = rtt_dsxx::Mat1<double>( ncy+1 );
    zf = rtt_dsxx::Mat1<double>( ncz+1 );


// Initialize the cell center locations.

    xc(0) = mdb.xmin;
    yc(0) = mdb.ymin;
    zc(0) = mdb.zmin;

    for( int i=1; i < ncx; i++ )
	xc(i) = xc(i-1) + (mdb.dx[i]+mdb.dx[i-1])/2.0;
    for( int i=1; i < ncy; i++ )
	yc(i) = yc(i-1) + (mdb.dy[i]+mdb.dy[i-1])/2.0;
    for( int i=1; i < ncz; i++ )
	zc(i) = zc(i-1) + (mdb.dz[i]+mdb.dz[i-1])/2.0;

// Initialize the face locations.


    xf(0) = mdb.xmin;
    yf(0) = mdb.ymin;
    zf(0) = mdb.zmin;

    for( int i=1; i < ncx+1; i++ )
	xf(i) = xf(i-1) + mdb.dx[i-1];
    for( int i=1; i < ncy+1; i++ )
	yf(i) = yf(i-1) + mdb.dy[i-1];
    for( int i=1; i < ncz+1; i++ )
	zf(i) = zf(i-1) + mdb.dz[i-1];

// Initialize the face areas.

    //    xA = rtt_dsxx::Mat1<double>( ncx+1 );
    //    yA = rtt_dsxx::Mat1<double>( ncy+1 );
    //    zA = rtt_dsxx::Mat1<double>( ncz+1 );

    //   for( int i=0; i < ncx+1; i++ )
    //	xA(i) = dy * dz;
    //   for( int i=0; i < ncy+1; i++ )
    //	yA(i) = dx * dz;
    //    for( int i=0; i < ncz+1; i++ )
    //	zA(i) = dx * dy;

// Uhh, initialize the "new" face locations, cell centers, deltas...

    for( int c=0; c < ncp; c++ )
    {
	int i = I(c), j = J(c), k = K(c);

	const double dx = mdb.dx[i];
	const double dy = mdb.dy[j];
	const double dz = mdb.dz[k];
	const double x = xc[i];
	const double y = yc[j];
	const double z = zc[k];

    // Initialize the cell deltas.

        dX(c) = dx;
        dY(c) = dy;
        dZ(c) = dz;

    // Intitalize the cell volumes.

	vc(c) = dX(c) * dY(c) * dZ(c);

    // Initialize the cell centers.

        xC(c) = x;
        yC(c) = y;
        zC(c) = z;

    // Now loop through all the faces, and initialize face locations using
    // deltas from cell center.

    // 0 == left,  1 == right       x variation
    // 2 == front  3 == back        y variation
    // 4 == bottom 5 == top         z variation

	xF(c,0) = x - dx/2.0;   xF(c,1) = x + dx/2.0;
	xF(c,2) = x;            xF(c,3) = x;
	xF(c,4) = x;            xF(c,5) = x;

	yF(c,0) = y;            yF(c,1) = y;
	yF(c,2) = y - dy/2.0;   yF(c,3) = y + dy/2.0;
	yF(c,4) = y;            yF(c,5) = y;

	zF(c,0) = z;            zF(c,1) = z;
	zF(c,2) = z;            zF(c,3) = z;
	zF(c,4) = z - dz/2.0;   zF(c,5) = z + dz/2.0;
    }

// Initialize the unit normals

    xhat(0) = 1.;
    xhat(1) = 0.;
    xhat(2) = 0.;
    yhat(0) = 0.;
    yhat(1) = 1.;
    yhat(2) = 0.;
    zhat(0) = 0.;
    zhat(1) = 0.;
    zhat(2) = 1.;

// Initialize the face normals

    for( int c=0; c < ncp; c++ ) {
        face_norms(c,0) = -xhat;
        face_norms(c,1) = xhat;
        face_norms(c,2) = -yhat;
        face_norms(c,3) = yhat;
        face_norms(c,4) = -zhat;
        face_norms(c,5) = zhat;
    }

}

// do this like the xF, yF and zF

void Mesh_XYZ::get_face_areas(Mesh_XYZ::fcdsf& fa) const
{
    for( int c=0; c < ncp; c++ )
    {
	int i = I(c), j = J(c), k = K(c);

	const double dx = mdb.dx[i];
	const double dy = mdb.dy[j];
	const double dz = mdb.dz[k];

	fa(i,j,k,0) = dy * dz;
	fa(i,j,k,1) = fa(i,j,k,0);
	fa(i,j,k,2) = dx * dz;
	fa(i,j,k,3) = fa(i,j,k,2);
	fa(i,j,k,4) = dx * dy;
	fa(i,j,k,5) = fa(i,j,k,4);
    }
}

void Mesh_XYZ::get_face_lengths(fcdsf& fl) const
{
    for( int c=0; c < ncp; c++ )
    {
	int i = I(c), j = J(c), k = K(c);

	const double dx = mdb.dx[i];
	const double dy = mdb.dy[j];
	const double dz = mdb.dz[k];

	fl(i,j,k,0) = dx;
	fl(i,j,k,1) = fl(i,j,k,0);
	fl(i,j,k,2) = dy;
	fl(i,j,k,3) = fl(i,j,k,2);
	fl(i,j,k,4) = dz;
	fl(i,j,k,5) = fl(i,j,k,4);
    }
}

void Mesh_XYZ::get_vertex_volumes(vcsf& vols) const
{
    gather( vols, vc, Mesh_XYZ::OpAssign() );
    vols /= 8.;
}

void Mesh_XYZ::get_node_volumes(ncsf& vols) const
{
    vcsf vertex_volumes(this);
    get_vertex_volumes(vertex_volumes);
    scatter( vols, vertex_volumes, Mesh_XYZ::OpAddAssign() );
}

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.cc
//---------------------------------------------------------------------------//
