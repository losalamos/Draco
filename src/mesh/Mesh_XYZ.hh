//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.hh
// RHS Linux User
// Thu Nov 20 22:08:34 1997
//---------------------------------------------------------------------------//
// @> A 3-d cartesian structured mesh facility.
//---------------------------------------------------------------------------//

#ifndef __mesh_Mesh_XYZ_hh__
#define __mesh_Mesh_XYZ_hh__

#include "Mesh_DB.hh"
#include "ds++/Mat.hh"
#include "ds++/SP.hh"
#include "c4/NodeInfo.hh"

struct XYZ_Mapper : public Mesh_DB, public C4::NodeInfo
{
    int nct;			// # of total cells in problem.
    int ncp;			// # of cells on this processor.

    int goff;			// global offset of this processor's cells.
    int nxy;

// Methods

    XYZ_Mapper( const Mesh_DB& mdb );

// Convert a local cell index to its (i,j,k) indexes in whole domain.

    int I(int x) const { x += goff; return (x % nxy) % ncx; }
    int J(int x) const { x += goff; return (x % nxy) / ncx; }
    int K(int x) const { x += goff; return  x / nxy; }

    int goffset( int i, int j, int k ) const { return ncx*(k*ncy + j) +i; }
    int goffset( int nc ) const
    {
	nc += goff;
	return goffset( I(nc), J(nc), K(nc) );
    }

    int local_cell_index( int i, int j, int k ) const
    {
	return goffset(i,j,k) - goff;
    }

    const Mesh_DB& get_Mesh_DB() const{ return *this; }
};

//===========================================================================//
// class Mesh_XYZ - A 3-d cartesian mesh class

// This is a 3-d cartesisan structured mesh.  The main purpose of having this
// class is in order for it to be instantiated by various test articles
// throughout the Draco system.  It could also be useful to Draco clients as
// a point of reference during the construction of their own mesh classes.
//===========================================================================//

class Mesh_XYZ : private XYZ_Mapper
{
  public:

// Face centered discontinuous scalar field
// Has a value on each face in each cell.
    
    class fcdsf : private XYZ_Mapper {
	friend class Mesh_XYZ;

	Mat2<double> data;

	fcdsf( const Mesh_XYZ *m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncp(), 6 )
	{}

      public:
	fcdsf( const SP<Mesh_XYZ>& m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncp(), 6 )
	{}

    // i, j, k == global xyz cell indicies
    // f == face index
    // c == local cell index

	double& operator()( int c, int f )       { return data(c,f); }
	double  operator()( int c, int f ) const { return data(c,f); }

	double& operator()( int i, int j, int k, int f )
	{
	    return data( local_cell_index(i,j,k), f );
	}
	double  operator()( int i, int j, int k, int f ) const 
	{
	    return data( local_cell_index(i,j,k), f );
	}
    };

  private:
    Mat1<double> xc, yc, zc;
    Mat1<double> xf, yf, zf;
    double       dx, dy, dz;
    Mat1<double> vc;
    Mat1<double> xA, yA, zA;

    fcdsf xF, yF, zF;

  public:
    typedef XYZ_Mapper Coord_Mapper;

    Mesh_XYZ( const Mesh_DB& mdb );

    const Mesh_DB& get_Mesh_DB() const { return XYZ_Mapper::get_Mesh_DB(); }
    const Mat1<double>& get_xc() const { return xc; }
    const Mat1<double>& get_yc() const { return yc; }
    const Mat1<double>& get_zc() const { return zc; }
    int get_ncp() const { return ncp; }

    double get_dx() const { return dx; }
    double get_dy() const { return dy; }
    double get_dz() const { return dz; }

    const Mat1<double>& get_xf() const { return xf; }
    const Mat1<double>& get_yf() const { return yf; }
    const Mat1<double>& get_zf() const { return zf; }

    const Mat1<double>& get_xA() const { return xA; }
    const Mat1<double>& get_yA() const { return yA; }
    const Mat1<double>& get_zA() const { return zA; }

    const Mat1<double>& get_vc() const { return vc; }

    class cell_array : public Mat1<double> {
      public:
	cell_array( const SP<Mesh_XYZ>& m ) : Mat1<double>( m->get_ncp() ) {}
    };
};

#endif                          // __mesh_Mesh_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Mesh_XYZ.hh
//---------------------------------------------------------------------------//
