//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.hh
// RHS Linux User
// Thu Nov 20 22:08:34 1997
//---------------------------------------------------------------------------//
// @> 
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

//    const Mesh_DB& get_Mesh_DB() const{ return *static_cast<Mesh_DB
//    *>(this); }
    const Mesh_DB& get_Mesh_DB() const{ return *this; }
};

//===========================================================================//
// class Mesh_XYZ - 

// 
//===========================================================================//

// class Mesh_XYZ : private Mesh_DB, private C4::NodeInfo
class Mesh_XYZ : private XYZ_Mapper//, private C4::NodeInfo
{
#if 0
    int nct;			// # of total cells in problem.
    int ncp;			// # of cells on this processor.

    int goff;			// global offset of this processor's cells.
    int nxy;

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
#endif

    Mat2<double> A;
    Mat1<double> xc, yc, zc;
    Mat1<double> xf, yf, zf;
    double       dx, dy, dz;
    Mat1<double> vc;
    Mat1<double> xA, yA, zA;

  public:
    Mesh_XYZ( const Mesh_DB& mdb );
//     Mesh_XYZ( const Mesh_XYZ& );
//     ~Mesh_XYZ();
//     Mesh_XYZ& operator=( const Mesh_XYZ& );

    typedef XYZ_Mapper Coord_Mapper;

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
