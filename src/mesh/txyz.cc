//----------------------------------*-C++-*----------------------------------//
// txyz.cc
// Geoffrey M. Furnish
// Wed May 13 09:53:44 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "mesh/Mesh_XYZ.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
using namespace std;

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    cout << "t1: passed\n";

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    SP<Mesh_XYZ> spm = new Mesh_XYZ( mdb );

    cout << "t2: passed" << endl;


    Mesh_XYZ::cell_array<double> x( spm ), y( spm ), z( spm );
    Mesh_XYZ::cell_array<int> xi( spm ), yi( spm ), zi( spm );
    Mesh_XYZ::fcdsf xf( spm ), yf( spm ), zf( spm );

    x = 1.;
    y = 2.;
    z = x + y;

    xi = 1;
    yi = 2;
    zi = xi + yi;

    xf = 1.;
    yf = xf;
    zf = xf + yf;
    zf = xf - yf;

    x = 1.;
    xf = x;

    Mesh_XYZ::guarded_cell_array<double> xgc( spm );

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
