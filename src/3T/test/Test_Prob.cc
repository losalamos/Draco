//----------------------------------*-C++-*----------------------------------//
// Test_Prob.cc
// Geoffrey M. Furnish
// Wed Nov 19 16:18:54 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/test/Test_Prob.hh"
#include "3T/test/Test_3T.hh"
#include "3T/test/Run_DB.hh"
#include "3T/test/XYZ_Quadratic.hh"

#include "mesh/Mesh_XYZ.hh"

#include "nml/Group.hh"

SP<Test_Prob> Test_Prob_allocator( int argc, char *argv[] )
{
    SP<Test_Prob> prob;

    NML_Group g( "test" );

    Run_DB rdb;
    rdb.setup_namelist( g );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    Quad_Params qpdb;
    qpdb.setup_namelist( g );

    g.readgroup( "test.in" );
    g.writegroup( "test.out" );

    SP<Mesh_XYZ> spm = new Mesh_XYZ( mdb );

// Theoretically we could parse argc, argv to figure out which test problem
// to initiate.  For now, however, we just hardwire one.

    prob = new Test_3T< Mesh_XYZ, XYZ_Quadratic >( spm, qpdb );

    return prob;
}

//---------------------------------------------------------------------------//
//                              end of Test_Prob.cc
//---------------------------------------------------------------------------//
