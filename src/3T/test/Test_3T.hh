//----------------------------------*-C++-*----------------------------------//
// Test_3T.hh
// Geoffrey M. Furnish
// Wed Nov 19 17:05:08 1997
//---------------------------------------------------------------------------//
// @> Test problem template.
//---------------------------------------------------------------------------//

#ifndef __3T_test_Test_3T_hh__
#define __3T_test_Test_3T_hh__

#include "Test_Prob.hh"
#include "Run_DB.hh"

class ADFile;

#include "3T/Diffusion_XYZ.hh"
#include "linalg/pcg_DB.hh"

//===========================================================================//
// class Test_3T - Template class for accomodating various formulations

// This class implements the Test_Prob abstraction, but is not itself a
// "complete" test.  Rather, it is a template for a test, and the specific
// test problem is provided as a template parameter.
//===========================================================================//

template<class MT, class Problem>
class Test_3T : public Test_Prob,
		private Run_DB,
		private Problem,
		private MT::Coord_Mapper
{
    dsxx::SP<MT> spm;
    typename MT::fcdsf Df;

    dsxx::SP< Diffusion_XYZ<MT> > spd;

    pcg_DB pcg_db;

    ADFile *adf;

    int verbose;

  public:
    typedef double NumT;
    typedef typename MT::ccsf cell_array_double;

    Test_3T( const dsxx::SP<MT>& spm_, const Run_DB& rdb,
             const Diffusion_DB& diffdb, const typename Problem::params& p,
             const pcg_DB& pcg_db_, int verbose_ );
    ~Test_3T();

    void run();

    void diag_out( const cell_array_double& data,
		   const char *name ) const;

    int get_ncp() const { return ncp; }
    int get_nct() const { return nct; }
    int get_goff() const { return goff; }
};

#endif                          // __3T_test_Test_3T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/Test_3T.hh
//---------------------------------------------------------------------------//
