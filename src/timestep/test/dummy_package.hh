//----------------------------------*-C++-*----------------------------------//
// dummy_package.hh
// John McGhee
// Thu Aug 27 07:48:41 1998
//---------------------------------------------------------------------------//
// @> A dummy package to exercize the field time-step advisors.
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __timestep_test_dummy_package_hh__
#define __timestep_test_dummy_package_hh__

#include "ds++/SP.hh"

// FORWARD REFERENCES

namespace rtt_timestep {

class ts_manager;
class field_ts_advisor;

}

//===========================================================================//
// class dummy_package - Exercizes the field time-step advisors.
//
// This class serves as an example of how any particular package can
// make use of the time-step manaager/advisor utility. It also exercizes
// the field time-step advisors.
// 
//===========================================================================//

class dummy_package {


// DATA
  private:

    rtt_timestep::ts_manager &tsm;
    dsxx::SP<rtt_timestep::field_ts_advisor> sp_te;
    dsxx::SP<rtt_timestep::field_ts_advisor> sp_ti;
    dsxx::SP<rtt_timestep::field_ts_advisor> sp_ri; 


// CREATORS

  public:
    dummy_package(rtt_timestep::ts_manager &tsm_);
    ~dummy_package();

// MANIPULATORS

// Method to advance the time-step

    void advance_state();

//ACCESSORS

    bool tests_passed() const;
};

#endif                          // __timestep_test_dummy_package_hh__

//---------------------------------------------------------------------------//
//                              end of dummy_package.hh
//---------------------------------------------------------------------------//
