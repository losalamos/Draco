//----------------------------------*-C++-*----------------------------------//
// IMCTEST_Man.hh
// Thomas M. Evans
// Wed Jun  3 10:36:11 1998
//---------------------------------------------------------------------------//
// @> IMCTEST_Man class header file.
//---------------------------------------------------------------------------//

#ifndef __imc_IMCTEST_Man_hh__
#define __imc_IMCTEST_Man_hh__

//===========================================================================//
// class IMCTEST_Man - 
//
// Purpose : Manager for running IMCTEST package.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"

template<class MT> class Particle;

IMCSPACE

template<class MT, class PT = Particle<MT> >
class IMCTEST_Man 
{
private:

public:
  // constructor
    IMCTEST_Man();
    IMCTEST_Man( const IMCTEST_Man& );
    ~IMCTEST_Man();
    IMCTEST_Man& operator=( const IMCTEST_Man& );
};

CSPACE

#endif                          // __imc_IMCTEST_Man_hh__

//---------------------------------------------------------------------------//
//                              end of imc/IMCTEST_Man.hh
//---------------------------------------------------------------------------//
