//----------------------------------*-C++-*----------------------------------//
// Source.hh
// Thomas M. Evans
// Thu May 14 08:45:49 1998
//---------------------------------------------------------------------------//
// @> Source class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Source_hh__
#define __imctest_Source_hh__

//===========================================================================//
// class Source - 
//
// Purpose : produce a particle for IMC transport
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Particle.hh"
#include "imctest/Particle_Buffer.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <fstream>
#include <string>

IMCSPACE

using std::ifstream;
using std::string;

template<class MT>
class Source
{
private:
  // volume source particle counters
    typename MT::CCSF_int vol_rnnum;
    typename MT::CCSF_int nvol;

  // surface source particle counters
    typename MT::CCSF_int ss_rnnum;
    typename MT::CCSF_int nss;

  // census file
    ifstream census;

  // total number of sources
    int nvoltot;
    int nsstot;
    int ncentot;

  // Particle Buffer
    Particle_Buffer<Particle<MT> > buffer;

public:
  // constructor
    Source(typename MT::CCSF_int &, typename MT::CCSF_int &,
	   typename MT::CCSF_int &, typename MT::CCSF_int &, string, int, 
	   int, int);

  // required services for Source
    inline SP<Particle> get_Particle(); 

  // Particle sources
    SP<Particle> get_census();
    SP<Particle> get_evol();
    SP<Particle> get_ss();
};

CSPACE

#endif                          // __imctest_Source_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Source.hh
//---------------------------------------------------------------------------//
