//----------------------------------*-C++-*----------------------------------//
// Source.hh
// Thomas M. Evans
// Thu May 14 08:45:49 1998
//---------------------------------------------------------------------------//
// @> Source class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Source_hh__
#define __imc_Source_hh__

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

#include "imc/Names.hh"
#include "imc/Particle.hh"
#include "imc/Particle_Buffer.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <fstream>
#include <string>

IMCSPACE

// draco components
using RNG::Rnd_Control;

// STL components
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

  // random number controller
    SP<Rnd_Control> rcon;

  // Particle Buffer
    Particle_Buffer<Particle<MT> > buffer;

public:
  // constructor
    Source(typename MT::CCSF_int &, typename MT::CCSF_int &,
	   typename MT::CCSF_int &, typename MT::CCSF_int &, string, int, 
	   int, int, SP<Rnd_Control>);

  // required services for Source
    inline SP<Particle<MT> > get_Particle(); 

  // Particle sources
    SP<Particle<MT> > get_census();
    SP<Particle<MT> > get_evol();
    SP<Particle<MT> > get_ss();
};

CSPACE

#endif                          // __imc_Source_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source.hh
//---------------------------------------------------------------------------//
