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

template<class MT, class PT=Particle<MT> >
class Source
{
public:
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
    Particle_Buffer<PT> buffer;

public:
  // constructor
    Source(typename MT::CCSF_int &, typename MT::CCSF_int &,
	   typename MT::CCSF_int &, typename MT::CCSF_int &, string, 
	   int, int, int, SP<Rnd_Control>, const Particle_Buffer<PT> &);

  // required services for Source
    inline SP<PT> get_Particle(); 

  // Particle sources
    SP<PT> get_census();
    SP<PT> get_evol();
    SP<PT> get_ss();
};

CSPACE

#endif                          // __imc_Source_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source.hh
//---------------------------------------------------------------------------//
