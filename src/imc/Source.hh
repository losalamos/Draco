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
#include "imc/Mat_State.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <string>

IMCSPACE

// draco components
using RNG::Rnd_Control;
using dsxx::SP;

// STL components
using std::string;
using std::ostream;

template<class MT, class PT=Particle<MT> >
class Source
{
public:
  // volume source particles: number and first random number stream per cell
    typename MT::CCSF_int vol_rnnum;
    typename MT::CCSF_int nvol;
    typename MT::CCSF_double ew_vol;
    typename MT::CCVF_double t4_slope;

  // surface source particles: number and first random number stream per cell
    typename MT::CCSF_int ss_rnnum;
    typename MT::CCSF_int nss;
    typename MT::CCSF_int fss;
    typename MT::CCSF_double ew_ss;
    string ss_dist;

  // census bank
    typename Particle_Buffer<PT>::Census census;

  // total number of sources
    int nvoltot;
    int nsstot;
    int ncentot;

  // nsrcdone_cell is the running number of source particles completed for a
  // particular source type in a particular cell. 
    int nsrcdone_cell;

  // cell currently under consideration
    int current_cell;

  // running totals of completed source particles, by type
    int nssdone;
    int nvoldone;
    int ncendone;

  // random number controller
    SP<Rnd_Control> rcon;

  // Particle Buffer
    Particle_Buffer<PT> buffer;

  // Material State
    SP<Mat_State<MT> > material;

public:
  // constructor
    Source(typename MT::CCSF_int &, typename MT::CCSF_int &,
	   typename MT::CCSF_double &, typename MT::CCVF_double &,
	   typename MT::CCSF_int &, typename MT::CCSF_int &, 
	   typename MT::CCSF_int &, typename MT::CCSF_double &,
	   typename Particle_Buffer<PT>::Census &,
	   string, int, int, SP<Rnd_Control>, 
	   const Particle_Buffer<PT> &, SP<Mat_State<MT> >);

  // required services for Source
    SP<PT> get_Source_Particle(double); 

  // Particle sources
    SP<PT> get_census(double);
    SP<PT> get_evol(double);
    SP<PT> get_ss(double);

  // accessors
    int get_nvoltot() const { return nvoltot; }
    int get_nsstot() const { return nsstot; }
    int get_ncentot() const { return ncentot; }

  // source diagnostic
    void print(ostream &) const;

  // bool conversion
    inline operator bool() const;
};

//---------------------------------------------------------------------------//
// inline functions for Source
//---------------------------------------------------------------------------//
// bool conversion operator

template<class MT, class PT>
inline Source<MT, PT>::operator bool() const
{
    return (ncentot != ncendone || nsstot != nssdone || nvoltot != nvoldone);
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// output ascii version of the source

template<class MT, class PT>
inline ostream& operator<<(ostream &output, const Source<MT,PT> &object)
{
    object.print(output);
    return output;
}


CSPACE

#endif                          // __imc_Source_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source.hh
//---------------------------------------------------------------------------//
