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
using std::ostream;
using std::string;

template<class MT, class PT=Particle<MT> >
class Source
{
public:
  // volume source particles: number and first random number stream per cell
    typename MT::CCSF_int vol_rnnum;
    typename MT::CCSF_int nvol;

  // surface source particles: number and first random number stream per cell
    typename MT::CCSF_int ss_rnnum;
    typename MT::CCSF_int nss;

  // census file
    ifstream census;

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

public:
  // constructor
    Source(typename MT::CCSF_int &, typename MT::CCSF_int &,
	   typename MT::CCSF_int &, typename MT::CCSF_int &, 
	   string, int, int, int, SP<Rnd_Control>, 
	   const Particle_Buffer<PT> &);

  // required services for Source
    inline SP<PT> get_Source_Particle(); 

  // Particle sources
    SP<PT> get_census();
    SP<PT> get_evol();
    SP<PT> get_ss();

  // source diagnostic
    void print(ostream &) const;
};

//---------------------------------------------------------------------------//
// inline functions for Source
//---------------------------------------------------------------------------//
// get a source particle

// template<class MT, class PT>
// inline SP<PT> Source<MT, PT>::get_Source_Particle()
// {
//     bool sampled = false;

//   // instantiate particle to return
//     SP<Particle<MT> > source_particle;

//   // do all surface source particles, one-by-one
//     while (!sampled && nssdone < nsstot)
//     {
// 	if (nsrcdone_cell < nss(current_cell))
// 	{
// 	    source_particle = get_ss();
// 	    sampled = true;
// 	    nsrcdone_cell++;
// 	    nssdone++;
// 	}
// 	else
// 	{
// 	    current_cell++;
// 	    nsrcdone_cell = 0;
// 	    if (current_cell > numcells)
// 	    {
// 		Check (nssdone == nsstot);
// 		current_cell = 1;
// 	    }
// 	}
//     }

//   // do all volume emission particles, one-by-one
//     while (!sampled && nvoldone < nvoltot)
//     {
// 	if (nsrcdone_cell < nvol(current_cell))
// 	{
// 	    source_particle = get_evol();
// 	    sampled = true;
// 	    nsrcdone_cell++;
// 	    nvoldone++;
// 	}
// 	else
// 	{
// 	    current_cell++;
// 	    nsrcdone_cell = 0;
// 	    if (current_cell > numcells)
// 	    {
// 		Check (nvoldone == nvoltot);
// 		current_cell = 1;
// 	    }
// 	}
//     }

//   // do all census particles, one-by-one
//     while (!sampled && ncendone < ncentot)
//     {
// 	    source_particle = get_census();
// 	    sampled = true;
// 	    ncendone++;
//     }

//     Check (sampled);

//     return source_particle;
// }


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
