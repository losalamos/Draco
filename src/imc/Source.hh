//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source.hh
 * \author Thomas M. Evans
 * \date   Thu May 14 08:45:49 1998
 * \brief  Source class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Source_hh__
#define __imc_Source_hh__

#include "Particle.hh"
#include "Particle_Buffer.hh"
#include "Mat_State.hh"
#include "Mesh_Operations.hh"
#include "mc/Topology.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <string>

namespace rtt_imc 
{

//===========================================================================//
/*!
 * \class Source
 
 * \brief Get source particles for IMC transport.

 * The Source class generates rtt_imc::Particle objects when queried with the
 * get_Source_Particle() function.  The Source has an overloaded bool()
 * operator that allows the user to query when the source is active or empty.

 */
// revision history:
// -----------------
// 0) original
// 1) 21-DEC-99: removed T^4 data and replaced it with the Mesh_Operations
//               class that does volume emission sampling using the tilt
// 2) 03-JUL-00: added Topology member data to map census particles from 
//               global cells to local cells
// 3) 24-AUG-00: added work-around capability for random number streams to
//               wrap around 2e9 seamlessly.
// 
//===========================================================================//

template<class MT, class PT = Particle<MT> >
class Source
{
  public:
    // Usefull typedefs
    typedef typename MT::template CCSF<int>      ccsf_int;
    typedef typename MT::template CCSF<double>   ccsf_double;
    typedef rtt_dsxx::SP<rtt_mc::Topology>       SP_Topology;
    typedef typename Particle_Buffer<PT>::Census PB_Census;
    typedef rtt_dsxx::SP<PB_Census>              SP_Census;
    typedef rtt_dsxx::SP<rtt_rng::Rnd_Control>   SP_Rnd_Control; 
    typedef rtt_dsxx::SP<Mat_State<MT> >         SP_Mat_State;
    typedef rtt_dsxx::SP<Mesh_Operations<MT> >   SP_Mesh_Op; 

  private:
    // Problem topology.
    SP_Topology topology;

    // volume source particles: number and first random number stream per
    // cell
    ccsf_int vol_rnnum;
    ccsf_int nvol;
    ccsf_double ew_vol;

    // surface source particles: number and first random number stream per
    // cell
    ccsf_int ss_rnnum;
    ccsf_int nss;
    ccsf_int fss;
    ccsf_double ew_ss;
    std::string ss_dist;

    // census bank
    SP_Census census;

    // total number of sources
    int nvoltot;
    int nsstot;
    int ncentot;

    // nsrcdone_cell is the running number of source particles completed for
    // a particular source type in a particular cell.
    int nsrcdone_cell;

    // cell currently under consideration
    int current_cell;

    // running totals of completed source particles, by type
    int nssdone;
    int nvoldone;
    int ncendone;

    // random number controller
    SP_Rnd_Control rcon;

    // Material State
    SP_Mat_State material;

    // Mesh_Operations class for performing volume emission sampling with T^4
    // tilts
    SP_Mesh_Op mesh_op;

  public:
    // Constructor.
    Source(ccsf_int &, ccsf_int &, ccsf_double &, ccsf_int &, ccsf_int &, 
	   ccsf_int &, ccsf_double &, SP_Census, std::string, int, int, 
	   SP_Rnd_Control, SP_Mat_State, SP_Mesh_Op, SP_Topology);

    // Required services for Source.
    rtt_dsxx::SP<PT> get_Source_Particle(double); 

    // Particle sources.
    rtt_dsxx::SP<PT> get_census(double);
    rtt_dsxx::SP<PT> get_evol(double);
    rtt_dsxx::SP<PT> get_ss(double);

    // Accessors.
    int get_nvoltot() const { return nvoltot; }
    int get_nsstot() const { return nsstot; }
    int get_ncentot() const { return ncentot; }

    // Source diagnostic.
    void print(std::ostream &) const;

    // Bool conversion.
    inline operator bool() const;

    // Get number of cells.
    int num_cells() const { return nvol.get_Mesh().num_cells(); }
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
inline std::ostream& operator<<(std::ostream &output, 
				const Source<MT,PT> &object) 
{
    object.print(output);
    return output;
}


} // end namespace rtt_imc

#endif                          // __imc_Source_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source.hh
//---------------------------------------------------------------------------//
