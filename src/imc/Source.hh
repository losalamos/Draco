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

 * The Source class generates the next rtt_imc::Particle object when queried
 * with the get_Source_Particle() function.  The Source has an overloaded
 * bool() operator that allows the user to query when the source is active or
 * empty.

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
// 4) 14-DEC-00: added accessor to get the total number of source particles 
//               contained in this source
// 5) 31-JUL-01: used INTEGER_MODULO_1E9 from rtt_mc for random num. stream
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
    // >>> DATA

    // Problem topology.
    SP_Topology topology;

    // Volume source particles: number and first random number stream per
    // cell.
    ccsf_int vol_rnnum;
    ccsf_int nvol;
    ccsf_double ew_vol;

    // Surface source particles: number and first random number stream per
    // cell.
    ccsf_int ss_rnnum;
    ccsf_int nss;
    ccsf_int fss;
    ccsf_double ew_ss;
    std::string ss_dist;

    // Census bank.
    SP_Census census;

    // Total number of sources.
    int nvoltot;
    int nsstot;
    int ncentot;

    // The running number of source particles completed for a particular
    // source type in a particular cell.
    int nsrcdone_cell;

    // Cell currently under consideration.
    int current_cell;

    // Running totals of completed source particles, by type.
    int nssdone;
    int nvoldone;
    int ncendone;

    // Random number controller.
    SP_Rnd_Control rcon;

    // Material State.
    SP_Mat_State material;

    // Mesh_Operations class for performing volume emission sampling with T^4
    // tilts.
    SP_Mesh_Op mesh_op;

  private:
    // >>> PRIVATE IMPLEMENTATION

    // Particle sources accessors.
    rtt_dsxx::SP<PT> get_census(double);
    rtt_dsxx::SP<PT> get_evol(double);
    rtt_dsxx::SP<PT> get_ss(double);

  public:
    // Constructor.
    Source(ccsf_int &, ccsf_int &, ccsf_double &, ccsf_int &, ccsf_int &, 
	   ccsf_int &, ccsf_double &, SP_Census, std::string, int, int, 
	   SP_Rnd_Control, SP_Mat_State, SP_Mesh_Op, SP_Topology);

    // Primary source service, get a source particle.
    rtt_dsxx::SP<PT> get_Source_Particle(double); 

    //! Get the number of volume emission sources contained in this source. 
    int get_nvoltot() const { return nvoltot; }

    //! Get the number of surface sources contained in this source object.
    int get_nsstot() const { return nsstot; }

    //! Get the number of census particles contained in this source object.
    int get_ncentot() const { return ncentot; }

    // Get the total number of source particles contained in this object.
    inline int get_num_source_particles() const;

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
/*!

 * \brief Boolean conversion operator for source object.

 * \return true if source particles exist; false if the source is empty

 */
template<class MT, class PT>
Source<MT, PT>::operator bool() const
{
    return (ncentot != ncendone || nsstot != nssdone || nvoltot != nvoldone);
}

//---------------------------------------------------------------------------//
/*!
 
 * \brief Get the total number of source particles in this source object.

 * This function returns the total number of source particles initially
 * contained in this source object. It is the sum of all three source
 * species: volume, surface source, and census.

 * \return the sum of the volume, surface source, and census particles
 * contained in this source object.
 
 */
template<class MT, class PT>
int Source<MT,PT>::get_num_source_particles() const
{
    return nvoltot + nsstot + ncentot;
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
