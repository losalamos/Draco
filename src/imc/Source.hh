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

#include "Mesh_Operations.hh"
#include "Global.hh"
#include "Opacity.hh"
#include "Mat_State.hh"
#include "mc/Particle_Stack.hh"
#include "mc/Topology.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <string>
#include <vector>

namespace rtt_imc 
{

// Forward declarations.
class Gray_Frequency;
class Multigroup_Frequency;
template<class MT> class Gray_Particle;
template<class MT> class Multigroup_Particle;

//===========================================================================//
/*!
 * \class Frequency_Sampling_Data
 *
 * \brief Frequency-dependent sampling struct.
 *
 * This class holds data that is used by source and source builder
 * specializations on frequency type (Gray_Frequency, Multigroup_Frequency).
 */
//===========================================================================//

template<class MT, class FT>
struct Frequency_Sampling_Data
{
    typedef rtt_dsxx::SP<MT> SP_Mesh;

    Frequency_Sampling_Data(SP_Mesh mesh) {/*...*/}
};

//---------------------------------------------------------------------------//
/*!
 * \brief Specialization of Frequency_Sampling_Data on Multigroup_Frequency.
 *
 * This class holds data needed by multigroup specializations of the
 * rtt_imc::Source and rtt_imc::Source_Builder classes.
 */
//---------------------------------------------------------------------------//

template<class MT>
struct Frequency_Sampling_Data<MT, Multigroup_Frequency>
{
    typedef rtt_dsxx::SP<MT>                   SP_Mesh;
    typedef typename MT::template CCSF<double> ccsf_double;

    // Probability that a volume emission source particle is emitted from an
    // opacity weighted Planckian as opposed to a straight Planckian.
    ccsf_double prob_of_straight_Planck_emission;

    // Surface source temperature, in keV, for each cell.  Needed for
    // sampling frequency group of a surface source particle in a multigroup
    // treatment.
    ccsf_double ss_temperature;

    Frequency_Sampling_Data(SP_Mesh mesh)
	: prob_of_straight_Planck_emission(mesh),
	  ss_temperature(mesh) {/*...*/}
};

//===========================================================================//
/*!
 * \class Source
 *
 * \brief Get source particles for IMC transport.
 *
 * The Source class generates the next rtt_imc::Particle object when queried
 * with the get_Source_Particle() function.  The Source has an overloaded
 * bool() operator that allows the user to query when the source is active or
 * empty.
 *
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
// 6) 05-FEB-02: updated for multigroup
// 7) 18-MAR-02: added ss_temperature to Frequency_Sampling_Data.
// 
//===========================================================================//

template<class MT, class FT, class PT>
class Source
{
  public:
    // Usefull typedefs
    typedef typename MT::template CCSF<int>                    ccsf_int;
    typedef typename MT::template CCSF<double>                 ccsf_double;
    typedef rtt_dsxx::SP<rtt_mc::Topology>                     SP_Topology;
    typedef typename rtt_mc::Particle_Containers<PT>::Census   PB_Census;
    typedef rtt_dsxx::SP<PB_Census>                            SP_Census;
    typedef rtt_dsxx::SP<rtt_rng::Rnd_Control>                 SP_Rnd_Control; 
    typedef rtt_dsxx::SP<Mat_State<MT> >                       SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >                      SP_Opacity;
    typedef rtt_dsxx::SP<Mesh_Operations<MT> >                 SP_Mesh_Op;
    typedef Gray_Particle<MT>                                  Gray_PT;
    typedef Multigroup_Particle<MT>                            MG_PT;
    typedef rtt_dsxx::SP<Gray_PT>                              SP_Gray_PT;
    typedef rtt_dsxx::SP<MG_PT>                                SP_MG_PT;
    typedef rtt_imc::global::Type_Switch<Gray_Frequency>       Switch_Gray;
    typedef rtt_imc::global::Type_Switch<Multigroup_Frequency> Switch_MG;
    typedef rtt_imc::global::Type_Switch<Gray_PT>              Switch_Gray_PT;
    typedef rtt_imc::global::Type_Switch<MG_PT>                Switch_MG_PT;
    typedef std::vector<double>                                sf_double;
    typedef typename rtt_imc::global::Type_Switch<FT>::Type    Dummy_Type;

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

    // Opacity
    SP_Opacity opacity;

    // Mesh_Operations class for performing volume emission sampling with T^4
    // tilts.
    SP_Mesh_Op mesh_op;

    // Probabilities for sampling a volume emission source that is emitted
    // from a straight Planckian.
    Frequency_Sampling_Data<MT,FT> freq_samp_data;

  private:
    // >>> PRIVATE IMPLEMENTATION

    // Particle sources accessors.
    rtt_dsxx::SP<PT> get_census(double);
    rtt_dsxx::SP<PT> get_evol(double);
    rtt_dsxx::SP<PT> get_ss(double);

    // Make a Gray_Particle for surface source.
    template<class Stop_Explicit_Instantiation>
    SP_Gray_PT make_ss_particle(Switch_Gray, Switch_Gray_PT, 
				const sf_double &, const sf_double &,
				const double, int, const rtt_rng::Sprng &,
				const double, const double);

    // Make a Multigroup_Particle for surface source.
    template<class Stop_Explicit_Instantiation>
    SP_MG_PT make_ss_particle(Switch_MG, Switch_MG_PT, 
			      const sf_double &, const sf_double &,
			      const double, int, const rtt_rng::Sprng &,
			      const double, const double);

    // Make a Gray_Particle for surface source.
    template<class Stop_Explicit_Instantiation>
    SP_Gray_PT make_vol_particle(Switch_Gray, Switch_Gray_PT, 
				 const sf_double &, const sf_double &,
				 const double, int, const rtt_rng::Sprng &,
				 const double, const double);

    // Make a Multigroup_Particle for surface source.
    template<class Stop_Explicit_Instantiation>
    SP_MG_PT make_vol_particle(Switch_MG, Switch_MG_PT, 
			       const sf_double &, const sf_double &,
			       const double, int, const rtt_rng::Sprng &,
			       const double, const double);

  public:
    // Constructor.
    Source(ccsf_int &, ccsf_int &, ccsf_double &, ccsf_int &, ccsf_int &, 
	   ccsf_int &, ccsf_double &, SP_Census, std::string, int, int, 
	   SP_Rnd_Control, SP_Mat_State, SP_Mesh_Op, SP_Topology, 
	   SP_Opacity, const Frequency_Sampling_Data<MT,FT> &);

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
template<class MT, class FT, class PT>
Source<MT,FT,PT>::operator bool() const
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
template<class MT, class FT, class PT>
int Source<MT,FT,PT>::get_num_source_particles() const
{
    return nvoltot + nsstot + ncentot;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// output ascii version of the source

template<class MT, class FT, class PT>
inline std::ostream& operator<<(std::ostream &output, 
				const Source<MT,FT,PT> &object) 
{
    object.print(output);
    return output;
}


} // end namespace rtt_imc

#endif                          // __imc_Source_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source.hh
//---------------------------------------------------------------------------//
