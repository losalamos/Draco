//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/DD_Source_Builder.hh
 * \author Todd J. Urbatsch
 * \date   Tue May  2 14:40:43 2000
 * \brief  Header file for DD_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_DD_Source_Builder_hh__
#define __imc_DD_Source_Builder_hh__

#include "Source_Builder.hh"
#include "ds++/Soft_Equivalence.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class DD_Source_Builder
 *
 * \brief DD_Source_Builder is a derived rtt_imc::Source_Builder class for
 * full domain decomposition problems.
 *
 * DD_Source_Builder is used to build rtt_imc::Source object in a full 
 * domain decomposition topology.
 *
 * As with the Source_Builder base class, DD_Source_Builder is templated on
 * Mesh Type (MT) and Particle Type (PT).  The PT defaults to
 * rtt_imc::Particle.
 * 
 */
// revision history:
// -----------------
// 0) original
// 1) 19 Jun 2000 : added function recalc_census_ew_after_comb to better
//                  conserve energy after doing our reproducible comb.
// 2) 26 Jun 2001 : relaxed check of ecentot in constructor
// 3) 07 Jan 2002 : moved constructor to header file so that automatic
//                  instantiation will work; updated to work with new
//                  Particle_Stack 
// 5) 05-FEB-2002 : updated for multigroup
//===========================================================================//

template<class MT, class FT, class PT>
class DD_Source_Builder : public Source_Builder<MT,FT,PT>
{
  public:
    // typedefs used in the inheritance chain
    typedef rtt_dsxx::SP<Source<MT,FT,PT> >                  SP_Source;
    typedef typename rtt_mc::Particle_Containers<PT>::Census Census;
    typedef rtt_dsxx::SP<Census>                             SP_Census;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >                    SP_Opacity;
    typedef rtt_dsxx::SP<Mat_State<MT> >                     SP_Mat_State;
    typedef rtt_dsxx::SP<MT>                                 SP_Mesh;
    typedef rtt_dsxx::SP<rtt_rng::Rnd_Control>               SP_Rnd_Control; 
    typedef rtt_dsxx::SP<rtt_mc::Topology>                   SP_Topology;
    typedef rtt_dsxx::SP<rtt_mc::Comm_Patterns>              SP_Comm_Patterns;
    typedef std::vector<int>                                 sf_int;
    typedef std::vector<double>                              sf_double;
    typedef std::vector<std::string>                         sf_string;
    typedef std::vector<std::vector<double> >                vf_double;
    typedef std::string                                      std_string;
    typedef typename MT::template CCSF<double>               ccsf_double;
    typedef typename MT::template CCSF<int>                  ccsf_int;
    typedef typename MT::template CCVF<double>               ccvf_double;       

  private:
    // Data fields unique or more properly defined for full DD topologies.

    // Number of local census particles.
    ccsf_int local_ncen;
    int local_ncentot;

    // Total global number of census particles
    int global_ncentot; 

    // Number of local volume emission particles.
    ccsf_int local_nvol;
    int local_nvoltot;

    // Total global number of volume emission particles
    int global_nvoltot; 

    // Number of local surface source particles.
    ccsf_int local_nss;
    int local_nsstot;

    // Total global number of surface source particles
    int global_nsstot; 

    // Total global energies: census, volume emission, surface source
    double global_ecentot;
    double global_evoltot;
    double global_esstot;

    // local energy loss for volume emission, surface source, and census.
    double local_eloss_vol;
    double local_eloss_ss;
    double local_eloss_cen; 
 
    // global energy loss for volume emission, surface source, and census.
    double global_eloss_vol;
    double global_eloss_ss;
    double global_eloss_cen; 
  
    // Calculate the local random number stream IDs
    void calc_fullDD_rn_fields(const int, ccsf_int &, int &, ccsf_int &, int);

    // Calculate the initial number of census particles for full DD topology.
    void calc_initial_ncen(ccsf_int &);

    // Recalculate the post-comb census energy-weights for full DD topology.
    void recalc_census_ew_after_comb(SP_Mesh, ccsf_int &, SP_Census, double &);

    // Calculate source numbers for full DD topology.
    void calc_source_numbers();

  public:
    // Constructor
    template<class IT>
    DD_Source_Builder(rtt_dsxx::SP<IT>, SP_Mesh, SP_Topology);

    // Build source.
    SP_Source build_Source(SP_Mesh, SP_Mat_State, SP_Opacity,
			   SP_Rnd_Control, SP_Comm_Patterns); 
	

    // Calculate the initial census.
    void calc_initial_census(SP_Mesh, SP_Mat_State, SP_Opacity,
			     SP_Rnd_Control);

    // IMPLEMENTATION OF BASE CLASS ACCESSORS

    //! Get global total intial census energy.
    inline double get_initial_census_energy() const;

    //! Get global energy loss in volume emission.
    double get_eloss_vol() const { return global_eloss_vol; }

    //! Get global total number of volume emission particles.
    int get_nvoltot() const { return global_nvoltot; }

    //! Get local number of volume emission particles / cell on processor.
    int get_nvol(int cell) const { return local_nvol(cell); }
    
    //! Get global energy loss in surface source.
    double get_eloss_ss() const { return global_eloss_ss; } 

    //! Get global total number of surface source particles.
    int get_nsstot() const { return global_nsstot; }

    //! Get global, total energy loss in census.
    double get_eloss_cen() const { return global_eloss_cen; }
    
    //! Get global, total number of post-comb census particles.
    int get_ncentot() const { return global_ncentot; }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

template<class MT, class FT, class PT>
double DD_Source_Builder<MT,FT,PT>::get_initial_census_energy() const
{
    double energy = ecentot;
    C4::gsum(energy);
    return energy;
}

//---------------------------------------------------------------------------//
// TEMPLATE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for DD_Source_Builder.
 */
template<class MT, class FT, class PT>
template<class IT>
DD_Source_Builder<MT,FT,PT>::DD_Source_Builder(rtt_dsxx::SP<IT> interface, 
					       SP_Mesh          mesh, 
					       SP_Topology      top)
    : Source_Builder<MT,FT,PT>(interface, mesh, top),
    local_ncen(mesh),
    local_ncentot(0),
    local_eloss_cen(0),
    global_eloss_cen(0),
    local_nvol(mesh),
    local_nvoltot(0),
    local_nss(mesh),
    local_nsstot(0),
    local_eloss_vol(0),
    global_eloss_vol(0),
    local_eloss_ss(0),
    global_eloss_ss(0)
{ 
    // at the beginning of the timestep, random number stream ID should be
    // the same on every processor.
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
    
    // Update the persistent census energy data, unless the census does not
    // exist yet, which is the case on the first IMC cycle.
    if (census)
    {
	// a check on the total census energy
	double ecentot_check = 0.0;

	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    // get local values of ecen from the interface
	    ecen(cell)     = interface->get_ecen(cell);
	    ecentot_check += ecen(cell);

	    // more updates may follow -- probably time-cumulative edits
	}

	// sum up census energy check from all processors
	C4::gsum(ecentot_check);

	// get total global census energy from the interface
	global_ecentot = interface->get_ecentot();

	// check consistency of energies and totals
	Check (rtt_mc::global::soft_equiv(global_ecentot, ecentot_check,
					  mesh->num_cells() * 1.0e-12));
    }
}

} // end namespace rtt_imc

#endif                          // __imc_DD_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/DD_Source_Builder.hh
//---------------------------------------------------------------------------//
