//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Source_Builder.hh
 * \author Thomas M. Evans and Todd J. Urbatsch
 * \date   Thu Dec  9 10:31:16 1999
 * \brief  Header file for Rep_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Rep_Source_Builder_HH
#define RTT_imc_Rep_Source_Builder_HH

#include "Source_Builder.hh"
#include "ds++/Soft_Equivalence.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Rep_Source_Builder
 *
 * \brief Rep_Source_Builder is a derived rtt_imc::Source_Builder class for
 * full replication problems.
 *
 * Rep_Source_Builder is used to build rtt_imc::Source objects in full
 * replication topologies.
 *
 * As with the Source_Builder base class, Rep_Source_Builder is templated on
 * Mesh Type (MT) and Particle Type (PT).  The PT defaults to
 * rtt_imc::Particle.
 */
// revision history:
// -----------------
// 0) original
// 1) 19 Jun 2000 : added function recalc_census_ew_after_comb to better
//                  conserve energy after doing our reproducible comb.
// 2) 26 Jul 2000 : modified calc_num_part_and_rn_fields to better distribute 
//                  leftover particles (after integer math) across all procs.
// 3) 26 Jun 2001 : relaxed check of ecentot in constructor
// 4) 07 Jan 2002 : moved constructor to header file so that automatic
//                  instantiation will work; updated to work with new
//                  Particle_Stack 
// 5) 05-FEB-2002 : updated for multigroup
// 6) 10-FEB-2003 : added scoping for the Source_Builder base class
//===========================================================================//

template<class MT, class FT, class PT>
class Rep_Source_Builder : public Source_Builder<MT,FT,PT>
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

  public:
    // Typedef of base class for scoping.
    typedef Source_Builder<MT,FT,PT> Base;

  private:
    // Protected base class members that are used in this class.  We are
    // placing CCSF's here because gcc has trouble with operator() when the
    // base class scoping operator is applied.
    using Base::ecen;
    using Base::ew_cen;
    using Base::evol;
    using Base::evol_net;
    using Base::mat_vol_src;
    using Base::ew_vol;
    using Base::ess;
    using Base::ss_face_in_cell;
    using Base::ew_ss;
    using Base::volrn;
    using Base::ssrn;

  private:
    // Data fields unique or more properly defined for full Replication
    // topologies.

    // Number of global census particles.
    ccsf_int global_ncen;
    int global_ncentot;

    // Number of local census particles.
    ccsf_int local_ncen;
    int local_ncentot;

    // Number of global volume emission particles.
    ccsf_int global_nvol;
    int global_nvoltot;

    // Number of local volume emission particles.
    ccsf_int local_nvol;
    int local_nvoltot;

    // Number of global surface source particles.
    ccsf_int global_nss;
    int global_nsstot;

    // Number of local surface source particles.
    ccsf_int local_nss;
    int local_nsstot;

    // Global energy loss for volume emission, surface source, and census.
    double global_eloss_vol;
    double global_eloss_ss;
    double global_eloss_cen;

  private:
    // >>> IMPLEMENTATION

    // Calculate the local random number stream IDs and number of particles.
    void calc_num_part_and_rn_fields(const ccsf_int &, const int, int &,
				     ccsf_int &, int &, ccsf_int &);

    // Calculate the initial number of census particles.
    void calc_initial_ncen(ccsf_int &);

    // Recalculate the post-comb census energy-weights for full replication..
    void recalc_census_ew_after_comb(SP_Mesh, ccsf_int &, SP_Census, double &);

    // Calculate source numbers.
    void calc_source_numbers();

  public:
    // Constructor.
    template<class IT>
    Rep_Source_Builder(rtt_dsxx::SP<IT>, SP_Mesh, SP_Topology);

    // Build source.
    SP_Source build_Source(SP_Mesh, SP_Mat_State, SP_Opacity,
			   SP_Rnd_Control, SP_Comm_Patterns);

    // Calculate the initial census.
    void calc_initial_census(SP_Mesh, SP_Mat_State, SP_Opacity,
			     SP_Rnd_Control);

    // >>> IMPLEMENTATION OF BASE CLASS ACCESSORS
    
    //! Get global total initial census energy on processor.
    double get_initial_census_energy() const { return Base::ecentot; }

    //! Get global energy loss in volume emission on processor.
    double get_eloss_vol() const { return global_eloss_vol; }

    //! Get global, total number of volume emission particles on processor.
    int get_nvoltot() const { return global_nvoltot; }

    //! Get global number of volume emission particles / cell on processor.
    int get_nvol(int cell) const { return global_nvol(cell); }
    
    //! Get global energy loss in surface source on processor.
    double get_eloss_ss() const { return global_eloss_ss; } 

    //! Get global, total number of surface source particles on processor.
    int get_nsstot() const { return global_nsstot; }

    //! Get global energy loss in census on processor.
    double get_eloss_cen() const { return global_eloss_cen; }
    
    //! Get global, total number of post-comb census particles on processor.
    int get_ncentot() const { return global_ncentot; }
};

//---------------------------------------------------------------------------//
// TEMPLATE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for Rep_Source_Builder.
 */
template<class MT, class FT, class PT>
template<class IT>
Rep_Source_Builder<MT,FT,PT>::Rep_Source_Builder(rtt_dsxx::SP<IT> interface, 
						 SP_Mesh          mesh, 
						 SP_Topology      top)
    : Source_Builder<MT,FT,PT>(interface, mesh, top),
      global_ncen(mesh),
      global_ncentot(0),
      local_ncen(mesh),
      local_ncentot(0),
      global_eloss_cen(0),
      global_nvol(mesh),
      global_nvoltot(0),
      local_nvol(mesh),
      local_nvoltot(0),
      global_nss(mesh),
      global_nsstot(0),
      local_nss(mesh),
      local_nsstot(0),
      global_eloss_vol(0),
      global_eloss_ss(0)
{ 
    using rtt_dsxx::soft_equiv;

    Check(Base::parallel_data_op.check_global_equiv(rtt_rng::rn_stream));

    // if the census does not exist yet then we build it --> if the census
    // does not exist it means that this is the initial IMC cycle
    if (Base::census)
    {
	// sweep through cells and fill in persistent data

	// checks of persistent data
	double ecentot_check = 0.0;
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    // get global values of ecen from the interface
	    ecen(cell)     = interface->get_ecen(cell);
	    ecentot_check += ecen(cell);
	    
	    // more may follow --> probably cumulative edits
	}
	
	// get ecentot
	Base::ecentot = interface->get_ecentot();
    
	// checks
	Check(soft_equiv(Base::ecentot, ecentot_check,
			 mesh->num_cells() * 1.e-12));
    }
}

} // end namespace rtt_imc

#endif                          // RTT_imc_Rep_Source_Builder_HH

//---------------------------------------------------------------------------//
//                              end of imc/Rep_Source_Builder.hh
//---------------------------------------------------------------------------//
