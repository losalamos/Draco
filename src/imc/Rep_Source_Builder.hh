//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Source_Builder.hh
 * \author Thomas M. Evans
 * \date   Thu Dec  9 10:31:16 1999
 * \brief  Header file for Rep_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Rep_Source_Builder_hh__
#define __imc_Rep_Source_Builder_hh__

#include "Source_Builder.hh"
#include "Source.hh"
#include "Particle_Buffer.hh"
#include "Opacity.hh"
#include "Mat_State.hh"
#include "mc/Topology.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>
#include <iostream>

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
// 
//===========================================================================//

template<class MT, class PT = Particle<MT> >
class Rep_Source_Builder : public Source_Builder<MT,PT>
{
  public:
    // typedefs used in the inheritance chain
    typedef dsxx::SP<Source<MT,PT> >                       SP_Source;
    typedef dsxx::SP<typename Particle_Buffer<PT>::Census> SP_Census;
    typedef dsxx::SP<Opacity<MT> >                         SP_Opacity;
    typedef dsxx::SP<Mat_State<MT> >                       SP_Mat_State;
    typedef dsxx::SP<MT>                                   SP_Mesh;
    typedef dsxx::SP<rtt_rng::Rnd_Control>                 SP_Rnd_Control;
    typedef dsxx::SP<rtt_mc::Topology>                     SP_Topology;
    typedef std::vector<int>                               sf_int;
    typedef std::vector<double>                            sf_double;
    typedef std::vector<std::string>                       sf_string;
    typedef std::vector<std::vector<double> >              vf_double;
    typedef std::string                                    std_string;
    typedef typename MT::CCSF_double                       ccsf_double;
    typedef typename MT::CCSF_int                          ccsf_int;
    typedef typename MT::CCVF_double                       ccvf_double;

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

    // Calculate the local random number stream IDs and number of particles.
    void calc_num_part_and_rn_fields(const ccsf_int &, const int, int &,
				     ccsf_int &, int &, ccsf_int &);

    // Calculate the initial number of census particles.
    void calc_initial_ncen(ccsf_int &);

    // Calculate source numbers.
    void calc_source_numbers();

  public:
    // Constructor.
    template<class IT>
    Rep_Source_Builder(dsxx::SP<IT>, SP_Mesh, SP_Topology);

    // Build source.
    SP_Source build_Source(SP_Mesh, SP_Mat_State, SP_Opacity,
			   SP_Rnd_Control);

    // Calculate the initial census.
    void calc_initial_census(SP_Mesh, SP_Mat_State, SP_Opacity,
			     SP_Rnd_Control);

    // IMPLEMENTATION OF BASE CLASS ACCESSORS

    //! Get global energy loss in volume emission on processor.
    double get_eloss_vol() const { return global_eloss_vol; }

    //! Get global, total number of volume emission particles on processor.
    int get_nvoltot() const { return global_nvoltot; }
    
    //! Get global energy loss in surface source on processor.
    double get_eloss_ss() const { return global_eloss_ss; } 

    //! Get global, total number of surface source particles on processor.
    int get_nsstot() const { return global_nsstot; }

    //! Get global energy loss in census on processor.
    double get_eloss_cen() const { return global_eloss_cen; }
    
    //! Get global, total number of post-comb census particles on processor.
    int get_ncentot() const { return global_ncentot; }
};

} // end namespace rtt_imc

#endif                          // __imc_Rep_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Rep_Source_Builder.hh
//---------------------------------------------------------------------------//
