//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Random_Walk.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 14 18:26:24 2003
 * \brief  Random_Walk class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Random_Walk_hh
#define rtt_imc_Random_Walk_hh

#include <algorithm>
#include <vector>
#include <utility>
#include "ds++/SP.hh"
#include "mc/Constants.hh"
#include "Opacity.hh"
#include "Hybrid_Diffusion.hh"
#include "Diffusion_Opacity.hh"

// Forward declarations.
namespace rtt_rng
{
class Sprng;
}

namespace rtt_imc
{
    
//===========================================================================//
/*!
 * \class Random_Walk_Sampling_Tables
 * 
 * \brief Holds data for sampling Random Walk.
 *
 * \sa imc/test/tstRandom_Walk.cc for examples.
 */
// revision history:
// -----------------
// 0) original
// 1) 3-SEP-03 : updated Random Walk conditions to by 5*eff_mfp, see 
//               do_a_random_walk()
// 
//===========================================================================//

class Random_Walk_Sampling_Tables
{
  private: 
    // >>> DATA

    // Abcissas for sampling tables (Dt/Ro^2).
    double a[41];
    double a_R[34];

    // Probability table for particle exiting sphere.
    double prob_exit[41];

    // Probability table for particle that goes to census without escaping
    // sphere as a function of [a][b] where a=Dtcen/Ro^2 and b=R1/Ro.
    double prob_R_cen[34][13];

    // Ratios of R1/R0 for prob_R_cen table.
    double b[13];

  private:
    // >>> IMPLEMENTATION

    // Set abcissas.
    void set_abcissas();

    // Set probability of exiting sphere.
    void set_prob_exit();

    // Set probability of staying in sphere as a function of R1.
    void set_prob_R_cen();

    // Disallow copy constructor.
    Random_Walk_Sampling_Tables(const Random_Walk_Sampling_Tables &);

    // Disallow assignment.
    const Random_Walk_Sampling_Tables& operator=(
	const Random_Walk_Sampling_Tables &);

  public:
    // Constructor.
    Random_Walk_Sampling_Tables();

    // Get probability of exiting sphere for a particle.
    double get_prob_exit(double t, double D, double Ro) const;
	
    // Get elapsed time particle spends inside a sphere.
    double get_elapsed_time(double D, double Ro, double ran) const;

    // Get radius of sphere for particle that goes to census.
    double get_radius(double t, double D, double Ro, double ran) const;
};
 
//===========================================================================//
/*!
 * \class Random_Walk
 * 
 * \brief Does hybrid diffusion/IMC with Fleck and Canfield's Random Walk
 * algorithm.
 *
 */
/*!
 * \example imc/test/tstRandom_Walk.cc
 * Tests Random_Walk class and Random Walk sampling tables.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class Random_Walk : public Hybrid_Diffusion
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Diffusion_Opacity<MT> > SP_Diffusion_Opacity;
    typedef std::vector<double>                  sf_double;
    typedef rtt_rng::Sprng                       Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>               SP_Rnd_Type;
    typedef rtt_dsxx::SP<MT>                     SP_Mesh;
    typedef std::pair<double, double>            pair_double;

  private:
    // >>> DATA

    // Data for sampling Random Walk.
    const Random_Walk_Sampling_Tables table;

    // Smart pointer to the mesh.
    SP_Mesh mesh;

    // Diffusion Opacity object.
    SP_Diffusion_Opacity diff_opacity;

    // Radius of random walk sphere.
    double rw_radius;
    bool   rw_set;

    // Minimum optical radius for doing a random walk step.
    double min_optical_radius;
    
  public:
    // Constructor.
    Random_Walk(SP_Mesh, SP_Diffusion_Opacity, double = 5.0);

    // Number of cells in the mesh.
    int num_cells() const { return mesh->num_cells(); }

    // See if we should do Random Walk on this step.
    template<class FT>
    bool do_a_random_walk(int, double, double, double, const Opacity<MT,FT> &);

    // Do a random walk.
    pair_double random_walk(sf_double &, sf_double &, double &, int, 
			    Rnd_Type, bool &);

    // Reset the random walk.
    void reset();
};

//---------------------------------------------------------------------------//
// TEMPLATE MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief See if we should do Random Walk on this step.
 *
 * Random Walk is valid during a step if:
 * - Ro > min_optical_radius * effective-mean-free-path AND
 * - Ro > distance-to-collision AND
 * - distance-to-census > distance-to-collision
 * .
 * Here, the effective-mean-free-path is the maximum of the effective
 * Rosseland and Planck mean free paths.
 *
 * \param cell current cell index
 * \param rw_radius Random Walk sphere radius
 * \param d_collision distance to collision (cm)
 * \param d_census distance to census (shakes)
 * \param opacity rtt_imc::Opacity object
 */
template<class MT>
template<class FT>
bool Random_Walk<MT>::do_a_random_walk(int                   cell,
				       double                rw_radius_in,
				       double                d_collision,
				       double                d_census,
				       const Opacity<MT,FT> &opacity)
{
    Require (cell > 0);
    Require (cell <= diff_opacity->num_cells());
    Require (rw_radius_in >= 0.0);
    Require (rw_set == false);

    // get the gray Planck and Rosseland opacities
    double planck = opacity.get_Planck(cell);
    double ross   = diff_opacity->get_Rosseland_opacity(cell);

    // calculate the effective mean free paths: NOTE, the gray scattering
    // opacity may be inconsistent with the Planck or Rosseland gray
    // absorption opacities; in multigroup, the gray scattering is done by
    // Rosseland integration; in gray, the gray scattering is provided
    // through the interface; we make no attempt to store consistent gray
    // scattering cross sections (ie, Rosseland and Planck scattering)

    // get the effective Rosseland mean free path
    double rmfp = diff_opacity->get_Rosseland_effmfp(cell);

    // calculate the effective Planck mean free paths
    double pmfp  = rtt_mc::global::huge;
    double denom = (1.0 - diff_opacity->get_fleck(cell)) * planck + 
	diff_opacity->get_gray_scattering(cell);
    if (denom > 0.0)
	pmfp = 1.0 / denom;

    // choose the mfp to use for random walk
    double mfp = min_optical_radius * std::max(rmfp, pmfp);

    // set the random walk radius for this step
    rw_radius = rw_radius_in;

    // check to see if random walk is valid

    // rw_radius must be greater than the minimum number of mfps
    if (rw_radius <= mfp)
	return false;

    // rw_radius must be greater than the distance to collision
    else if (rw_radius <= d_collision)
	return false;

    // rw_radius is greater than dist-to-collision so rw is possible;
    // however, if d_census is the limiting event then don't do random walk 
    else if (d_census <= d_collision)
	return false;

    // set the rw_set flag to true (random walk is on)
    rw_set = true;

    // return random walk boolean flag
    return true;
}

} // end namespace rtt_imc

#endif                          // rtt_imc_Random_Walk_hh

//---------------------------------------------------------------------------//
//                              end of imc/Random_Walk.hh
//---------------------------------------------------------------------------//
