//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Random_Walk.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 14 18:26:24 2003
 * \brief  Random_Walk class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Random_Walk_hh
#define rtt_imc_Random_Walk_hh

#include "Hybrid_Diffusion.hh"
#include "Diffusion_Opacity.hh"
#include "ds++/SP.hh"
#include <vector>

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
    Random_Walk(SP_Mesh, SP_Diffusion_Opacity, double = 1.0);

    // Number of cells in the mesh.
    int num_cells() const { return mesh->num_cells(); }

    // See if we should do Random Walk on this step.
    bool do_a_random_walk(int, double, double, double);

    // Do a random walk.
    double random_walk(sf_double &, sf_double &, double &, int, 
		       Rnd_Type, bool &);

    // >>> QUERY FUNCTIONS
};

} // end namespace rtt_imc

#endif                          // rtt_imc_Random_Walk_hh

//---------------------------------------------------------------------------//
//                              end of imc/Random_Walk.hh
//---------------------------------------------------------------------------//
