//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Random_Walk_Sub_Tally.hh
 * \author Thomas M. Evans
 * \date   Tue Jun 24 13:43:28 2003
 * \brief  Random_Walk_Sub_Tally class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Random_Walk_Sub_Tally_hh
#define rtt_imc_Random_Walk_Sub_Tally_hh

#include "ds++/Assert.hh"

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Random_Walk_Sub_Tally
 * \brief Tally random walk related statistics and quantities.
 *
 * The Random_Walk_Sub_Tally class is a \ref sub_tally compliant class.  It
 * is designed to be used in the rtt_imc::Tally class.  It provides
 * accumulators and accessors for edit quantities specific to random walk
 * phenomenon. 
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Random_Walk_Sub_Tally 
{
  private:
    // >>> DATA

    // Number of random walks.
    int n_random_walks;

    // Radius of random walk spheres.
    double sphere_radius;
    int    n_spheres_generated;

    // Length of random walk steps.
    double rw_step_length;

  public:
    // Constructor.
    Random_Walk_Sub_Tally();

    // >>> ACCUMULATORS
    
    // Accumulate number of random walks.
    inline void accum_n_random_walks(const int n = 1);

    // Accumulate random walk sphere radii.
    inline void accum_sphere_radii(const double radius);

    // Accumulate random walk step length.
    inline void accum_step_length(const double length);

    // >>> ACCESSORS

    //! Get number of random walks.
    int get_accum_n_random_walks() const { return n_random_walks; }

    //! Get accumulated random walk sphere radii.
    double get_accum_sphere_radii() const { return sphere_radius; }

    //! Get number of random walk spheres generated.
    int get_accum_n_spheres() const { return n_spheres_generated; }

    //! Get accumulated random walk step lengths.
    double get_accum_step_lengths() const { return rw_step_length; }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate the number of random walks made by particles.
 */
void Random_Walk_Sub_Tally::accum_n_random_walks(const int n)
{ 
    Check (n >= 0);
    n_random_walks += n;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate the number and radii of generated random walk spheres.
 */
void Random_Walk_Sub_Tally::accum_sphere_radii(const double radius)
{
    Check(radius >= 0.0);
    sphere_radius += radius;
    n_spheres_generated++;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Accumulate random walk step length.
 */
void Random_Walk_Sub_Tally::accum_step_length(const double length)
{
    Check (length >= 0.0);
    rw_step_length += length;
}

} // end namespace rtt_imc

#endif // rtt_imc_Random_Walk_Sub_Tally_hh

//---------------------------------------------------------------------------//
//              end of imc/Random_Walk_Sub_Tally.hh
//---------------------------------------------------------------------------//
