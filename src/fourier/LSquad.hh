//----------------------------------*-C++-*----------------------------------//
// LSquad.hh
// William D. Hawkins
// Mon Jun 14 10:07:47 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __quad_LSquad_hh__
#define __quad_LSquad_hh__

// Level-Symmetric 3-D Angular Quadrature Class for Isotropic Scattering
// Returns quadrature points and weights for S2 - S24
// LSquad.hh: Declaration of the LSquad class
// LSquad.cc: Member function definitions

#include <vector>

using std::vector;

namespace rtt_quad
{

const double I_PI = 3.14159265358979323846;

//===========================================================================//
// class LSquad - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class LSquad
{
  public:
    // Default Constructor
    LSquad(int = 2, double = 4.0*I_PI);
    LSquad(int levels_, int dirs_, int odirs_, double norm_,
	   const vector<double> &mu_, const vector<double> &eta_,
	   const vector<double> &xi_, const vector<double> &wt_);

    // Value "get" Functions
    // Return normalization value
    double getNorm() const { return norm; }
    // Return mu value for a specific direction
    double getMu(int) const;
    // Return eta value for a specific direction
    double getEta(int) const;
    // Return xi value for a specific direction
    double getXi(int) const;
    // Return weight value for a specific direction
    double getWt(int) const;
    const vector<double>& get_mu() const { return mu; }
    const vector<double>& get_eta() const { return eta; }
    const vector<double>& get_xi() const { return xi; }
    const vector<double>& get_wt() const { return wt; }
    // Return number of quadrature levels
    int getLevels() const { return levels; }
    // Return number of directions in set
    int getDirs() const { return dirs; }
    // Return number of directions in an octant
    int getOctDirs() const { return odirs; }

    // Vector "get" Functions
    // Returns the omega vector for a specific direction
    vector<double> getOmega(int) const;

    // Print Functions
    void printQuad() const;      // Print quadrature table

    // Compute Functions
    void computeQuadSet(int, double); // Compute the quadrature set

    // Destructor
    ~LSquad();

  private:
    int levels;         // Number of quadrature levels (even values 2-24)
    int dirs;           // Number of quadrature directions
    int odirs;          // Number of quadrature directions in an octant
    double norm;        // Normalization value
    vector<double> mu;  // Computed values of mu
    vector<double> eta; // Computed values of eta
    vector<double> xi;  // Computed values of xi
    vector<double> wt;  // Computed weights

    // Set and Check Functions
    void setLevels(int);         // Set number of quadrature levels
    void setNorm(double);        // Set normalization value
};

} // end namespace rtt_quad

#endif                          // __quad_LSquad_hh__

//---------------------------------------------------------------------------//
//                              end of quad/LSquad.hh
//---------------------------------------------------------------------------//
