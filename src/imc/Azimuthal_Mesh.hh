//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Azimuthal_Mesh.hh
 * \author Mike Buksas
 * \date   Mon Jun 23 12:56:30 2003
 * \brief  Header file for Azimuthal_Mesh
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Azimuthal_Mesh_hh
#define rtt_imc_Azimuthal_Mesh_hh 

#include <vector>

namespace rtt_imc
{

class Surface_Tracking_Interface;

//===========================================================================//
/*!
 * \class Azimuthal_Mesh
 * \brief Stores azimuthal angular bin information
 *
 * This class implements an azimuthally symmetric angular mesh with it's axis
 * of symmetry oriented along the z axis of a Cartesian coordinate
 * system. The information necessary to describe this mesh is the cosines of
 * the boundaries between the angular bins. 
 *
 * Note that all of the bin arguments in the interface functions are 1-based.
 *
 * \sa Azimuthal_Mesh.cc for detailed descriptions.
 *
 */
/*! 
 * \example imc/test/imc_test.cc 
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Azimuthal_Mesh 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! construct from a vector of cosines
    Azimuthal_Mesh(const std::vector<double>& cosines);

    //! construct from a Surface_Tracking_Interface
    Azimuthal_Mesh(const Surface_Tracking_Interface& interface);

    //! destructor
    ~Azimuthal_Mesh() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Azimuthal_Mesh
    Azimuthal_Mesh& operator=(const Azimuthal_Mesh &rhs);

    // ACCESSORS

    double get_lower_cosine(int bin) const; //!< Cosine of lower boundary angle
    double get_upper_cosine(int bin) const; //!< Cosine of upper boundary angle

    int size() const { return bins; } //!< Get the number of azimuthal bins 

    //! Determine whrether the given direction is in the given bin
    bool is_in_bin(const std::vector<double>& direction, int bin) const;

    //! Find the first bin that contains the given direction.
    int find_bin(const std::vector<double>& direction) const;

  private:

    // DATA
    
    int bins;  //! Number of bins
    std::vector<double> bin_cosines;  //! Cosines of bin boundaries incl -1,1

    // IMPLEMENTATION

    bool check_cosines() const;

};

} // end namespace rtt_imc

#endif // rtt_imc_Azimuthal_Mesh_hh

//---------------------------------------------------------------------------//
//              end of imc/Azimuthal_Mesh.hh
//---------------------------------------------------------------------------//
