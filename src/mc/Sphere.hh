//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphere.hh
 * \author Mike Buksas
 * \date   Mon Jun 16 16:14:46 2003
 * \brief  Implements a spherical surface for surface tallies
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Sphere_hh
#define rtt_mc_Sphere_hh

#include "Surface.hh"
#include <vector>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Sphere
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
 *
 * \sa Sphere.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example mc/test/mc_test.cc 
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Sphere : public Surface
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    Sphere(double center, double radius);

    //! copy constructor
    Sphere(const Sphere &rhs);

    //! destructor
    ~Sphere() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Sphere
    Sphere& operator=(const Sphere &rhs);

    // ACCESSORS

    double distance_to(std::vector<double> position,
		       const std::vector<double>& direction) const;

    double distance_to(std::vector<double> position,
		       const std::vector<double>& direction,
		       bool is_inside) const;

    bool is_inside(std::vector<double> position) const;
    bool is_inside(std::vector<double> position,
		   const std::vector<double> direction) const;

    double surface_area() const;
    double volume() const;

  private:

    // IMPLEMENTATION

    // DATA

    double center;
    double radius, radius_2;


};

} // end namespace rtt_mc

#endif // rtt_mc_Sphere_hh

//---------------------------------------------------------------------------//
//              end of mc/Sphere.hh
//---------------------------------------------------------------------------//
