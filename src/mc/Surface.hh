//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Surface.hh
 * \author Mike Buksas
 * \date   Mon Jun 16 16:07:44 2003
 * \brief  Abstract base class for surfaces.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Surface_hh
#define rtt_mc_Surface_hh

#include <vector>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Surface
 * \brief An abstract base classes for surfaces
 *
 * Class Surface is an abstract base class which provides virtual functions
 * required for surface tracking.
 *
 * The constructor, copy operators and destructor are all protected to
 * provect unintended object slicing.
 *
 */
/*! 
 * 
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Surface 
{
  protected:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    Surface() { /* ... */ }

    //! copy constructor
    Surface(const Surface &rhs);

    //! destructor
    virtual ~Surface() { /* ... */ }

  public:

    // MANIPULATORS
    
    //! Assignment operator for Surface
    Surface& operator=(const Surface &rhs);

    // ACCESSORS

    //! Distance from a point, along a direction, to surface
    virtual double distance_to(std::vector<double> position,
			       const std::vector<double>& direction) const = 0;

    //! Distance from a point, with offical status, along a direction, to surface
    virtual double distance_to(std::vector<double> position,
			       const std::vector<double>& direction,
			       bool is_inside) const = 0;


    //! Determine if a point is inside the surface
    virtual bool is_inside(std::vector<double> position) const = 0;

    //! Determine if a point and direction is inside w/ distance_to > 0
    virtual bool is_inside(std::vector<double> position,
			   const std::vector<double> direction) const = 0;

    //! Compute the surface area of the surface
    virtual double surface_area() const = 0;

    //! Compute the volume contained by the surface
    virtual double volume() const = 0;
    
		

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

};

} // end namespace rtt_mc

#endif // rtt_mc_Surface_hh

//---------------------------------------------------------------------------//
//              end of mc/Surface.hh
//---------------------------------------------------------------------------//
