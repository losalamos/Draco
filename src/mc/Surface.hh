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
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
 *
 * \sa Surface.cc for detailed descriptions.
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

class Surface 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    Surface() { /* ... */ }

    //! copy constructor
    Surface(const Surface &rhs);

    //! destructor
    virtual ~Surface() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Surface
    Surface& operator=(const Surface &rhs);

    // ACCESSORS

    virtual double distance_to(std::vector<double> position,
			       const std::vector<double>& direction) = 0;

    virtual double distance_to(std::vector<double> position,
			       const std::vector<double>& direction,
			       bool is_inside) = 0;

    virtual bool is_inside(std::vector<double> position) = 0;

    virtual double surface_area() = 0;
    virtual double volume() = 0;
    
		

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
