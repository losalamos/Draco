//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_tracker.hh
 * \author Mike Buksas
 * \date   Thu Jun 19 11:33:00 2003
 * \brief  Computes and tallies surface crossings for a collection of surfaces.  
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Surface_tracker_hh
#define rtt_imc_Surface_tracker_hh

#include <vector>

#include "mc/Surface.hh"
#include "ds++/SP.hh"


namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Surface_tracker
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
 *
 * \sa Surface_tracker.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
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

class Surface_tracker 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_dsxx::SP<rtt_mc::Surface> SP_Surface;
    typedef std::vector<SP_Surface>::const_iterator surface_iterator;

    // CREATORS
    
    //! default constructors
    Surface_tracker(const std::vector<SP_Surface>& surfaces);

    //! copy constructor
    Surface_tracker(const Surface_tracker &rhs);

    //! destructor
    ~Surface_tracker() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Surface_tracker
    Surface_tracker& operator=(const Surface_tracker &rhs);

    void initialize_status(const std::vector<double>& position,
			   const std::vector<double>& direction);

    void tally_crossings(const std::vector<double>& position,
			 const std::vector<double>& direction,
			 double distance,
			 double initial_ew,
			 double sigma);

    // ACCESSORS:

    bool get_inside(int i) const { return is_inside[i]; }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

    std::vector<SP_Surface> surfaces;
    std::vector<bool> is_inside;

};

} // end namespace rtt_imc

#endif // rtt_imc_Surface_tracker_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_tracker.hh
//---------------------------------------------------------------------------//
