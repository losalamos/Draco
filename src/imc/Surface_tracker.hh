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

class Surface_Tally;

//===========================================================================//
/*!
 * \class Surface_tracker
 * \brief
 *
 * Surface_tracker detects and tallies occasions when a path segment in
 * Cartesian space crosses any of a collection of abstract surfaces. It also
 * maintains the offical "inside / outside" status of a particle for each
 * surface and uses this status in the correct determination of the streaming
 * distance to each surface.
 *
 * Surfaces are held through a vector of smart pointers to the abstract base
 * class Surface, which are provided in the constructor. This class uses the
 * functions 
 *    Surface::is_inside(vector, const vector&) and
 *    Surface::distance_to(vector, const vector&, bool)
 *
 * The inside/outside status for a new particle is initialized via a call to
 * Surface_tracker::initialize_status(const vector&, const vector&). This
 * call should be performed with the intial position and direction of each
 * particle that surface crossing tallies are desired for.
 *
 * \sa Surface_tracker.cc for detailed descriptions.
 *
/*! 
 * \example imc/test/imc_test.cc 
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) 25 June 1993  First complete version.
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
			 double sigma,
			 Surface_Tally&);

    // ACCESSORS:

    bool get_inside(int i) const { return is_inside[i]; }

  private:

    // DATA

    std::vector<SP_Surface> surfaces;
    std::vector<bool> is_inside;

};

} // end namespace rtt_imc

#endif // rtt_imc_Surface_tracker_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_tracker.hh
//---------------------------------------------------------------------------//
