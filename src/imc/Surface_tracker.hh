//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_tracker.hh
 * \author Mike Buksas
 * \date   Thu Jun 19 11:33:00 2003
 * \brief  Header file for Surface_tracker
 * \note   Copyright © 2003 The Regents of the University of California.
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

class Surface_Sub_Tally;

//===========================================================================//
/*!
 * \class Surface_tracker
 * \brief Detects particle paths which cross a collection of surfaces.
 *
 * Surface_tracker detects and tallies occasions when a path segment in
 * Cartesian space crosses any of a collection of abstract surfaces. It also
 * maintains the offical "inside / outside" status of a particle for each
 * surface and uses this status in the correct determination of the streaming
 * distance to each surface.
 *
 * Surfaces are held through a vector of smart pointers to the abstract base
 * class Surface, which are provided in the constructor. This class uses the
 * functions:
 * \code   
 *    Surface::is_inside(vector position, const vector& direction)
 *    Surface::distance_to(vector position, const vector& direction, bool is_inside)
 * \endcode
 *
 * The inside/outside status for a new particle is initialized via a call to
 * Surface_tracker::initialize_status(const vector&, const vector&). This
 * call should be performed with the intial position and direction of each
 * particle that surface crossing tallies are desired for.
 *
 * The surface tracker also maintains a list of global surface indices for
 * it's surfaces. It will generally hold fewer than the total number of
 * surfaces on domain-decomosed meshes. We also allow the surface tracker to
 * hold no surfaces. 
 *
 * \sa Surface_tracker.cc for detailed descriptions.
 */
/*! 
 * \example imc/test/tstSurface_tracker.cc 
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
    typedef std::vector<bool>::iterator bool_iterator;

    // CREATORS
    
    //! Construct with a list of surfaces
    Surface_tracker(const std::vector<SP_Surface>& surfaces);

    //! Construct with a list of surfaces and tally indices
    Surface_tracker(const std::vector<SP_Surface>& surfaces,
		    const std::vector<int>& tally_indices);

    //! copy constructor
    Surface_tracker(const Surface_tracker &rhs);

    //! destructor
    ~Surface_tracker() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Surface_tracker
    Surface_tracker& operator=(const Surface_tracker &rhs);

    void initialize_status(const std::vector<double> &position,
			   const std::vector<double> &direction);

    void tally_crossings_implicit_abs(const std::vector<double> &position,
				      const std::vector<double> &direction,
				      double distance,
				      double initial_ew,
				      double sigma,
				      Surface_Sub_Tally&);

    void tally_crossings_analog_abs(const std::vector<double> &position,
				    const std::vector<double> &direction,
				    double distance,
				    double ew,
				    Surface_Sub_Tally&);

    // ACCESSORS:

    inline bool get_inside(int surface) const;

    bool get_surface_in_cell(int cell) const;

    int surfaces() const { return surface_list.size(); }

  private:

    // DATA

    std::vector<SP_Surface> surface_list;
    std::vector<int> tally_indices;
    std::vector<bool> is_inside;

    // IMPLEMENTATION

    inline int get_data_index(int surface) const;

};

//---------------------------------------------------------------------------//
// Inline functions
//---------------------------------------------------------------------------//

bool Surface_tracker::get_inside(int surface) const
{
    int index = get_data_index(surface);
    
    Require ( index >= 0 );
    Require ( index < is_inside.size() );

    return is_inside[index];

}

//---------------------------------------------------------------------------//
int Surface_tracker::get_data_index(int surface) const
{
    Check(surface > 0);

    // Find the surface index in the vector of indices:
    std::vector<int>::const_iterator surface_iterator = 
	std::find(tally_indices.begin(), tally_indices.end(), surface);

    Require ( surface_iterator != tally_indices.end());

    return static_cast<int>(surface_iterator - tally_indices.begin());

}
    


} // end namespace rtt_imc

#endif // rtt_imc_Surface_tracker_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_tracker.hh
//---------------------------------------------------------------------------//
