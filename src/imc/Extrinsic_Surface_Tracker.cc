//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Surface_Tracker.cc
 * \author Mike Buksas
 * \date   Mon Jul 14 16:19:43 2003
 * \brief  Implementation file for Extrinsic_Surface_Tracker
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Extrinsic_Surface_Tracker.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a list of surfaces and surface-in-cell data.
 *
 * This constructor assumes a one-to-one correspondence between local and
 * global surface indices.  In other words, this constructor requires a list
 * of \b all global surfaces in the problem
 * 
 * \param surfaces_ list of all surfaces in the problem
 * \param surface_areas_ list of surface areas for each surface
 * \param surface_in_cell_data_ cell sized list that contains a boolean entry
 * for each cell indicating that one or more surfaces intersect the cell
 */
Extrinsic_Surface_Tracker::Extrinsic_Surface_Tracker(
    const std::vector<SP_Surface> &surfaces,
    const std::vector<double>     &surface_areas_,
    const std::vector<bool>       &surface_in_cell_data_)
    : Surface_tracker(surfaces, surface_areas_),
      surface_in_cell_data(surface_in_cell_data_),
      number_cells(surface_in_cell_data_.size())
{
    
    Check( number_cells > 0 );
    
} 

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a list of surfaces, tally indices, and
 * surface-in-cell data.
 *
 * The surface list includes only those surfaces that intersect at least one
 * cell on the local mesh.  Thus, we require a map of local surface index to
 * global surface index.
 * 
 * \param num_surfaces_ number of global surfaces in the problem
 * \param surfaces_ list of surfaces that are \b local to the mesh
 * \param tally_indices_ map of local surface index to global surface index
 * \param surface_areas_ list of surface areas for each surface
 * \param surface_in_cell_data_ cell sized list that contains a boolean entry
 * for each cell indicating that one or more surfaces intersect the cell
 */
Extrinsic_Surface_Tracker::Extrinsic_Surface_Tracker(
    const int                      num_surfaces_,
    const std::vector<SP_Surface> &surfaces,
    const std::vector<int>        &tally_indices,
    const std::vector<double>     &surface_areas_,
    const std::vector<bool>       &surface_in_cell_data_)
    : Surface_tracker(num_surfaces_, surfaces, tally_indices, surface_areas_),
      surface_in_cell_data(surface_in_cell_data_),
      number_cells(surface_in_cell_data_.size())
{

    Check ( number_cells > 0 );

}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//

void Extrinsic_Surface_Tracker::tally_crossings_implicit_abs(
    const std::vector<double> &position,
    const std::vector<double> &direction,
    int cell,
    double distance,
    double initial_ew,
    double sigma,
    Surface_Sub_Tally& tally)
{

    Check (cell > 0);  Check (cell <= number_cells);

    if ( surface_in_cell_data[cell-1] ) 
    {
	Surface_tracker::tally_crossings_implicit_abs(
	    position, direction, distance, initial_ew, sigma, tally);
    }

}

//---------------------------------------------------------------------------//

void Extrinsic_Surface_Tracker::tally_crossings_analog_abs(
    const std::vector<double> &position,
    const std::vector<double> &direction,
    int cell,
    double distance,
    double ew,
    Surface_Sub_Tally& tally)
{

    Check (cell > 0);  Check (cell <= number_cells);

    if ( surface_in_cell_data[cell-1] ) 
    {
	Surface_tracker::tally_crossings_analog_abs(
	    position, direction, distance, ew, tally);
    }

}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Extrinsic_Surface_Tracker.cc
//---------------------------------------------------------------------------//
