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

Extrinsic_Surface_Tracker::Extrinsic_Surface_Tracker(
    const std::vector<SP_Surface> &surfaces,
    const std::vector<bool>       &surface_in_cell_data_)
    : Surface_tracker(surfaces),
      surface_in_cell_data(surface_in_cell_data_),
      number_cells(surface_in_cell_data_.size())
{
    
    Check( number_cells > 0);
    
} 

Extrinsic_Surface_Tracker::Extrinsic_Surface_Tracker(
    const std::vector<SP_Surface> &surfaces,
    const std::vector<int>        &tally_indices,
    const std::vector<bool>       &surface_in_cell_data_)
    : Surface_tracker(surfaces, tally_indices),
      surface_in_cell_data(surface_in_cell_data_),
      number_cells(surface_in_cell_data_.size())
{

    Check ( number_cells > 0);

}


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
