//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Surface_Tracker.hh
 * \author Mike Buksas
 * \date   Mon Jul 14 16:19:43 2003
 * \brief  Header file for Extrinsic_Surface_Tracker
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Extrinsic_Surface_Tracker_hh
#define rtt_imc_Extrinsic_Surface_Tracker_hh

#include "Surface_tracker.hh"
#include "mc/Surface.hh"
#include "ds++/SP.hh"

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Extrinsic_Surface_Tracker
 * \brief Implements surface traking vis-a-vis \c Surface_tracker in the
 * presence of an RZ wedge mesh.
 *
 * The surfaces referred to here are extrinsic in the sense of being overlaid
 * onto a mesh. This class keeps (and requires at construction) a boolean
 * vector of the local mesh size (synonymous with the global size for
 * replicated meshes). True indicates that the corresponding cell in the mesh
 * is intersected by at least one of the surfaces.
 *
 * The surface-in-cell data serves two purposes: It is provided as a service
 * to the transport algorithm where it may affect the implementation of
 * hybrid imc/diffusion methods. It is used internally to short-circut the
 * tally crossings algorithms. 
 *
 */
// revision history:
// -----------------
// 0) (Mon Jul 14 16:19:43 2003) Mike Buksas: original
// 1) 05-JAN-2004: added surface areas
// 
//===========================================================================//

class Extrinsic_Surface_Tracker : public Surface_tracker
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_dsxx::SP<rtt_mc::Surface> SP_Surface;
    typedef std::vector<SP_Surface>::const_iterator surface_iterator;
    typedef std::vector<bool>::iterator bool_iterator;

    // Full constructor.
    Extrinsic_Surface_Tracker(const int num_global_surfaces,
			      const std::vector<SP_Surface> &surfaces,
			      const std::vector<int> &tally_indices,
			      const std::vector<double> &surface_areas,
			      const std::vector<bool> &surface_in_cell_data);

    // Construct with default indices.
    Extrinsic_Surface_Tracker(const std::vector<SP_Surface> &surfaces,
			      const std::vector<double> &surface_areas,
			      const std::vector<bool> &surface_in_cell_data);

    //! Copy constructor (the long doxygen description is in the .cc file).
    Extrinsic_Surface_Tracker(const Extrinsic_Surface_Tracker &rhs);

    //! Destructor.
    ~Extrinsic_Surface_Tracker() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Extrinsic_Surface_Tracker.
    Extrinsic_Surface_Tracker& operator=(const Extrinsic_Surface_Tracker& rhs);

    // ACCESSORS

    inline bool surface_in_cell(int cell);

    //! Process implicit streamings across surfaces in a particlar cell.
    void tally_crossings_implicit_abs(const std::vector<double> &position,
				      const std::vector<double> &direction,
				      int cell,
				      double distance,
				      double initial_ew,
				      double sigma,
				      Surface_Sub_Tally&);

    //! Process analog streamings across surfaces in a particlar cell.
    void tally_crossings_analog_abs(const std::vector<double> &position,
				    const std::vector<double> &direction,
				    int cell,
				    double distance,
				    double ew,
				    Surface_Sub_Tally&);
  private:

    // DATA

    std::vector<bool> surface_in_cell_data;
    int number_cells;

};

//---------------------------------------------------------------------------//
// Inline functions
//---------------------------------------------------------------------------//
/*!
 * \brief Return true if there is a surface in the given cell on the local
 * mesh.
 */
bool Extrinsic_Surface_Tracker::surface_in_cell(int cell)
{
    Check (cell >0 );  Check (cell <= number_cells);
    return surface_in_cell_data[cell-1];
}

} // end namespace rtt_imc

#endif // rtt_imc_Extrinsic_Surface_Tracker_hh

//---------------------------------------------------------------------------//
//              end of imc/Extrinsic_Surface_Tracker.hh
//---------------------------------------------------------------------------//
