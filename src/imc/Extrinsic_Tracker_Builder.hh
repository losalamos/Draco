//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Tracker_Builder.hh
 * \author Mike Buksas
 * \date   Thu Jul 17 13:16:13 2003
 * \brief  Header file for Extrinsic_Tracker_Builder
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Extrinsic_Tracker_Builder_hh
#define rtt_imc_Extrinsic_Tracker_Builder_hh 

#include "ds++/SP.hh"
#include "mc/Surface.hh"
#include "Extrinsic_Surface_Tracker.hh"
#include "Surface_Tracking_Interface.hh"

namespace rtt_mc { class Sphere; class RZWedge_Mesh;}

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Extrinsic_Tracker_Builder
 * \brief This class constructs an Extrinsic_Surface_Tracker
 *
 * This class is responsible for making sure that the geometry of the tally
 * surfaces and the information about mesh cells intersected by these
 * surfaces are consistent. It is constructed with a smart pointer to a class
 * which implements the abstract interface Surface_Tracking_Interface. The
 * derived class must then provide the necessary data to the constructor
 * through this interface.
 *
 */
/*! 
 * \example imc/test/tstExtrinsic_Tracker
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) (Thu Jul 17 13:16:13 2003) Mike Buksas: original
// 
//===========================================================================//

class Extrinsic_Tracker_Builder 
{
  public:

    // NESTED CLASSES AND TYPEDEFS
    
    typedef Surface_Tracking_Interface Interface;

    // CREATORS
    
    //! constructor
    Extrinsic_Tracker_Builder(const rtt_mc::RZWedge_Mesh &mesh,
			      rtt_dsxx::SP<Surface_Tracking_Interface> interface);

    //! destructor
    ~Extrinsic_Tracker_Builder() { /* ... */ }

    // ACCESSORS

    //! Build and return the surface tracker
    rtt_dsxx::SP<Extrinsic_Surface_Tracker> build_tracker();

    //! Get the surface containing status for a cell.
    inline bool get_cell_status (int cell) const;

    //! Get the number of surfaces on this mesh
    int get_local_surfaces() const { return surfaces.size(); }

    //! Get the number of global surfaces
    int get_global_surfaces() const { return global_surface_number; }


  private:

    // DATA
    int global_surface_number;
    int local_surfaces;
    int number_of_cells;

    const rtt_mc::RZWedge_Mesh& mesh;
    rtt_dsxx::SP<Extrinsic_Surface_Tracker> tracker;

    std::vector<rtt_dsxx::SP<rtt_mc::Surface> > surfaces;
    std::vector<int> surface_indices;
    std::vector<bool> surface_in_cell;

    // IMPLEMENTATION

    // General
    void process_surface(const rtt_mc::Surface_Descriptor &descriptor);
    void add_surface_to_list(rtt_dsxx::SP<rtt_mc::Surface> surface);
    bool check_point(const rtt_mc::Surface& surface, double x, double z);

    // Sphere-centric
    void process_sphere(const rtt_mc::Surface_Descriptor &descriptor); 
    bool sphere_intersects_cell(const rtt_mc::Sphere& sphere, int cell);
    bool check_intersections(const rtt_mc::Sphere& sphere);


};

//---------------------------------------------------------------------------//
// Inline functions
//---------------------------------------------------------------------------//

bool Extrinsic_Tracker_Builder::get_cell_status(int cell) const
{
    Check(cell > 0); Check(cell <= number_of_cells);

    return surface_in_cell[cell-1];

}

} // end namespace rtt_imc

#endif // rtt_imc_Extrinsic_Tracker_Builder_hh

//---------------------------------------------------------------------------//
//              end of imc/Extrinsic_Tracker_Builder.hh
//---------------------------------------------------------------------------//
