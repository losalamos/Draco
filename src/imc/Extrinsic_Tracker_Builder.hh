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

#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "mc/Surface.hh"
#include "Extrinsic_Surface_Tracker.hh"
#include "Surface_Tracking_Interface.hh"

// Forward declarations
namespace rtt_mc 
{ 

class Sphere; 
class Surface_Descriptor;
class RZWedge_Mesh;

}

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
 * Extrinsic_Tracker_Builder and Extrinsic_Surface_Tracker test.
 */
// revision history:
// -----------------
// 0) (Thu Jul 17 13:16:13 2003) Mike Buksas: original
// 1) 20-AUG-2003 : added templating on MT; this requires some specialization
//                  in the sphere_intersects_cell function
// 
//===========================================================================//

template<class MT>
class Extrinsic_Tracker_Builder 
{
  public:

    // NESTED CLASSES AND TYPEDEFS
    
    typedef Surface_Tracking_Interface    Interface;
    typedef rtt_dsxx::SP<Interface>       SP_Interface;
    typedef Extrinsic_Surface_Tracker     Tracker;
    typedef rtt_dsxx::SP<Tracker>         SP_Tracker;
    typedef rtt_dsxx::SP<rtt_mc::Surface> SP_Surface;

    // CREATORS
    
    //! Construct from mesh and pointer to interface.
    Extrinsic_Tracker_Builder(const MT &, SP_Interface interface);

    //! Construct from mesh and interface.
    Extrinsic_Tracker_Builder(const MT &, const Interface& interface);

    //! Destructor.
    ~Extrinsic_Tracker_Builder() { /* ... */ }

    // ACCESSORS

    //! Build and return the surface tracker on an MT.
    SP_Tracker build_tracker();

    //! Get the surface containing status for a cell.
    inline bool get_cell_status (int cell) const;

    //! Get the number of surfaces on this mesh
    int get_local_surfaces() const { return surfaces.size(); }

    //! Get the number of global surfaces
    int get_global_surfaces() const { return global_surface_number; }

  private:

    // DATA

    // Const reference to a mesh.
    const MT &mesh;

    // Data used to build surface tracker.
    int                                         number_of_cells;
    int                                         global_surface_number;
    std::vector<rtt_dsxx::SP<rtt_mc::Surface> > surfaces;
    std::vector<int>                            surface_indices;
    std::vector<bool>                           surface_in_cell;

    // IMPLEMENTATION

    // Construction
    void construction_implementation(const Interface& interface);

    // General
    void process_surface(const rtt_mc::Surface_Descriptor &descriptor);
    void add_surface_to_list(SP_Surface surface);
    bool check_point(const rtt_mc::Surface& surface, double x, double z);

    // Sphere-centric
    void process_sphere(const rtt_mc::Surface_Descriptor &descriptor); 
    bool sphere_intersects_cell(const rtt_mc::Sphere& sphere, int cell);
    bool check_intersections(const rtt_mc::Sphere& sphere);
};

//---------------------------------------------------------------------------//
// SPECIALIZED MEMBER FUNCTIONS
// 
// For each mesh type a specialization of sphere_intersects_cell() must be
// provided.
//---------------------------------------------------------------------------//

// Specialization on RZWedge_Mesh.
template<>
bool Extrinsic_Tracker_Builder<rtt_mc::RZWedge_Mesh>::sphere_intersects_cell(
    const rtt_mc::Sphere& sphere, int cell);

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

template<class MT>
bool Extrinsic_Tracker_Builder<MT>::get_cell_status(int cell) const
{
    Check(cell > 0); Check(cell <= number_of_cells);

    return surface_in_cell[cell-1];
}

} // end namespace rtt_imc

#endif // rtt_imc_Extrinsic_Tracker_Builder_hh

//---------------------------------------------------------------------------//
//              end of imc/Extrinsic_Tracker_Builder.hh
//---------------------------------------------------------------------------//
