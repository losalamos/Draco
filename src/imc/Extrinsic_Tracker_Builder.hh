//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Extrinsic_Tracker_Builder.hh
 * \author Mike Buksas
 * \date   Thu Jul 17 13:16:13 2003
 * \brief  
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Extrinsic_Tracker_Builder_hh
#define rtt_imc_Extrinsic_Tracker_Builder_hh 

#include "ds++/SP.hh"
#include "mc/Surface.hh"
#include "Extrinsic_Surface_Tracker.hh"

namespace rtt_mc { class Sphere; class RZWedge_Mesh;}

namespace rtt_imc
{

// Forward declerations:


//===========================================================================//
/*!
 * \class Extrinsic_Tracker_Builder
 * \brief This class constructs an Extrinsic_Surface_Tracker
 *
 * This class is responsible for making sure that the geometry of the tally
 * surfaces and the information about mesh cells intersected by these
 * surfaces are consistent.
 *
 * The interface of this class only supports adding spheres to the list of
 * tally surfaces. Cylinders (if implemented) will require a new public
 * interface member and support for the detection of intersections with cells
 * in the mesh.
 *
 * Once the Extrinsic_Surface_Tracker is created, this class stores a smart
 * pointer to it and avoids creating it again. Adding additional surfaces is
 * not allowed after creating the surface tracker. Thus this class has two
 * states:
 *   - a) "accumulating surfaces" 
 *   - b) "provide surface tracker". 
 *
 * The class is in state a when constructed. Calling the build_tracker()
 * method moves it to state b. It is not possible to switch back to state a.
 *
 */
/*! 
 * \example imc/test/imc_test.cc 
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

    // CREATORS
    
    //! constructor
    Extrinsic_Tracker_Builder(const rtt_mc::RZWedge_Mesh&);

    //! copy constructor (the long doxygen description is in the .cc file)
    Extrinsic_Tracker_Builder(const Extrinsic_Tracker_Builder &rhs);

    //! destructor
    ~Extrinsic_Tracker_Builder() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Extrinsic_Tracker_Builder
    Extrinsic_Tracker_Builder& operator=(const Extrinsic_Tracker_Builder &rhs);
    
    //! Add a sphere to the list of surfaces
    void add_sphere(double z, double r); 

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

    bool sphere_intersects_cell(const rtt_mc::Sphere& sphere, int cell);
    bool check_point(const rtt_mc::Sphere& sphere, double x, double z);
    bool check_intersections(const rtt_mc::Sphere& sphere);
    void add_sphere_to_list(rtt_dsxx::SP<rtt_mc::Sphere> sphere);

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
