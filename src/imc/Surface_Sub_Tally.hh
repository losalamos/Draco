//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Sub_Tally.hh
 * \author Mike Buksas
 * \date   Mon Jun 23 15:33:15 2003
 * \brief  Header file for Surface_Sub_Tally
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Surface_Sub_Tally_hh
#define rtt_imc_Surface_Sub_Tally_hh 

#include "ds++/SP.hh"
#include <vector>

namespace rtt_imc
{

class Azimuthal_Mesh;

//===========================================================================//
/*!
 * \class Surface_Sub_Tally
 * \brief Records angular dependence and direction of particle tracings
 *
 * This class accumulates all of the surface crossing information (for a
 * single energy group in a multi-group context). Specifically, the total
 * energy weight, per azimuthal bin, per surface, per direction
 * (inward/outward).
 *
 * The azimuthal bins are defined by a Azimuthal_Mesh object. A smart pointer
 * to the object and the number of surfaces to tally over are requied at
 * construction. 
 *
 * \sa Surface_Sub_Tally.cc for detailed descriptions.
 *
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

class Surface_Sub_Tally 
{
  public:

    // CREATORS
    
    //! Constructor
    Surface_Sub_Tally(rtt_dsxx::SP<Azimuthal_Mesh> az_mesh, int surfaces);

    //! Destructor
    ~Surface_Sub_Tally();

    // MANIPULATORS
    
    //! Adds a crossing of given surface with given weight and direction.
    void add_to_tally(int surface, const std::vector<double>& direction, 
		      bool outward, double ew);

    // ACCESSORS

    //! Access the outward weight tally for a given surface #
    const std::vector<double>& get_outward_weight_tally(int surface) const;

    //! Access the inward weight tally for a given surface #
    const std::vector<double>& get_inward_weight_tally (int surface) const;

    //! Access the outward count tally for a given surface #
    const std::vector<int>& get_outward_count_tally(int surface) const;

    //! Access the inward count tally for a given surface #
    const std::vector<int>& get_inward_count_tally (int surface) const;

    //! Access the entire weight tally.
    inline const std::vector<std::vector<double> >& get_weight_tally() const; 

    //! Access the entire count tally.
    inline const std::vector<std::vector<int> >& get_count_tally() const; 

    //! Access the weight tallied for a surface and bin
    inline double weight(int surface, bool outward, int bin) const;

    //! Access the crossing countfor a surface and bin
    inline int crossings(int surface, bool outward, int bin) const;

    //! Get the number of surfaces
    int get_number_surfaces() const { return surfaces; } 
							 
    //! Get the number of bins in the mesh
    int get_mesh_size() const { return mesh_size; }

  private:

    // DATA

    rtt_dsxx::SP<Azimuthal_Mesh>      azimuthal_mesh;
    std::vector<std::vector<double> > weight_tally;
    std::vector<std::vector<int>    > count_tally;
    int mesh_size;   //!< number of bins in the azimuthal mesh
    int surfaces;    //!< number of surfaces
    int tallies;     //!< number of tallies

    // Implementation

    inline int get_surface_index(int surface, bool is_outward) const;


};


//---------------------------------------------------------------------------//
// Inline accessors:
//---------------------------------------------------------------------------//

const std::vector<std::vector<double> >& Surface_Sub_Tally::get_weight_tally() 
    const
{ return weight_tally; }

const std::vector<std::vector<int> >& Surface_Sub_Tally::get_count_tally() 
    const
{ return count_tally; }

//---------------------------------------------------------------------------//
double Surface_Sub_Tally::weight(int surface, bool outward, int bin) const
{

    Check (bin > 0); Check(bin <= mesh_size);
    int surface_index = get_surface_index(surface, outward);
    int bin_index = bin - 1;

    return weight_tally[surface_index][bin_index];
}

//---------------------------------------------------------------------------//
int Surface_Sub_Tally::crossings(int surface, bool outward, int bin) const
{

    Check (bin > 0); Check(bin <= mesh_size);
    int surface_index = get_surface_index(surface, outward);
    int bin_index = bin - 1;

    return count_tally[surface_index][bin_index];
}

//---------------------------------------------------------------------------//
// Inline implementation
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Takes a 1-based surface index and outward bool and computes
 * internal surface index
 * 
 * \param surface 1-based surface number
 * \param outward boolean true if outward direction is desires
 * \return internally used index for the tally
 */
int Surface_Sub_Tally::get_surface_index(int surface, bool is_outward) const
{
    Check(surface > 0); Check(surface <= surfaces);
    int surface_index = 2 * (surface - 1) + static_cast<int>(is_outward);
    Ensure(surface_index >= 0); Ensure(surface_index < tallies);
    return surface_index;
}


} // end namespace rtt_imc

#endif // rtt_imc_Surface_Sub_Tally_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_Sub_Tally.hh
//---------------------------------------------------------------------------//
