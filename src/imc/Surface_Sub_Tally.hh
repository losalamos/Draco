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

    //! Copy Constructor
    Surface_Sub_Tally(const Surface_Sub_Tally &rhs);

    //! Destructor
    ~Surface_Sub_Tally();

    // MANIPULATORS
    
    //! Assignment operator for Surface_Sub_Tally
    Surface_Sub_Tally& operator=(const Surface_Sub_Tally &rhs);

    //! Adds a crossing of given surface with given weight and direction.
    void add_to_tally(int surface, const std::vector<double>& direction, 
		      bool outward, double ew);

    // ACCESSORS

    //! Access the outward tally for a given surface #
    const std::vector<double>& get_outward_tally(int surface) const;

    //! Access the inward tally for a given surface #
    const std::vector<double>& get_inward_tally (int surface) const;

    //! Access the entire tally.
    const std::vector<std::vector<double> >& get_tally() const { return tally; }

    //! Get the number of surfaces
    int get_number_surfaces() const { return surfaces; } 
							 
    //! Get the number of bins in the mesh
    int get_mesh_size() const { return mesh_size; }

  private:

    // DATA

    rtt_dsxx::SP<Azimuthal_Mesh> azimuthal_mesh;
    std::vector<std::vector<double> > tally;
    int mesh_size;   //!< number of bins in the azimuthal mesh
    int surfaces;    //!< number of surfaces
    int tallies;     //!< number of tallies


};

} // end namespace rtt_imc

#endif // rtt_imc_Surface_Sub_Tally_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_Sub_Tally.hh
//---------------------------------------------------------------------------//
