//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Sub_Tally.hh
 * \author Mike Buksas
 * \date   Mon Jun 23 15:33:15 2003
 * \brief  Contains tally information for surface crossings
 *
 * Long description.
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
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
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

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    Surface_Sub_Tally(rtt_dsxx::SP<Azimuthal_Mesh> az_mesh, int surfaces);

    //! copy constructor
    Surface_Sub_Tally(const Surface_Sub_Tally &rhs);

    //! destructor
    ~Surface_Sub_Tally();

    // MANIPULATORS
    
    //! Assignment operator for Surface_Sub_Tally
    Surface_Sub_Tally& operator=(const Surface_Sub_Tally &rhs);

    //! Adds a crossing of given surface with given weight and direction.
    void add_to_tally(int surface, const std::vector<double>& direction, 
		      bool outward, double ew);

    // ACCESSORS

    const std::vector<double>& get_outward_tally(int surface) const;
    const std::vector<double>& get_inward_tally (int surface) const;
    const std::vector<std::vector<double> > get_tally() const { return tally; }

    int get_number_surfaces() const { return surfaces; }
    int get_mesh_size() const { return mesh_size; }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

    rtt_dsxx::SP<Azimuthal_Mesh> the_mesh;
    std::vector<std::vector<double> > tally;
    int mesh_size;   //! number of bins in the azimuthal mesh
    int surfaces;    //! number of surfaces
    int tallies;     //! number of tallies


};

} // end namespace rtt_imc

#endif // rtt_imc_Surface_Sub_Tally_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_Sub_Tally.hh
//---------------------------------------------------------------------------//
