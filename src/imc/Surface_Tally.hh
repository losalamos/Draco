//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Tally.hh
 * \author Mike Buksas
 * \date   Mon Jun 23 15:33:15 2003
 * \brief  Contains tally information for surface crossings
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Surface_Tally_hh
#define rtt_imc_Surface_Tally_hh 

#include "ds++/SP.hh"
#include <vector>

namespace rtt_imc
{

class Azimuthal_Mesh;

//===========================================================================//
/*!
 * \class Surface_Tally
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
 *
 * \sa Surface_Tally.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
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

class Surface_Tally 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    Surface_Tally(rtt_dsxx::SP<Azimuthal_Mesh> az_mesh);

    //! copy constructor
    Surface_Tally(const Surface_Tally &rhs);

    //! destructor
    ~Surface_Tally();

    // MANIPULATORS
    
    //! Assignment operator for Surface_Tally
    Surface_Tally& operator=(const Surface_Tally &rhs);

    void add_to_tally(const std::vector<double>& direction, 
		      bool outward, 
		      double ew);

    // ACCESSORS

    const std::vector<double>& get_outward_tally() { return outward; }
    const std::vector<double>& get_inward_tally() { return inward; }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

    rtt_dsxx::SP<Azimuthal_Mesh> the_mesh;
    int mesh_size;

    std::vector<double> inward;
    std::vector<double> outward;

};

} // end namespace rtt_imc

#endif // rtt_imc_Surface_Tally_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_Tally.hh
//---------------------------------------------------------------------------//
