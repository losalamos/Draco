//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Azimuthal_Mesh.hh
 * \author Mike Buksas
 * \date   Mon Jun 23 12:56:30 2003
 * \brief  Class representing a mesh of aximuthal bins centered on the z axis.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Azimuthal_Mesh_hh
#define rtt_imc_Azimuthal_Mesh_hh 

#include <vector>

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Azimuthal_Mesh
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at:
 *
 * Doxygen tutorial: http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * Doxygen keywords: http://www.stack.nl/~dimitri/doxygen/commands.html
 *
 * \sa Azimuthal_Mesh.cc for detailed descriptions.
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

class Azimuthal_Mesh 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! default constructors
    Azimuthal_Mesh(const std::vector<double>& cosines);

    //! copy constructor
    Azimuthal_Mesh(const Azimuthal_Mesh &rhs);

    //! destructor
    ~Azimuthal_Mesh() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Azimuthal_Mesh
    Azimuthal_Mesh& operator=(const Azimuthal_Mesh &rhs);

    // ACCESSORS

    double get_lower_cosine(int bin) const;
    double get_upper_cosine(int bin) const;

    int size() const { return bins; }
    bool is_in_bin(const std::vector<double>& direction, int bin) const;
    int find_bin(const std::vector<double>& direction) const;

  private:

    // DATA
    
    int bins;
    std::vector<double> bin_cosines;

    // IMPLEMENTATION

    bool check_cosines(const std::vector<double>& cosines) const;

};

} // end namespace rtt_imc

#endif // rtt_imc_Azimuthal_Mesh_hh

//---------------------------------------------------------------------------//
//              end of imc/Azimuthal_Mesh.hh
//---------------------------------------------------------------------------//
