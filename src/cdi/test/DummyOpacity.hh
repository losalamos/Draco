//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 13 16:11:55 2000
 * \brief  DummyOpacity header file (derived from the Opacity class)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyOpacity_hh__
#define __cdi_DummyOpacity_hh__

#include <vector>
#include <string>

#include "../Opacity.hh"

namespace rtt_dummy_opacity
{

using std::string;
using std::vector;
 
//===========================================================================//
/*!
 * \class DummyOpacity
 *
 * \brief This is a concrete class derived from the cdi/Opacity
 *        class.  It has no real functionality since it is only used
 *        to test the CDI package. 
 *
 * \sa Short Abstract goes here
 */

/*!
 * \example cdi/test/tCDI.cc
 *
 * This test code provides an example of how to use CDI to access an
 * user defined opacity class.  We have created an opacity class
 * called dummyOpacity that is used in the creation of a CDI object.
 * The CDI object is then used to obtain obacity data (via
 * dummyOpacity).
 *
 * The test code also provides a mechanism to test the CDI independent 
 * of any "real" data objects.
 *
 */

// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class DummyOpacity : public rtt_cdi::Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    /*!
     * \brief A list of material identifiers.
     *
     * \sa Normally matIDs would contain a list of material
     * identifiers that specify what materials are available in the
     * data file.  When the concrete opacity class is constructed the
     * user specified material identifier is checked against this list
     * to ensure that the requested material is available.
     */
    vector<int> matIDs; // empty for class DummyOpacity.
    
  public:

    // CREATORS
    
    // defaulted DummyOpacity( );
    // defaulted DummyOpacity(const DummyOpacity &rhs);
    // defaulted ~DummyOpacity();

    // MANIPULATORS
    
    //defaulted DummyOpacity& operator=(const DummyOpacity &rhs);

    // ACCESSORS

    /*!
     * \brief This function returns a dummy string since there is no
     *        data file associated with this class.
     *
     * \sa Normally this function returns the path and name of the
     *     data file that is supplying the opacity data.  For this
     *     dummyOpacity object there is no associated data file so we
     *     return an explanitory string.
     */
    string const getDataFilename() const
    { 
	return "no data file associated with this class"; 
    }

    /*!
     * \brief This function returns a list of material identifiers
     *        that are available in the assocated data file.
     *
     * \sa Again, there is no assocated data file or data so this
     *     function will return an empty vector.
     */
    vector<int> const getMatIDs() const
    {
	return matIDs;
    }

    /*!
     * \breif Returns a single opacity value for the prescribed
     *        temperature and density.
     *
     * This opacity object does not do any type of look up or
     * interpolation.  It simply returns a dummy value.
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    double getGray( const double targetTemperature, 
		    const double targetDensity );
    /*!
     * \breif Returns a vector of the opacity values for each energy
     *        group for the prescribed temperature and density.
     *
     * The opacity object doesn not do any type of look up or
     * interpolation.  It simply returns a vector of dummy values.
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material.
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector for this class has length 3 and all entries are
     *         set equal to zero.
     */
    vector<double> getMG( const double targetTemperature, 
			  const double targetDensity );

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_dummy_opacity

#endif // __cdi_DummyOpacity_hh__

//---------------------------------------------------------------------------//
//             end of cdi/test/DummyOpacity.hh
//---------------------------------------------------------------------------//
