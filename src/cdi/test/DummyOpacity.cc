//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 13 16:11:55 2000
 * \brief  Implementation file for DummyOpacity test class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyOpacity.hh"

namespace rtt_dummy_opacity
{

    /*!
     * \breif Returns a single opacity value for the prescribed
     *        temperature and density.
     *
     * This opacity object does not do any type of look up or
     * interpolation.  It simply returns a dummy value (-1.0).
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    double DummyOpacity::getGray( const double targetTemperature, 
				  const double targetDensity )
	{
	    return -1.0;
	}

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
     *         vector for this class has length 3 and the entries are
     *         set equal to { 1.0, 2.0, 3.0 }.
     */    
    vector<double> DummyOpacity::getMG( const double targetTemperature, 
					const double targetDensity ) 
	{
	    vector<double> dummyMGOpacity(3);
	    for (int i=0; i<3; ++i)
		dummyMGOpacity[i] = static_cast<double>(i+1);
	    return dummyMGOpacity;
	}
    
} // end namespace rtt_dummy_opacity


//---------------------------------------------------------------------------//
//                          end of DummyOpacity.cc
//---------------------------------------------------------------------------//
