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
     * \brief Constructor for DummyOpacity object.
     *
     * \sa The DummyOpacity constructor initialized the value of
     *     dummyFilename. 
     */
    DummyOpacity::DummyOpacity()
	: dummyFilename ( "no data file associated with this class" )
	{
	    // empty
	}

    /*!
     * \breif Returns a single opacity value for the prescribed
     *        temperature and density.
     *
     * This opacity object does not do any type of look up or
     * interpolation.  It simply returns a dummy value that is a
     * function of the targetTemperature and the targetDensity.
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    double DummyOpacity::getGrayRosseland(
	const double targetTemperature, 
	const double targetDensity ) 
	{
	    return targetTemperature + targetDensity/10000;
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
     *         are a function of the targetTemperature and the
     *         targetDensity. 
     */    
    vector<double> DummyOpacity::getMGRosseland( 
	const double targetTemperature, 
	const double targetDensity,
	const std::string skey ) 
	{
	    vector<double> dummyMGOpacity(3);
	    for (int i=0; i<3; ++i)
		dummyMGOpacity[i] = (i+1)*1000.0 + targetTemperature + 
		    targetDensity/10000;
	    return dummyMGOpacity;
	}

    /*!
     * \breif Returns a single opacity value for the prescribed
     *        temperature and density.
     *
     * This opacity object does not do any type of look up or
     * interpolation.  It simply returns a dummy value that is a
     * function of the targetTemperature and the targetDensity.
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    double DummyOpacity::getGrayPlank(
	const double targetTemperature, 
	const double targetDensity ) 
	{
	    return targetTemperature + targetDensity/10000;
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
     *         are a function of the targetTemperature and the
     *         targetDensity. 
     */    
    vector<double> DummyOpacity::getMGPlank( 
	const double targetTemperature, 
	const double targetDensity ) 
	{
	    vector<double> dummyMGOpacity(3);
	    for (int i=0; i<3; ++i)
		dummyMGOpacity[i] = (i+1)*1000.0 + targetTemperature + 
		    targetDensity/10000;
	    return dummyMGOpacity;
	}
    
} // end namespace rtt_dummy_opacity


//---------------------------------------------------------------------------//
//                          end of DummyOpacity.cc
//---------------------------------------------------------------------------//
