//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/GrayOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:02:21 2001
 * \brief  GrayOpacity class header file (an abstract class)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_GrayOpacity_hh__
#define __cdi_GrayOpacity_hh__

#include <vector>
#include <string>

#include "OpacityEnums.hh"

namespace rtt_cdi
{

//===========================================================================//
/*!
 * \class GrayOpacity
 *
 * \brief This is a pure virtual class that defines a standard
 *  interface for all derived GrayOpacity objects. 
 *
 * \sa Any derived GrayOpacity object must provide as a minumum the
 * functionality outlined in this routine.  This functionality
 * includes access to the data grid and the ability to return
 * interpolated opacity values.
 */

/*!
 * \example cdi/test/tDummyOpacity.cc
 * \example cdi/test/tCDI.cc
 *
 */
//===========================================================================//

class GrayOpacity
{
    // DATA

    // There is no data for a pure virtual object.  This class
    // provides an interface and does not preserve state.

  public:

    // ---------- //
    // Destructor //
    // ---------- //

    /*!
     * \brief Default GrayOpacity() destructor.
     *
     * This is required to correctly release memory when any
     * object derived from GrayOpacity is destroyed.
     */
    virtual ~GrayOpacity() {};

    // --------- //
    // Accessors //
    // --------- //

    /*!
     * \brief Opacity accessor that returns a single opacity that 
     *     corresponds to the provided temperature and density.
     *
     * \parameter targetTemperature The temperature value for which an
     *     opacity value is being requested (keV).
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested (g/cm^3).
     *
     * \return A single interpolated opacity (cm^2/g).
     */
    virtual double getOpacity( 
	const double targetTemperature,
	const double targetDensity ) const = 0; 

    /*!
     * \brief Opacity accessor that returns a vector of opacities that
     *     correspond to the provided vector of temperatures and a
     *     single density value. 
     *
     * \parameter targetTemperature A vector of temperature values for
     *     which opacity values are being requested (keV).
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested (g/cm^3).
     *
     * \return A vector of opacities (cm^2/g).
     */
    virtual std::vector< double > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const = 0; 
    
    /*!
     * \brief Opacity accessor that returns a vector of opacities
     *     that correspond to the provided vector of densities and a
     *     single temperature value. 
     *
     * \parameter targetTemperature The temperature value for which an 
     *     opacity value is being requested (keV).
     *
     * \parameter targetDensity A vector of density values for which
     *     opacity values are being requested (g/cm^3).
     *
     * \return A vector of opacities (cm^2/g).
     */
    virtual std::vector< double > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity ) const = 0; 

    /*!
     * \brief Returns a string that describes the EnergyPolicy.
     *     Currently this will return either "mg" or "gray."
     */ 
     virtual const std::string& getEnergyPolicyDescriptor() const = 0;

    /*!
     * \brief Returns a "plain English" description of the opacity
     *     data that this class references. (e.g. "Gray Rosseland
     *     Scattering".) 
     */
    virtual const std::string& getDataDescriptor() const = 0;

    /*!
     * \brief Returns the name of the associated data file (if any).
     */
    virtual const std::string& getDataFilename() const = 0;

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table (keV).
     */
    virtual const std::vector<double>& getTemperatureGrid() const = 0;

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table (g/cm^3).
     */
    virtual const std::vector<double>& getDensityGrid() const = 0;

    /*!
     * \brief Returns the size of the temperature grid.
     */
    virtual int getNumTemperatures() const = 0;

    /*! 
     * \brief Returns the size of the density grid.
     */
    virtual int getNumDensities() const = 0;

}; // end of class GrayOpacity


} // end namespace rtt_cdi

#endif // __cdi_GrayOpacity_hh__

//---------------------------------------------------------------------------//
// end of cdi/GrayOpacity.hh
//---------------------------------------------------------------------------//
