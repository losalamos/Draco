//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/MultigroupOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 14:58:55 2001
 * \brief  MultigroupOpacity class header file (an abstract class)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_MultigroupOpacity_hh__
#define __cdi_MultigroupOpacity_hh__

#include <vector>
#include <string>

#include "OpacityEnums.hh"

namespace rtt_cdi
{

//===========================================================================//
/*!
 * \class MultigroupOpacity
 *
 * \brief This is a pure virtual class that defines a standard
 * interface for all derived MultigroupOpacity objects.
 *
 * \sa  Any derived MultigroupOpacity object must provide as a minimum 
 * the functionality outlined in this routine.  This functionality
 * includes access to the data grid and the ability to return
 * interpolated opacity values.
 */

/*!
 * \example cdi/test/tDummyOpacity.cc
 * \example cdi/test/tCDI.cc
 *
 */
//===========================================================================//

class MultigroupOpacity
{
    // DATA

    // There is no data for a pure virtual object.  This class
    // provides an interface and does not preserve state.

  public:

    // ---------- //
    // Destructor //
    // ---------- //

    /*!
     * \brief Default Opacity() destructor.
     *
     * This is required to correctly release memory when any
     * object derived from MultigroupOpacity is destroyed.
     */
     virtual ~MultigroupOpacity() {};

    // --------- //
    // Accessors //
    // --------- //

    /*!
     * \brief Opacity accessor that returns a vector of opacities (a
     *     single opacity for each group) that correspond to the
     *     provided temperature and density. 
     *
     * \parameter targetTemperature The temperature value for which
     *     these opacity values are being requested (keV).
     *
     * \parameter targetDensity The density value for which these opacity 
     *     values are being requested (g/cm^3)
     *
     * \return A vector of opacities (a single opacity for each group).
     */
    virtual std::vector<double> getOpacity( 
	const double targetTemperature,
	const double targetDensity ) const = 0; 

    /*!
     * \brief Opacity accessor that returns a vector of multigroup
     *     opacity vectors that correspond to the provided vector of
     *     temperatures and a single density value.
     *
     * \parameter targetTemperature A vector of temperature values for
     *     which opacity values are being requested (keV).
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested (g/cm^3).
     *
     * \return A vector of multigroup opacity vectors (cm^2/g).
     */
    virtual std::vector< std::vector<double> > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const = 0; 
    
    /*!
     * \brief Opacity accessor that returns a vector of multigroup
     *     opacity vectors that correspond to the provided vector of
     *     densities and a single temperature value.
     *
     * \parameter targetTemperature The temperature value for which an 
     *     opacity value is being requested (keV).
     *
     * \parameter targetDensity A vector of density values for which
     *     opacity values are being requested (g/cm^3).
     *
     * \return A vector of multigroup opacity vectors (cm^2/g).
     */
    virtual std::vector< std::vector<double> > getOpacity( 
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
     *     opacity data table.
     */
    virtual const std::vector<double>& getTemperatureGrid() const = 0;

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
    virtual const std::vector<double>& getDensityGrid() const = 0;

    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  
     */
    virtual const std::vector<double>& getGroupBoundaries() const = 0;
    
    /*!
     * \brief Returns the size of the temperature grid.
     */
    virtual int getNumTemperatures() const = 0;

    /*! 
     * \brief Returns the size of the density grid.
     */
    virtual int getNumDensities() const = 0;

    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
    virtual int getNumGroupBoundaries() const = 0;

    /*!
     * \brief Returns the number of energy groups 
     * ( getNumGroupBoundaries() - 1 ).
     */
    virtual int getNumGroups() const = 0;

}; // end of class MultigroupOpacity


} // end namespace rtt_cdi

#endif // __cdi_MultigroupOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi/MultigroupOpacity.hh
//---------------------------------------------------------------------------//
