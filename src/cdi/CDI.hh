//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:06 2000
 * \brief  CDI class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_CDI_hh__
#define __cdi_CDI_hh__

#include <vector>
#include <string>

#include "ds++/SP.hh"

namespace rtt_cdi
{

    // forward declaration (we don't include Opacity.hh in CDI.hh -- but
    // we do in CDI.cc).
    //    class Opacity;

//===========================================================================//
/*!
 * \class CDI
 *
 * \brief This class provides a Common Data Interface (CDI) to Atomic, 
 *        Nuclear and Equation of State (EOS) data.
 *
 * \sa The client must first instantiate concrete Opacity, Nuclear and EOS 
 * classes that are derived from abstrat classes found in the CDI
 * package.  A CDI object is then created using these concrete classes 
 * as constructor parameters.  Each CDI object will provide access to
 * data for one material.
 *
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

class CDI 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    /*!
     * \brief Smart pointer to the opacity object.
     *
     * spOpacity is a smart pointer that links a CDI object to an
     * opacity object (any type of opacity - Gandolf, EOSPAC,
     * Analytic, etc.).  The pointer is established in the CDI
     * constructor. 
     *
     */
    //    const rtt_dsxx::SP<Opacity> spOpacity;
    
    // future expansion
    // const rtt_dsxx::SP<EOS> spEOS;
    
  public:

    // CREATORS
    
    /*!
     * \brief The CDI object instantiates a CDI object by hooking
     *        itself to Opacity, Nuclear, and EOS Data objects.
     *
     * \sa Currently, CDI only interfaces with IPCRESS opacity data
     * (accessed through the GandolfOpacity class).  Because of this
     * only one constructor is currently available.
     *
     * \param _spOpaicty A smart pointer object to an opacity class.
     *                   The opacity class must be derived from the
     *                   abstract class found in the CDI package.
     * \return A CDI object.  A CDI object will be able to access the 
     *         data for a single material
     *
     * We may need to consider keyword arguments in the constructor to 
     * avoid requiring a large number of constructors.  See
     * B. Stroustrup, "The Design and Evolution of C++," Section 6.5.1.
     */
    CDI() {};
    //    CDI( const rtt_dsxx::SP<Opacity> _spOpacity );
    
    // defaulted CDI(const CDI &rhs);

    /*!
     * \brief Destructor for CDI objects.
     *
     * \sa We include a destructor for the CDI class so that, if
     *     another object inherits from CDI, the derived object
     *     correctly destroys the CDI base class.
     */
    virtual ~CDI();


    // MANIPULATORS
    
    // defaulted CDI& operator=(const CDI &rhs);

    // ACCESSORS

    /*!
     * \breif Returns the opacity data filename.
     */
//     virtual std::string getOpacityDataFilename() const;

    /*!
     * \breif Returns a single gray Rosseland opacity value for the
     *        prescribed temperature and density.
     *
     * \sa The opacity object that this CDI object links to only knows how 
     * to access the data for one material.  The material
     * identification is specified in the construction of the opacity
     * object. 
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material. 
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
//     virtual double getGrayRosselandOpacity( 
// 	const double targetTemperature, 
// 	const double targetDensity ) const;

    /*!
     * \breif Returns a vector of Rosseland opacity values
     *        for each energy group for the prescribed temperature and
     *        density. 
     *
     * \sa The opacity object that this CDI object links to only knows how 
     * to access the data for one material.  The material
     * identification is specified in the construction of the opacity
     * object. 
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     * \param skey Optional parameter used to specify if the returned
     *        value is scattering ("rsmg"), absorption ("ramg") or
     *        total ("rtmg") opacity.  The default is total.
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector has ngroups entries.  The number of groups is
     *         specified by the data file. 
     */
//     virtual std::vector<double> getMGRosselandOpacity( 
// 	const double targetTemperature,
// 	const double targetDensity,
//  	const std::string skey = "rtmg" ) const ;

    /*!
     * \breif Returns a single gray Plank opacity value for the
     *        prescribed temperature and density.
     *
     * \sa The opacity object that this CDI object links to only knows how 
     * to access the data for one material.  The material
     * identification is specified in the construction of the opacity
     * object. 
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material. 
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
//     virtual double getGrayPlankOpacity( 
// 	const double targetTemperature, 
// 	const double targetDensity ) const;

    /*!
     * \breif Returns a vector of Plank opacity values
     *        for each energy group for the prescribed temperature and
     *        density. 
     *
     * \sa The opacity object that this CDI object links to only knows how 
     * to access the data for one material.  The material
     * identification is specified in the construction of the opacity
     * object. 
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector has ngroups entries.  The number of groups is
     *         specified by the data file. 
     */
//     virtual std::vector<double> getMGPlankOpacity( 
// 	const double targetTemperature,
// 	const double targetDensity ) const ;

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif                          // __cdi_CDI_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/CDI.hh
//---------------------------------------------------------------------------//
