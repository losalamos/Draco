//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/Opacity.hh
 * \author Kelly Thompson
 * \date   Thu Jun 23 13:55:06 2000
 * \brief  Opacity header file (an abstract class).
 */
//---------------------------------------------------------------------------//
// $Ids: Opacity.hh,v 1.4 2000/07/24 22:49:20 kellyt Exp $
//---------------------------------------------------------------------------//

#ifndef __cdi_Opacity_hh__
#define __cdi_Opacity_hh__

#include <vector>
#include <string>

namespace rtt_cdi
{
 
//===========================================================================//
/*!
 * \class Opacity
 *
 * \brief This is an abstract base class that CDI opacity objects must be
 *        derived from.
 *
 * \sa Any opacity object that is to be used in conjunction with the CDI
 * must be derived from this abstract base class.  The CDI class must
 * know about the virtual functions found in this class.  Any derived
 * object must define all of these functions.
 *
 */
// revision history:
// -----------------
// 0) original
//===========================================================================//

class Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

  public:

    // CREATORS
    
    // defaulted Opacity() 
    // defaulted Opacity(const Opacity &rhs);
    
    /*!
     * \brief The opacity destructor must be defined here so that we
     *        do not create a memory leak (the virtual opacity class
     *        must be destroyed when the concrete base class is
     *        destroyed). 
     */
    virtual ~Opacity() {}

    // ACCESSORS
    
    /*!
     * \brief Returns the opacity data filename.
     */
    virtual const std::string& getDataFilename() const = 0;

    /*!
     * \brief Returns a single gray Rosseland opacity value for the
     *        prescribed temperature and density.
     *
     * \sa This opacity object only knows how to access the data for one 
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * We do not mark this access const because currently the accessor 
     * getGrayRosseland() in the GandolfOpacity class does modify
     * itself.  It loads the data table from file the first time it is 
     * called.
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    virtual double getGrayRosseland( 
	const double targetTemperature, 
	const double targetDensity ) = 0;

    /*!
     * \brief Returns a vector of the opacity values for each energy
     *        group for the prescribed temperature and density.
     *
     * The opacity object only knows how to access the data for one
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * We do not mark this access const because currently the accessor 
     * getGrayRosseland() in the GandolfOpacity class does modify
     * itself.  It loads the data table from file the first time it is 
     * called.
     *
     * \param targetTemperature The temperature (in keV) of the
     *        material.
     * \param targetDensity The density (in g/cm^3) of the material.
     * \param skey Optional parameter used to specify if the returned
     *        value is scattering ("rsmg"), absorption ("ramg") or
     *        total ("rtmg") opacity.  The default is total.
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector has ngroups entries.  The number of groups is
     *         specified by the data file. 
     */
    virtual std::vector<double> getMGRosseland( 
	const double targetTemperature, 
	const double targetDensity,
	const std::string skey = "rtmg" ) = 0;

    /*!
     * \brief Returns a single gray Plank opacity value for the
     *        prescribed temperature and density.
     *
     * \sa This opacity object only knows how to access the data for one 
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * We do not mark this access const because currently the accessor 
     * getGrayRosseland() in the GandolfOpacity class does modify
     * itself.  It loads the data table from file the first time it is 
     * called.
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    virtual double getGrayPlank( 
	const double targetTemperature, 
	const double targetDensity ) = 0;

    /*!
     * \brief Returns a vector of the opacity values for each energy
     *        group for the prescribed temperature and density.
     *
     * The opacity object only knows how to access the data for one
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * We do not mark this access const because currently the accessor 
     * getGrayPlank() in the GandolfOpacity class does modify
     * itself.  It loads the data table from file the first time it is 
     * called.
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
    virtual std::vector<double> getMGPlank( 
	const double targetTemperature, 
	const double targetDensity ) = 0;

    virtual int getNumTemperatures() const = 0;
    virtual std::vector<double> getTemperatureGrid() const = 0;
    virtual int getNumDensities() const = 0;
    virtual std::vector<double> getDensityGrid() const = 0;
    virtual int getNumGroupBoundaries() const = 0;
    virtual std::vector<double> getGroupBoundaries() const = 0;

    //  protected:     

    // DATA

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif // __cdi_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/Opacity.hh
//---------------------------------------------------------------------------//
