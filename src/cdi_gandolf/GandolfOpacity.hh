//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  GandolfOpacity class header file (derived from cdi/Opacity)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfOpacity_hh__
#define __cdi_gandolf_GandolfOpacity_hh__

#include <vector>
#include <string>
#include <cmath>

#include "ds++/SP.hh"
#include "cdi/Opacity.hh"

namespace rtt_cdi_gandolf
{

    // forward declaration (we don't include GandolfFile.hh in this
    // header -- but we do in GandolfFile.cc.
    class GandolfFile;

//===========================================================================//
/*!
 * \class GandolfOpacity
 *
 * \brief This is a concrete class derived from cdi/Opacity.  This
 *        class allows to client to access the data in IPCRESS files
 *        via the Gandolf libraries.
 *
 * \sa This class is designed to be used in conjuction with the CDI.
 * The client code will create a GandolfOpacity object and use this
 * object as an argument during the CDI instantiation.  The purpose of 
 * this class is to provide a mechanism for accessing data in IPCRESS
 * files and works by calling the Gandolf library provided by X-5.
 * The GandolfOpacity constructor expects two arguments: the IPCRESS
 * file name and a material identifier.  Once constructed this object
 * allows the client to access any data found in the IPCRESS file for
 * that one material.  The client code will need to create a separate
 * GandolfOpacity object for each material that it need information
 * about. 
 *
 */

/*!
 * \example cdi_gandolf/test/tCDIGandolf.cc
 *
 * Example of GandolfOpacity usage independent of CDI.  In this
 * example we construct a GandolfOpacity object for the material
 * Aluminum (matID=10001 in our example IPCRESS file).  We then use
 * the GandolfOpacity object to compute a Rosseland Gray opacity
 * value for a specified material temperature and density.  In a
 * similar fashion we also request the GandolfOpacity object to return 
 * a vector of multigroup opacities for a specified temperature and
 * density. 
 *
 */

// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class GandolfOpacity : public rtt_cdi::Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    /*!
     * \brief GandolfFile object identifies the physical data file.
     */
    const rtt_dsxx::SP<GandolfFile> spGandolfFile;

    /*!
     * \brief Number of materials found in the data file.
     */
    int numMaterials;

    /*!
     * \brief The material ID that this instance of GandolfOpacity 
     *        represents.
     */
    const int matID;

    /*!
     * \brief List of keywords for current material found in the data
     *        file. 
     */
    std::vector<std::string> vkeys;
    /*!
     * \brief Number of keys available for the current material found
     *        in the data file.
     */
    int numKeys; 

    /*!
     * \breif Number of temperatures found in the data file.
     */
    int numTemps;

    /*!
     * \breif Number of densities found in the data file. 
     */
    int numDensities;

    /*!
     * \breif Number of energy group boundaries found in the data
     *        file. 
     *        ( numGroups = numGroupBoundaries - 1 )
     */
    int numGroupBoundaries; 

    /*!
     * \breif Number of gray opacities found in the data file.
     *        ( = numTemps * numDensities )
     */
    int numGrayOpacities; 

    /*!
     * \breif Number of multigroup intensities found in the data file
     *        ( = numTemps * numDensities * numGroups )
     */
    int numMGOpacities;
    
    /*!
     * \brief The log(temperature) grid found in the data file.  
     */
    std::vector<double> logTemperatures;

    /*!
     * \brief The log(density) grid found in the data file.
     */
    std::vector<double> logDensities;

    /*!
     * \brief The energy group boundaries found in the data file.
     */
    std::vector<double> groupBoundaries;

    /*!
     * \brief The log(gray opacity) table found in the data file.
     */
    std::vector<double> logGrayOpacities;

    /*!
     * \brief The log(multigroup opacity) table found in the data file.
     */
    std::vector<double> logMGOpacities;

    /*!
     * \brief Have we loaded the gray Rosseland opacity table yet?
     */
    bool grayRosselandTableLoaded;

    /*!
     * \brief Have we loaded the multigroup Rosseland opacity table yet?
     */
    bool mgRosselandTableLoaded;

    /*!
     * \brief Have we loaded the gray Plank opacity table yet?
     */
    bool grayPlankTableLoaded;

    /*!
     * \brief Have we loaded the multigroup Plank opacity table yet?
     */
    bool mgPlankTableLoaded;

  public:

    // CREATORS
    
    /*!
     * \brief Standard GandolfOpacity constructor.
     * 
     * \sa This is the standard GandolfOpacity constructor.  The
     * GandolfFile object must refer to an existing object and the
     * material identifier must exist in the data file.  This object
     * is usually instantiated as a smart pointer (especially if it is
     * to be used in conjunction with the CDI class). 
     *
     * example constructor:
     *
     * rtt_dsxx::SP spGF 
     *     = new rtt_cdi_gandolf::GandolfFile( filename );
     *
     * rtt_dsxx::SP spGanOp
     *     = new rtt_cdi_gandolf::GandolfOpacity( spGF, matid )
     *
     * \param _spGandolfFile A smart pointer to a GandolfFile object.
     *     This object links the data file to the GandolfOpacity object.
     *
     * \param _matid The material identifier is a 5 digit integer that 
     *     the client must specify.  This material identifier must
     *     exist in the IPCRESS file. 
     */
    GandolfOpacity( const rtt_dsxx::SP<GandolfFile> _spGandolfFile, 
		    const int _matid );

    // (defaulted) GandolfOpacity(const GandolfOpacity &rhs);
    // (defaulte) ~GandolfOpacity()

    // MANIPULATORS
    
    //defaulted GandolfOpacity& operator=(const GandolfOpacity &rhs);

    // ACCESSORS

    /*!
     * \brief Returns the IPCRESS data filename.
     */
    const std::string& getDataFilename() const;

    /*!
     * \brief Returns opacity value(s) for the prescribed
     *        temperature(s) and density(ies).
     *
     * This opacity object only knows how to access the data for one 
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * Additionally, the default behavior of this function is to
     * return a Rosseland total (scatter plus absorption) opacity
     * value.   
     *
     * \param targetTemperature The temperature(s) (in keV) of the 
     *        material. Expecting a double or a vector<double>.
     * \param targetDensity The density(ies) (in g/cm^3) of the
     *        material.  Expecting a double or a vector<double>.
     * \param skey Optional parameter used to specify if the returned
     *        value is scattering ("rsgray"), absorption ("ragray") or
     *        total ("rgray") opacity.  The the accessor returns the
     *        total opacity ("rgray") by default.
     * \return Gray opacity value(s) for the current material at
     *         targetTemperature(s) keV and targetDensity(ies)
     *         g/cm^3.  If a vector was used for targetTemperature or
     *         targetDensity then a vector of opacities will be
     *         returned.  If both targetTemperature and targetDensity
     *         are vectors then the resulting opacity vector has index 
     *         ordering (it*nt+id) where nt is
     *         targetTemperature.size(), it is the temperature index
     *         and id is the density index.
     */
    // Currenlty this accessor modifies the GandolfOpacity object so
    // it is not declared as const.  The opacity accessors load the
    // data when they are first called.  After the first call to this
    // accessor the tabulated data is retrieved from the
    // GandolfOpacity object and not the data file.
    double getGrayRosseland( 
	const double targetTemperature, 
	const double targetDensity,
	const std::string skey );
    std::vector<double> getGrayRosseland(
	const std::vector<double> targetTemperatures,
	const double targetDensity,
	const std::string skey );
    std::vector<double> getGrayRosseland(
	const double targetTemperatures,
	const std::vector<double> targetDensity,
	const std::string skey );
    std::vector<double> getGrayRosseland(
	const std::vector<double> targetTemperatures,
	const std::vector<double> targetDensity,
 	const std::string skey );

    /*!
     * \brief Returns a vector of the opacity values for each energy
     *        group for the prescribed temperature and density.
     *
     * The opacity object only knows how to access the data for one
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * Additionally, the default behavior of this function is to
     * return a Rosseland opacity value.  
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material.
     * \param targetDensity The density (in g/cm^3) of the material.
     * \param skey Optional parameter used to specify if the returned
     *        value is scattering ("rsmg"), absorption ("ramg") or
     *        total ("rtmg") opacity.  The default is total.
     *
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector has ngroups entries.  The number of groups is
     *         specified by the data file. 
     */

    // Currenlty this accessor modifies the GandolfOpacity object so
    // it is not declared as const.  The opacity accessors load the
    // data when they are first called.  After the first call to this
    // accessor the tabulated data is retrieved from the
    // GandolfOpacity object and not the data file.

    // The default value for "skey" is set in cdi/Opacity.hh.
    
    std::vector<double> getMGRosseland( 
	const double targetTemperature, 
	const double targetDensity,
	const std::string skey );
    std::vector<double> getMGRosseland( 
	const std::vector<double> targetTemperature, 
	const double targetDensity,
	const std::string skey );
    std::vector<double> getMGRosseland( 
	const double targetTemperature, 
	const std::vector<double> targetDensity,
	const std::string skey );
    std::vector<double> getMGRosseland( 
	const std::vector<double> targetTemperature, 
	const std::vector<double> targetDensity,
	const std::string skey );
 
    /*!
     * \brief Returns opacity value(s) for the prescribed
     *        temperature(s) and density(ies).
     *
     * This opacity object only knows how to access the data for one 
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * Additionally, the default behavior of this function is to
     * return a Plank total opacity value (scatter plus absorption)
     * opacity value.
     *
     * \param targetTemperature The temperature(s) (in keV) of the 
     *        material. Expecting a double or a vector<double>.
     * \param targetDensity The density(ies) (in g/cm^3) of the
     *        material.  Expecting a double or a vector<double>.
     * \param skey Optional parameter used to specify if the returned
     *        value is scattering ("psgray"), absorption ("pagray") or
     *        total ("pgray") opacity.  The the accessor returns the
     *        total opacity ("pgray") by default.
     * \return Gray opacity value(s) for the current material at
     *         targetTemperature(s) keV and targetDensity(ies)
     *         g/cm^3.  If a vector was used for targetTemperature or
     *         targetDensity then a vector of opacities will be
     *         returned.  If both targetTemperature and targetDensity
     *         are vectors then the resulting opacity vector has index 
     *         ordering (it*nt+id) where nt is
     *         targetTemperature.size(), it is the temperature index
     *         and id is the density index.
     */
    // Currenlty this accessor modifies the GandolfOpacity object so
    // it is not declared as const.  The opacity accessors load the
    // data when they are first called.  After the first call to this
    // accessor the tabulated data is retrieved from the
    // GandolfOpacity object and not the data file.
    double getGrayPlank( 
	const double targetTemperature, 
	const double targetDensity,
	const std::string skey );
    std::vector<double> getGrayPlank( 
	const std::vector<double> targetTemperature,
	const double targetDensity,
	const std::string skey );
    std::vector<double> getGrayPlank( 
	const double targetTemperature,
	const std::vector<double> targetDensity,
	const std::string skey );
    std::vector<double> getGrayPlank( 
	const std::vector<double> targetTemperature,
	const std::vector<double> targetDensity,
	const std::string skey );

    /*!
     * \brief Returns a vector of the opacity values for each energy
     *        group for the prescribed temperature and density.
     *
     * The opacity object only knows how to access the data for one
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * Additionally, the default behavior of this function is to
     * return a Plank opacity value.  
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
    // Currenlty this accessor modifies the GandolfOpacity object so
    // it is not declared as const.  The opacity accessors load the
    // data when they are first called.  After the first call to this
    // accessor the tabulated data is retrieved from the
    // GandolfOpacity object and not the data file.
    std::vector<double> getMGPlank( 
	const double targetTemperature, 
	const double targetDensity,
	const std::string skey );
    std::vector<double> getMGPlank( 
	const std::vector<double> targetTemperature, 
	const double targetDensity,
	const std::string skey );
    std::vector<double> getMGPlank( 
	const double targetTemperature, 
	const std::vector<double> targetDensity,
	const std::string skey );
    std::vector<double> getMGPlank( 
	const std::vector<double> targetTemperature, 
	const std::vector<double> targetDensity,
	const std::string skey );

    /*!
     * \breif Returns the number of temperature points in the data
     *        grid.
     */
    int getNumTemperatures() const { return numTemps; };

    /*!
     * \brief Returns a list of temperature (keV) that make up the
     *        temperature grid in the IPCRESS file.
     */
    std::vector<double> getTemperatureGrid() const;

    /*!
     * \breif Returns the number of density points in the data
     *        grid.
     */
    int getNumDensities() const { return numDensities; };

    /*!
     * \brief Returns a list of densities (g/cm^3) that make up the 
     *        density grid in the IPCRESS file.
     */
    std::vector<double> getDensityGrid() const;

    /*!
     * \breif Returns the number of group boundary points in the data
     *        grid.
     */
    int getNumGroupBoundaries() const { return numGroupBoundaries; };

    /*!
     * \brief Returns a list of Group Boundaries (keV) that make up the 
     *        data grid in the IPCRESS file.
     */
    std::vector<double> getGroupBoundaries() const 
    { 
	return groupBoundaries; 
    };

  private:
    
    // IMPLEMENTATION

    /*!
     * \breif Generic routine used to retrieve gray opacity data.
     *
     * \param skey This string identifies the model for the gray
     *        opacity data.  Valid values are "Rosseland" or "Plank". 
     * \param grayTableLoaded If the gray table has alredy been loaded 
     *        then we don't need to access the data on disk.
     * \param otherTableLoaded If we are loading a new grayTable on 
     *        top of another table this value will be "true".  If this 
     *        is the case, our routine verifies that the temperature
     *        and density grids for the this table are identical to
     *        those of the previously loaded table.
     * \param targetTemperature The temperature (in keV) of the
     *        material.
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray Opacity value corresponding to targetTemperature
     *        and targetDensity using model skey.
     */
    double getGray( const std::string &skey,
		    bool &grayTableLoaded,
		    const double targetTemperature,
		    const double targetDensity );

    /*!
     * \breif Generic routine used to retrieve multigroup opacity data.
     *
     * \param skey This string identifies the model for the gray
     *        opacity data.  Valid values are "Rosseland" or "Plank". 
     * \param mgTableLoaded If the multigroup table has alredy been loaded 
     *        then we don't need to access the data on disk.
     * \param otherTableLoaded If we are loading a new mgTable on 
     *        top of another table this value will be "true".  If this 
     *        is the case, our routine verifies that the temperature
     *        and density grids for the this table are identical to
     *        those of the previously loaded table.
     * \param targetTemperature The temperature (in keV) of the
     *        material.
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return A vector of opacity values corresponding to
     *         targetTemperature and targetDensity using model skey.
     *         This vector will be numGroups long.
     */
    std::vector<double> getMG( const std::string &skey,
			       bool &mgTableLoaded,
			       const double targetTemperature,
			       const double targetDensity );
    
    // These two functions are delcared "static" because we only need
    // one copy of these functions -- not one copy per instance of
    // GandolfOpacity. 

    /*! 
     * \brief This function returns "true" if "key" is found in the list
     *        of "keys".
     */
    template < typename T >
    static bool key_available( T key, std::vector<T> keys ); 
    
    /*!
     * \brief This function compares two double precision vectors.  If 
     *        the vectors are equal within some tolerance then the
     *        function returns "true".
     */
    static bool isSame( const std::vector<double> &v1, 
			const std::vector<double> &v2 );

    /*!
     * \brief Has any table been loaded?
     */
    bool anyTableLoaded() const { return grayRosselandTableLoaded &&
				      mgRosselandTableLoaded &&
				      grayPlankTableLoaded &&
				      mgPlankTableLoaded; };

}; // end of class GandolfOpacity

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfOpacity.hh
//---------------------------------------------------------------------------//
