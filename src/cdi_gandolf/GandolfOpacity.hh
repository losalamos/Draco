//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  Gandolf opacity header file (derived from cdi/Opacity)
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

  public:

    // CREATORS
    
    /*!
     * \brief Standard GandolfOpacity constructor.
     * 
     * \sa This is the standard GandolfOpacity constructor.  The
     * filename must refer to an IPCRESS data file and the material
     * identifier must exist in the data file.  This object is usually 
     * instantiated as a smart pointer (especially if it is to be used 
     * in conjunction with the CDI class).  
     *
     * rtt_dsxx::SP spGanOpAl;
     * spGanOpAl = new 
     *     Gandolf rtt_cdi_gandolf::GandolfOpacity( data_file, matid )
     *
     * \param _data_filename The name of the IPCRESS data file.  While 
     *                       this class can use any filename the F77
     *                       Gandolf library requires the name to have
     *                       80 characters or less. 
     *
     * \param _matid The material identifier is a 5 digit integer that 
     *               the client must specify.  This material
     *               identifier must exist in the IPCRESS file.
     */
//     GandolfOpacity( const std::string& _data_filename, 
// 		    const int _matid );

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
     * \breif Returns a single opacity value for the prescribed
     *        temperature and density.
     *
     * This opacity object only knows how to access the data for one 
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * Additionally, the default behavior of this function is to
     * return a Rosseland opacity value.  
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    double getGrayRosseland( 
	const double targetTemperature, const double targetDensity );

    /*!
     * \breif Returns a vector of the opacity values for each energy
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
     *
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector has ngroups entries.  The number of groups is
     *         specified by the data file. 
     */
    std::vector<double> getMGRosseland( 
	const double targetTemperature, const double targetDensity );
 
  private:
    
    // IMPLEMENTATION

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

};



} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfOpacity.hh
//---------------------------------------------------------------------------//
