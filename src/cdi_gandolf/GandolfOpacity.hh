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

#include "cdi/Opacity.hh"
#include "GandolfWrapper.hh"

namespace rtt_cdi_gandolf
{

using std::string;
using std::vector;

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
     * \brief IPCRESS data filename
     */
    const string dataFilename;

    /*!
     * \brief List of material identifiers found in the data file.
     */
    vector<int> matIDs;

    /*!
     * \brief The material ID that this instance of GandolfOpacity 
     *        represents.
     */
    int matID;

    /*!
     * \brief List of keywords for current material found in the data
     *        file. 
     */
    char keys[maxKeys][key_length];

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
     * \brief The temperature grid found in the data file.
     */
    vector<double> temperatures;

    /*!
     * \brief The density grid found in the data file.
     */
    vector<double> densities;

    /*!
     * \brief The energy group boundaries found in the data file.
     */
    vector<double> groupBoundaries;

    /*!
     * \brief The gray opacity table found in the data file.
     */
    vector<double> grayOpacities;

    /*!
     * \brief The multigroup opacity table found in the data file.
     */
    vector<double> MGOpacities;

  public:

    // CREATORS
    
    /*!
     * \brief GandolfOpacity constructor only requiring data filename.
     *
     * \sa This GandolfOpacity constructor only requires the data
     *     filename.   This constructor will provide the client
     *     with an object that allows them to see a list of
     *     available material ID tags.  This constructor should be
     *     used very rarely.  The primary constructor for this
     *     class is GandolfOpacity( string, int );
     *
     * \param _data_filename The name of the IPCRESS data file.  While 
     *                       this class can use any filename the F77
     *                       Gandolf library requires the name to have
     *                       80 characters or less. 
     */
    GandolfOpacity( string _data_filename );

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
    GandolfOpacity( string _data_filename, int _matid );

    // (defaulted) GandolfOpacity(const GandolfOpacity &rhs);
    // (defaulte) ~GandolfOpacity()

    // MANIPULATORS
    
    //defaulted GandolfOpacity& operator=(const GandolfOpacity &rhs);

    // ACCESSORS

    /*!
     * \brief Returns the IPCRESS data filename.
     */
    string getDataFilename() { return dataFilename; }

    /*!
     * \brief Returns a list of material identifiers available in the
     *        IPCRESS file.
     */
    vector<int> getMatIDs() 
    { 
	return matIDs; 
    }

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
    double getGray( const double targetTemperature,
		    const double targetDensity );

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
    vector<double> getMG( const double targetTemperature, 
	 		  const double targetDensity );
 
  private:
    
    // IMPLEMENTATION
};

/*! 
 * \brief This function returns "true" if "key" is found in the list
 *        of "keys".
 */
 bool key_available( char key[], char keys[][key_length], int numKeys ); 

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfOpacity.hh
//---------------------------------------------------------------------------//
