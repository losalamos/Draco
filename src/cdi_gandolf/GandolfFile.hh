//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfFile.hh
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Header file for GandolfFile class
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfFile_hh__
#define __cdi_gandolf_GandolfFile_hh__

#include <string>
#include <vector>

namespace rtt_cdi_gandolf
{
 
//===========================================================================//
/*!
 * \class GandolfFile
 *
 * \brief This class controls access the the physical IPCRESS data
 *        file for GandolfOpacity.  Only one GandolfFile object should 
 *        exist for each data file.  Several GandolfOpacity objects
 *        will access the same GandolfFile object (one per material
 *        found in the data file).
 *
 * This class is designed to be used in conjuction with
 * GandolfOpacity.  The client code should create a GandolfFile object 
 * and that object is passed to the GandolfOpacity constructor to
 * create a link between the opacity object and the IPCRESS data file.
 */

/*!
 * \example cdi_gandolf/test/tGandolfFile.cc
 *
 * Example of GandolfFile useage independent of GandolfOpacity or CDI.
 *
 */
//===========================================================================//

class GandolfFile 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
    /*!
     * \brief IPCRESS data filename
     */
    const std::string dataFilename;

    /*!
     * \brief Number of materials found in the data file.  This is not 
     *     a const value because it will be set after accessing the
     *     the data file.
     */
    size_t numMaterials;

    /*!
     * \brief A list of material IDs found in the data file.
     */
    std::vector<int> matIDs;

  public:

    // CREATORS
    
    /*!
     * \brief Standard GandolfFile constructor.
     *
     *    This is the standard GandolfFile constructor.  This object
     *    is typically instantiated as a smart pointer.
     *
     * \param gandolfDataFilename A string that contains the name of
     *     the Gandolf data file in IPCRESS format.  The f77 Gandolf
     *     vendor library expects a name with 80 characters or less.
     *     If the filename is longer than 80 characters the library
     *     will not be able to open the file.
     */
    GandolfFile( const std::string& gandolfDataFilename );

    // (defaulted) GandolfFile(const GandolfFile &rhs);
    // (defaulted) ~GandolfFile();

    // MANIPULATORS
    
    // (defaulted) GandolfFile& operator=(const GandolfFile &rhs);

    // ACCESSORS

    /*!
     * \brief Returns the IPCRESS data filename.
     */
    std::string const & getDataFilename() const 
    { 
	return dataFilename;
    }
    
    /*!
     * \brief Returns the number of materials found in the data file.
     */
    size_t getNumMaterials() const
    {
	return numMaterials;
    }

    /*!
     * \brief Returns a list of material identifiers found in the data 
     *     file.
     */
    const std::vector<int>& getMatIDs() const
    {
	return matIDs;
    }

    /*!
     * \brief Indicate if the requested material id is available in
     *        the data file.
     */
    bool materialFound( int matid ) const;

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfFile_hh__

//---------------------------------------------------------------------------//
// end of cdi_gandolf/GandolfFile.hh
//---------------------------------------------------------------------------//
