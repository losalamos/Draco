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

#include "ds++/Assert.hh"

namespace rtt_cdi_gandolf
{
 
//===========================================================================//
/*!
 * \class GandolfFile
 *
 */
// revision history:
// -----------------
// 0) original
// 
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
     * \brief Number of materials found in the data file.
     */
    int numMaterials;

    /*!
     * \brief A list of material IDs found in the data file.
     */
    std::vector<int> matIDs;

  public:

    // CREATORS
    
    GandolfFile( const std::string& gandolfDataFilename );
    // (defaulted) GandolfFile(const GandolfFile &rhs);
    // (defaulted) ~GandolfFile();

    // MANIPULATORS
    
    // (defaulted) GandolfFile& operator=(const GandolfFile &rhs);

    // ACCESSORS

    const std::string& getDataFilename() const 
    { 
	return dataFilename;
    }
    
    int getNumMaterials() const
    {
	return numMaterials;
    }

    const std::vector<int>& getMatIDs() const
    {
	return matIDs;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi_gandolf

#endif                          // __cdi_gandolf_GandolfFile_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_gandolf/GandolfFile.hh
//---------------------------------------------------------------------------//
