//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfFile.cc
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Implementation file for GandolfFile class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfFile.hh"
#include "GandolfWrapper.hh"
#include "GandolfException.hh"

namespace rtt_cdi_gandolf
{
    /*!
     * \brief The standard GandolfFile constructor.
     */
    GandolfFile::GandolfFile( const std::string& gandolfDataFilename )
	: dataFilename( gandolfDataFilename )
	{
	    // Gandolf will only look at the first "maxDataFilenameLength"
	    // characters of the data filename.  We must require that the
	    // given name is less than "maxDataFilenameLength" characters.
	    if ( dataFilename.length() >
		 wrapper::maxDataFilenameLength )
		throw gmatidsException( -1 );
	 
	    // This call to Gandolf validates the datafile and if
	    // successful returns a list of material identifiers for which 
	    // the materials that exist in the data file.
	    int errorCode = 0;
	    wrapper::wgmatids( dataFilename, matIDs, wrapper::maxMaterials,
			       numMaterials, errorCode ); 
	    
	    if ( errorCode != 0 )
		throw gmatidsException( errorCode );

	}

    /*!
     * \brief Indicate if the requested material id is available in
     *        the data file.
     */
    bool GandolfFile::materialFound( int matid ) const
	{
	    // Loop over all available materials.  If the requested
	    // material id matches on in the list then return true.
	    // If we reach the end of the list without a match return
	    // false. 
	    for ( int i=0; i<numMaterials; ++i )
		if ( matid == matIDs[i] ) return true;
	    return false;
	    
	} // end of materialFound()

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
// end of GandolfFile.cc
//---------------------------------------------------------------------------//
