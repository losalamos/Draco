//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfOpacity.hh"
#include "GandolfWrapper.hh"

#include <string>

namespace rtt_cdi_gandolf
{
/*!
 * \brief Constructor for Gandolf Opacity reader (as a part of CDI).
 */
GandolfOpacity::GandolfOpacity( string _data_filename )
    : dataFilename ( _data_filename )
	  
    {
	cout << "In GandolfOpacity::GandolfOpacity()" << endl;
	
	// Assert that the opacity data file exists and that it can be read.
	// Let Gandolf do this for now.
	//
	// ofstream infile( gandolfFilename );
	// Insist(!infile,"Could not open Gandalf data file for reading.");
	
	const int kmat = 2;
	int amatids[kmat], nmat=0, ier=0;

	// initialize amatids to zero.
	for ( int i=0; i<kmat; ++i ) {
	    // amatids[kmat]=0;
	    cout << "GandolfOpacity::GandolfOpacity() amatids[" 
		 << i << "] = " << amatids[i] << endl;
	}

	gmatids( dataFilename, amatids, kmat, nmat, ier );

// 	cout << "GandolfOpacity::GandolfOpacity()"
// 	     << "  back from call to gmatids" << endl
// 	     << "  we found ier = " << ier << endl
// 	     << "           nmat = " << nmat << endl
// 	     << "           amatids[0] = " << amatids[0] << endl 
// 	     << endl;

	// copy amatids into the vector matIDs
	matIDs.resize(nmat);

	// use sdnolen's example to copy amatids into matIDs (STL
	// mechanism).
	for ( int i=0; i<nmat; ++i ) {
	    matIDs[i] = amatids[i];
	    cout << "GandolfOpacity::GandolfOpacity() matIDs[" 
		 << i << "] = " << matIDs[i] << endl;
	}
	cout << endl;
    } // end GandolfOpacity::GandolfOpacity()


/*!
 * \brief Return a Rosseland Mean Gray Opacity value for the user
 *        specified temperature and density.
 */
 double GandolfOpacity::getGray( const double temp, const double density )
     {
	 double grayOpacity = 0.0;
	 
	 return grayOpacity;
     }

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
//                              end of GandolfOpacity.cc
//---------------------------------------------------------------------------//
