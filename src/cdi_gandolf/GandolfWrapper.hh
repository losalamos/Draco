//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfWrapper.hh
 * \author Kelly Thompson
 * \date   Thu Jul 13 15:31:56 2000
 * \brief  Header file for GandolfWrapper
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfWrapper_hh__
#define __cdi_gandolf_GandolfWrapper_hh__

#include "cdi_gandolf/config.h"

#include <string>
#include <cstring>
#include <vector>

namespace rtt_cdi_gandolf {

    namespace wrapper {

	using std::string;
	using std::vector;
 
	//====================================================================
	/*! 
	 * \brief Accessing the Gandolf library routines.
	 *
	 * The Gandolf routines are written in FORTRAN and provide a
	 * primative interface for extrating data.  These routines
	 * require maximum data sizes to be specifed for memory
	 * allocation in the FORTRAN library.  These maximums are are
	 * created here and should be available any where in the
	 * rtt_cdi_gandolf namespace.  
	 */
	//====================================================================

	// maxDataFilenameLength and key_length are size descriptors
	// expected by Gandolf.  There is no interactive way of asking
	// Gandolf for these values.  I have included these two
	// variables here so that if Gandolf ever changes its
	// definitions we can accomodate those changes easily.
	
	/*!
	 * \brief The maximum length of the data filename.  This
	 * 	  length is set by the Gandolf libraries.
	 */
	const size_t maxDataFilenameLength = 80;

	/*!
	 * \brief The length of each descriptor key (set by Gandolf).
	 */
	const size_t key_length = 24;  
	
	// maxMaterials and maxKeys are size descriptors that define
	// the pre-allocated sizes of the matIDs and keys arrays.
	// Ideally we would set both of these variables to the actual
	// number of materials in the data file and the actual number
	// of keys per material in the file.  Unfortunately, Gandolf
	// does not provide a mechanism for gaining this information
	// interactively.
	//
	// The only alternative would be to set numMaterials and
	// numKeys to a dummy value (e.g. 1) and then call the Gandolf
	// function gmatids().  If the number of materials found in
	// the IPCRESS file exceeds this size matids() will return an
	// error code but it will not tell you how many materials were
	// found.  We could create a loop that "seeks" the correct
	// number of materials by repeatedly calling gmatids() until
	// numMaterials is equal the actual number of materials in the
	// dataFile.

	/*!
	 * \brief Maximum number of materials allowed in the IPCRESS file.
	 */
        const size_t maxMaterials = 128;  // other codes use this magic value.
	
	/*!
	 * \brief Maximum number of data keys per material.
	 */
	const size_t maxKeys = 26;
	
	//====================================================================
	/*! 
	 * \brief C++ Gandolf wrapper routines.
	 *
	 * The Gandolf routines are written in FORTRAN.  The following
	 * are C++ prototypes that mimic the F77 Gandolf functions.
	 * Each of these routines flattens the data types and then
	 * calls the Gandolf library's F77 functions.  
	 */
	//====================================================================
	
	/*!
	 * \brief Retrieve a list of material identifiers assocaited with the
	 *        specified data file. 
	 *
	 * \param fname  The name of the IPCRESS data file.
	 * \param matids A list of material identifiers associated with fname.
	 * \param kamt   The maximum number of materials for which adequate
	 *               memory has been allocated.
	 * \param nmat   Actual number of materials found in the IPCRESS data
	 *               file.
	 * \param ier    Returned error code.  A value of zero indicates
	 *               sucess.
	 * \return       ier (also matids and nmat).
	 */
	int wgmatids( const string& fname, vector<int>& matids, 
		      const int kmat, int &nmat );
	/*!
	 * \brief Retrieve a list of keys that specify the types of
	 *        information available for the spacified material in the
	 *        IPCRESS data file. 
	 *
	 * \param fname  The name of the IPCRESS data file.
	 * \param matid  The material identifier for the material we are
	 *               querying. 
	 * \param keys   A list of character string identifiers.  These
	 *               identifiers specify what information is available for
	 *               the specified material. 
	 * \param kkeys  The maximum number of keys for which adequate memory
	 *               has been allocated.
	 * \param nkeys  Actual number of keys found for the material.
	 * \param ier    Returned error code.  A value of zero indicates
	 *               sucess.
	 * \return       ier (also returns keys and nkeys).
	 */
	int wgkeys( const string &fname, const int &matid, 
		    vector<string> &vkeys, const int &kkeys, size_t &nkeys );
	
	/*!
	 * \brief Retrieves the size of the data grid including the number of
	 *        temperature, density, energy group boundary, gray opacity
	 *        and multigroup opacity data points.
	 *
	 * \param fname  The name of the IPCRESS data file.
	 * \param matid  The material identifier for the material we are
	 *               querying. 
	 * \param nt     The number of temperature bins used in the data grid.
	 * \param nrho   The number of density bins used in the data grid.
	 * \param nhnu   The number of energy group boundaries.
	 * \param ngray  The number of gray opacity data points.
	 * \param nmg    The number of multigroup opacity data points.
	 * \param ier    Returned error code.  A value of zero indicates
	 *               sucess.
	 * \return       ier (also returns nt, nrho, nhnu, ngray, and nmg )
	 */
	int wgchgrids( const string &fname, const int &matid, int &nt,
		       int &nrho, int &nhnu, int &ngray, int &nmg );
 
	/*!
	 * \brief Retrieves the gray opacity data grid including the
	 *        temperature and density bin values.
	 *
	 * \param fname  The name of the IPCRESS data file.
	 * \param matid  The material identifier for the material we are
	 *               querying. 
	 * \param key    A character string identifier that specifies the type 
	 *               of data to extract from the data file (e.g. "rgray" = 
	 *               Rosseland gray opacities).
	 * \param temps  The temperature grid in keV.
	 * \param kt     The maximum number of temperatures for which adequate 
	 *               memory has been allocated.
	 * \param rhos   The density grid in g/cm^3.
	 * \param krho   The maximum number of dnesities for which adequate
	 *               memory has been allocated.
	 * \param gray   The gray opacity values in cm^2/g.
	 * \param kgray  The maximum number of opacities for which adequate
	 *               memory has been allocated.
	 * \param ier    Returned error code.  A value of zero indicates
	 *               sucess.
	 * \return       ier (also returns temps, rhos, and gray)
	 */
	int wggetgray( const string &fname,   
		       const int &matid,      const string key, 
		       vector<double> &temps, const int &kt,
		       vector<double> &rhos,  const int &krho,
		       vector<double> &gray,  const int &kgray );	

	/*!
	 * \brief Returns a gray opacity value based on user specified values
	 *        for temperature and density.  This routine interpolates from 
	 *        the data from ggetgray() to find the desired opacity.
	 *
	 * \param temps  The log of the temperature grid in keV.
	 * \param nt     The number of temperature bins used in the data grid.
	 * \param rhos   The log of the density grid in g/cm^3.
	 * \param nrho   The number of density bins used in the data grid.
	 * \param gray   The log of the gray opacity values in cm^2/g.
	 * \param ngray  The number of gray opacity data points.
	 * \param tlog   The log of the desired temperature.
	 * \param rlog   The log of the desired density.
	 * \return       The interpolated opacity.
	 */
	double wgintgrlog( const vector<double> &temps, const int &nt,
			   const vector<double> &rhos,  const int &nrho,
			   const vector<double> &gray,  const int &ngray,
			   const double &tlog, const double &rlog);
	
	/*!
	 * \brief Retrieves the multigroup opacity data grid including the
	 *        temperature bin, density bin and energy boundary values.
	 *
	 * \param fname  The name of the IPCRESS data file.
	 * \param matid  The material identifier for the material we are
	 *               querying. 
	 * \param key    A character string identifier that specifies the type 
	 *               of data to extract from the data file (e.g. "rgray" = 
	 *               Rosseland gray opacities).
	 * \param temps  The temperature grid in keV.
	 * \param kt     The maximum number of temperatures for which adequate 
	 *               memory has been allocated.
	 * \param rhos   The density grid in g/cm^3.
	 * \param krho   The maximum number of densities for which adequate
	 *               memory has been allocated.
	 * \param hnus   The energy group boundaries in keV.
	 * \param khnu   The maximum number of energy group boundaries for
	 *               which adequate memory has been allocated.
	 * \param data   The multigroup opacity values in cm^2/g.
	 * \param kdata  The maximum number of opacities for which adequate
	 *               memory has been allocated.
	 * \param ier    Returned error code.  A value of zero indicates
	 *               sucess.
	 * \return       ier ( also returns temps, rhos, hnus, and data.)
	 */
	int wggetmg( const string &fname,   
		     const int &matid,      const string key, 
		     vector<double> &temps, const int &kt,
		     vector<double> &rhos,  const int &krho,
		     vector<double> &hnus,  const int &khnu,
		     vector<double> &data,  const int &kdata );
	
	/*!
	 * \brief Returns a vector of multigroup opacity values based on user
	 *        specified values for temperature and density.  This routine
	 *        interpolates from the data obtained from ggetmg() to find
	 *        the desired opacities. 
	 *
	 * \param temps  The log of the temperature grid in keV.
	 * \param nt     The number of temperature bins used in the data grid.
	 * \param rhos   The log of the density grid in g/cm^3.
	 * \param nrho   The number of density bins used in the data grid.
	 * \param nhnu   The number of energy group boundaries in the data grid.
	 * \param data   The log of the multigroup opacity values in cm^2/g.
	 * \param ndata  The number of multigroup opacity data points.
	 * \param tlog   The log of the desired temperature.
	 * \param rlog   The log of the desired density.
	 * \return       A vector of interpolated multigroup opacity values.
	 *
	 */
	vector<double> 
	    wgintmglog( const vector<double> &temps, const int &nt,
			const vector<double> &rhos,  const int &nrho,
			const int &nhnu,
			const vector<double> &data,  const int &ndata,
			const double &tlog, const double &rlog );
	
    } // end namepsace wrapper
} // end namespace rtt_cdi_gandolf

#ifndef rtt_cdi_gandolf_stub

// Function prototypes for Gandolf F77 subroutines.
//---------------------------------------------------------------------------//

extern "C" {
    
    void gmatids( const char *cfname, int *matids, int &ckmat,
		  int &nmat, int &ier );
    
    //     void gkeys( const char* cfname, int &matid, 
    // 		const char* keys,
    // 		int &kkeys, int &nkeys, int &ier );
    
    void gkeys( const char* cfname, int &matid, 
		const char* bjkeys[],
		int &kkeys, int &nkeys, int &ier );
    
    //     void gkeys( const char* cfname, int &matid, 
    // 		char (*keys)[rtt_cdi_gandolf::wrapper::key_length],
    // 		int &kkeys, int &nkeys, int &ier );
    
    void gchgrids( const char *cfname, int &matid, int &nt,
		   int &nrho, int &nhnu, 
		   int &ngray, int &nmg, int &ier );
    
    void ggetgray( const char *cfname,  int &matid, const char *key, 
		   double *temps, int &kt,    int &nt, 
		   double *rhos,  int &krho,  int &nrho,
		   double *gray,  int &kgray, int &ngray,
		   int &ier );
    
    void gintgrlog( double *temps, int &nt,
		    double *rhos,  int &nrho,
		    double *gray,  int &ngray,
		    double &tlog, double &rlog, double &ans );
    
    void ggetmg( const char *cfname,   int &matid, const char *key, 
		 double *temps,  int &kt,    int &nt,
		 double *rhos,   int &krho,  int &nrho,
		 double *hnus,   int &khnu,  int &nhnu,
		 double *data,   int &kdata, int &ndata,
		 int &ier );

    void gintmglog( double *temps, int &nt,
		    double *rhos,  int &nrho,
		    int &nhnu,
		    double *data,  int &ndata,
		    double &tlog,  double &rlog, 
		    double *ansmg );

} // end of extern "C" block

#endif // rtt_cdi_gandolf_stub

#endif // __cdi_gandolf_GandolfWrapper_hh__

//---------------------------------------------------------------------------//
//                     end of cdi/GandolfWrapper.hh
//---------------------------------------------------------------------------//
