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

#include <string>
#include <cstring>
#include <vector>

namespace rtt_cdi_gandolf {

using std::string;
using std::vector;
 
//===========================================================================//
/*! 
 * \brief Accessing the Gandolf library routines.
 *
 * \sa The Gandolf routines are written in FORTRAN and provide a
 * primative interface for extrating data.  These routines require
 * maximum data sizes to be specifed for memory allocation in the
 * FORTRAN library.  These maximums are are created here and should be 
 * available any where in the rtt_cdi_gandolf namespace.
 */
//===========================================================================//

/*!
 * \brief The maximum length of the data filename.  This length is set 
 *        by the Gandolf libraries.
 */
 const int maxDataFilenameLength = 80;

 /*!
  * \brief The length of each descriptor key (set by Gandolf).
  */
 const int key_length = 24;  

 /*!
  * \brief Maximum number of materials allowed in the IPCRESS file.
  */
 const int maxMaterials = 10; 

 /*!
  * \brief Maximum number of data keys per material.
  */
 const int maxKeys = 25;

 /*!
  * \brief The maximum number of temperatures in the data grid.
  */
 const int maxTemps = 10;

 /*!
  * \brief The maximum number of densities in the data grid.
  */
 const int maxDensities = 10;

 /*!
  * \brief The maximum number of group energy boundaries in the data
  *        grid.  The maximum number of energy groups is 
  *        ( maxGroupBoundaries - 1).
  */
 const int maxGroupBoundaries = 35;

 /*!
  * \brief The maximum number of gray opacities in the data file.
  */
 const int maxGrayOpacities = maxTemps * maxDensities;

 /*!
  * \brief The maximum number of multigroup opacities in the data
  *        file.
 */
 const int maxMGOpacities = 
     maxGrayOpacities * ( maxGroupBoundaries - 1 );

//===========================================================================//
/*! 
 * \brief C++ Gandolf wrapper routines.
 *
 * \sa The Gandolf routines are written in FORTRAN.  The following are 
 * C++ prototypes that mimic the F77 Gandolf functions.  Each of these 
 * routines flattens the data types and then calls the Gandolf
 * library's F77 functions.
 */
//===========================================================================//

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
 * \return       matids, nmat and ier.
 *
 */
    void gmatids( const string &fname , vector<int> &matids, 
		  const int kmat, int &nmat, int &ier );
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
 * \return       keys, nkeys, ier.
 *
 */
    void gkeys( const string &fname, const int &matid, 
		vector<string> &vkeys,
		const int kkeys, int &nkeys, int &ier );

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
 * \return       nt, nrho, nhnu, ngray, nmg.
 *
 */
    void gchgrids( const string &fname, const int &matid,
		   int &nt, int &nrho, int &nhnu, int &ngray, int &nmg,
		   int &ier );
 
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
 * \param nt     The number of temperature bins used in the data grid.
 * \param rhos   The density grid in g/cm^3.
 * \param krho   The maximum number of dnesities for which adequate
 *               memory has been allocated.
 * \param nrho   The number of density bins used in the data grid.
 * \param gray   The gray opacity values in cm^2/g.
 * \param kgray  The maximum number of opacities for which adequate
 *               memory has been allocated.
 * \param ngray  The number of gray opacity data points.
 * \param ier    Returned error code.  A value of zero indicates
 *               sucess.
 * \return       temps, nt, rhos, nrho, gray, ngray.
 *
 */
    void ggetgray( const string &fname,   const int &matid, const string key, 
		   vector<double> &temps, const int &kt,    int &nt, 
		   vector<double> &rhos,  const int &krho,  int &nrho,
		   vector<double> &gray,  const int &kgray, int &ngray,
		   int &ier );

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
 * \param ans    The interpolated gray opacity value.
 * \return       ans.
 *
 */
    void gintgrlog( const vector<double> &temps, const int &nt,
		    const vector<double> &rhos,  const int &nrho,
		    const vector<double> &gray,  const int &ngray,
		    const double &tlog, const double &rlog, double &ans );

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
 * \param nt     The number of temperature bins used in the data grid.
 * \param rhos   The density grid in g/cm^3.
 * \param krho   The maximum number of densities for which adequate
 *               memory has been allocated.
 * \param nrho   The number of density bins used in the data grid.
 * \param hnus   The energy group boundaries in keV.
 * \param khnu   The maximum number of energy group boundaries for
 *               which adequate memory has been allocated.
 * \param nhnu   The actual number of energy group boundaries found in 
 *               the data file.
 * \param data   The multigroup opacity values in cm^2/g.
 * \param kdata  The maximum number of opacities for which adequate
 *               memory has been allocated.
 * \param ndata  The number of multigroup opacity data points found in 
 *               the data file.
 * \param ier    Returned error code.  A value of zero indicates
 *               sucess.
 * \return       temps, nt, rhos, nrho, hnus, nhnu, data, ndata, ier.
 *
 */
    void ggetmg( const string &fname,   const int &matid, const string key, 
		 vector<double> &temps, const int &kt,    int &nt,
		 vector<double> &rhos,  const int &krho,  int &nrho,
		 vector<double> &hnus,  const int &khnu,  int &nhnu,
		 vector<double> &data,  const int &kdata, int &ndata,
		 int &ier );
    
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
 * \param ansmg  The interpolated multigroup opacity values.
 * \return       ansmg.
 *
 */
    void gintmglog( const vector<double> &temps, const int &nt,
		    const vector<double> &rhos,  const int &nrho,
		    const int &nhnu,
		    const vector<double> &data,  const int &ndata,
		    const double &tlog, const double &rlog, 
		    vector<double> &ansmg );

/*!
 * \brief copies the source string into the target c-string.
 *
 * \param source       The string that needs to be copied into a c-string.
 * \param target       The c-string that is being filled.
 * \param targetLength The length of of the c-string "target[]".
 * \return             target[]
 *
 */
    void string2char( const string &source, char target[], 
		      int targetLength );	

} // end namespace rtt_cdi_gandolf



// Handle machine specific FORTRAN name linkage.
//---------------------------------------------------------------------------//

#if defined(sun) || defined(__sun) || defined(__sgi) || defined(__linux)
    
#define extc_gmatids gmatids_
#define extc_gkeys gkeys_    
#define extc_gchgrids gchgrids_
#define extc_ggetgray ggetgray_
#define extc_gintgrlog gintgrlog_
#define extc_ggetmg ggetmg_
#define extc_gintmglog gintmglog_

// Add defines for gandolf_integer and gandolf_double?

#endif
    


// Function prototypes for Gandolf F77 subroutines.
//---------------------------------------------------------------------------//

// The Gandolf library was compiled with -i8 so we must use "long int" 
// values.

extern "C" {

    void extc_gmatids( char *cfname, long int *matids, long int &ckmat,
		       long int &nmat, long int &ier );

    // key_length is specified to be 24 by the Gandolf standard.  This 
    // variable is set in the rtt_cdi_gandolf namespace but since this 
    // "extern C" block is outside of that namespace we must specify
    // this length manually.
    void extc_gkeys( char *cfname, long int &matid, 
		     char keys[][rtt_cdi_gandolf::key_length],
		     long int &kkeys, long int &nkeys, long int &ier );

    void extc_gchgrids( char *cfname, long int &matid, long int &nt,
			long int &nrho, long int &nhnu, 
			long int &ngray, long int &nmg, long int &ier );

    void extc_ggetgray( char *cfname,  long int &matid, char *key, 
			double *temps, long int &kt,    long int &nt, 
			double *rhos,  long int &krho,  long int &nrho,
			double *gray,  long int &kgray, long int &ngray,
			long int &ier );

    void extc_gintgrlog( double *temps, long int &nt,
			 double *rhos,  long int &nrho,
			 double *gray,  long int &ngray,
			 double &tlog, double &rlog, double &ans );

    void extc_ggetmg( char *cfname,   long int &matid, char *key, 
		      double *temps,  long int &kt,    long int &nt,
		      double *rhos,   long int &krho,  long int &nrho,
		      double *hnus,   long int &khnu,  long int &nhnu,
		      double *data,   long int &kdata, long int &ndata,
		      long int &ier );

    void extc_gintmglog( double *temps, long int &nt,
		         double *rhos,  long int &nrho,
			 long int &nhnu,
			 double *data,  long int &ndata,
			 double &tlog,  double &rlog, 
			 double *ansmg );

} // end of extern "C" block

#endif // __cdi_gandolf_GandolfWrapper_hh__

//---------------------------------------------------------------------------//
//                     end of cdi/GandolfWrapper.hh
//---------------------------------------------------------------------------//
