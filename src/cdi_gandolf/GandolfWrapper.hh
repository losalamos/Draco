//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfWrapper.hh
 * \author Kelly Thompson
 * \date   Thu Jul 13 15:31:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfWrapper_hh__
#define __cdi_gandolf_GandolfWrapper_hh__

#include <string>
#include <vector>

#include "GandolfOpacity.hh"

// function prototypes
//---------------------------------------------------------------------------//

namespace rtt_cdi_gandolf {

    using std::string;
    using std::vector;
 
    void gmatids( const string &fname , int matids[], 
		  const int kmat,int &nmat, int &ier );

    void gkeys( const string &fname, const int &matid, 
		char keys[][key_length],
		const int kkeys, int &nkeys, int &ier );

    void gchgrids( const string &fname, const int &matid,
		   int &nt, int &nrho, int &nhnu, int &ngray, int &nmg,
		   int &ier );
 
    void ggetgray( const string &fname,   const int &matid, char *key, 
		   vector<double> &temps, const int &kt,    int &nt, 
		   vector<double> &rhos,  const int &krho,  int &nrho,
		   vector<double> &gray,  const int &kgray, int &ngray,
		   int &ier );

    void gintgrlog( const vector<double> &temps, const int &nt,
		    const vector<double> &rhos,  const int &nrho,
		    const vector<double> &gray,  const int &ngray,
		    const double &tlog, const double &rlog, double &ans );

    void ggetmg( const string &fname,   const int &matid, char *key, 
		 vector<double> &temps, const int &kt,    int &nt,
		 vector<double> &rhos,  const int &krho,  int &nrho,
		 vector<double> &hnus,  const int &khnu,  int &nhnu,
		 vector<double> &data,  const int &kdata, int &ndata,
		 int &ier );
    
    void gintmglog( const vector<double> &temps, const int &nt,
		    const vector<double> &rhos,  const int &nrho,
		    const int &nhnu,
		    const vector<double> &data,  const int &ndata,
		    const double &tlog, const double &rlog, 
		    vector<double> &ansmg );

// I may need to create a "long int" version of the calling routines.
//     void gmatids( const std::string &fname , long int matids[], 
// 		  const long int kmat, long int &nmat, long int &ier );

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
    void extc_gkeys( char *cfname, long int &matid, char keys[][24],
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
