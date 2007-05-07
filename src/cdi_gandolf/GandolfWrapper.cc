//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfWrapper.cc
 * \author Kelly Thompson
 * \date   Thu Jul 13 15:31:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfWrapper.hh"

#ifdef rtt_cdi_gandolf_stub
#warning "USING GANDOLF STUB"
#endif

#include "ds++/Assert.hh"

#include <cstring>

namespace rtt_cdi_gandolf {
    
    namespace wrapper {
	
	using std::string;
	
	/*!
	 * \brief Converts a const string into a const char * that is padded
	 *        with white space.
	 *
	 * \param source The data in this string will be returned as a const
	 *               char * and padded with white space up to a length 
	 *               specified by n. 
	 * \param c1     A character string that has been allocated to length 
	 *               n by the calling routine.
	 * \param n      The length of c1 and the length of the returned
	 *               const char * c-string.
	 * \return A const char * of length n is returned.  It contains a
	 *         copy of source and is padded with white space.
	 */
	const char *s2ccwp( const string &source, char *c1, int n )
	    {
		// Create a string to hold the needed amount of padding.
		string padding(n-source.size(),' ');
		// copy the source string into a form that can modified.
		string s1(source);
		// append the requested amount of white space padding.
		s1.append(padding);
		// copy the string into the c-string.
		std::copy(s1.begin(),s1.end(),c1);
		return c1;
	    }
	
	//----------------------------------------//
	//                gmatids                 //
	//----------------------------------------//
	
	int wgmatids( const std::string& fname, vector<int>& matids, 
		      const int const_kmat, int &nmat ) 
	    {
#ifndef rtt_cdi_gandolf_stub
		// I could change this subroutine so that it identifies
		// nmat=kmat by repeatedly calling gmatids_().
		
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		// copy filename into a const char * array;
		char cfname[maxDataFilenameLength];
		const char * ccfname = s2ccwp( fname, cfname,
					       maxDataFilenameLength );
		
		// we don't know the value of nmat until runtime so we
		// must dynamically allocate a_matids.
		int *a_matids = new int [ const_kmat ];
		
		// remove const-ness.
		int kmat = const_kmat;

		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		int errorCode = 0;
		gmatids( ccfname, a_matids, kmat, nmat, errorCode );
		
		// ----------------------------------------
		// Copy the data back into C++ data types
		// ----------------------------------------
		
		// resize and update the vector matids fromt he array version.
		matids.resize( nmat );
		std::copy( a_matids, a_matids+nmat, matids.begin() );
		
		// Free up dynamic memory and return.
		delete [] a_matids;
		
		return errorCode;
		
#else // ifndef rtt_cdi_gandolf_stub
                return 1; // requested file not found.
#endif
	    } // end of gmatids
	
	
	//----------------------------------------//
	//                gkeys                   //
	//----------------------------------------//
	
	int wgkeys( const std::string &fname, const int &const_matid, 
		    vector<string> &vkeys, const int &const_kkeys, int &nkeys )
	    {
#ifndef rtt_cdi_gandolf_stub
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		// copy filename into a const char * array;
		char cfname[maxDataFilenameLength];
		const char * ccfname = s2ccwp( fname, cfname,
					       maxDataFilenameLength );
		
		// remove const-ness
		int matid = const_matid;
		int kkeys = const_kkeys;
		int ier = 0;
		
		// we do not know the value of numKeys until after we call 
		// gkeys() so we create the character array keys[][] to be 
		// maxKeys long.  This array will later be copied into the
		// vector vkeys that is returned to the calling program.
		
		// char keys[maxKeys][key_length];
		// char (*keys)[key_length] = new char[maxKeys][key_length];
		// delete [] keys;
		
		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		
		// This declaration doesn't guarantee that we have enough
		// memory for maxKeys * key_length characters.	    
		const char *bjkeys[maxKeys];
		
//  	    gkeys( ccfname, matid, *keys, kkeys, nkeys, ier );
		gkeys( ccfname, matid, bjkeys, kkeys, nkeys, ier );
		
		// 	    for ( int i=0; i < maxKeys; ++i )
		// 		std::cout << "bjkeys[" << i << "] = " << bjkeys[i] <<
		// 		    std::endl;
		
		
		// ----------------------------------------
		// Copy the data back into C++ data types
		// ----------------------------------------
		
		// Resize vkeys and copy the data from the char array 
		// into the vector of strings.
		vkeys.resize( nkeys );
		char key[key_length];
		for ( int i=0; i<nkeys; ++i )
		    {
			// copy all 24 characters of keys[i] into key.
			std::strncpy( key, bjkeys[i], key_length );
			// kill trailing whitespace.
			std::strtok( key, " " );
			// store the keyword in a vector.
			vkeys[i].assign( key, 0, std::strlen(key) );
		    }

		return ier;
		
#else //  rtt_cdi_gandolf_stub
                return 1; // requested file not found.
#endif
	    } // end of gkeys
	
	
	//----------------------------------------//
	//                gchgrids                //
	//----------------------------------------//
	
	int wgchgrids( const std::string &fname, const int &const_matid, 
		       int &nt, int &nrho, int &nhnu, int &ngray, 
		       int &nmg )
	    {
#ifndef rtt_cdi_gandolf_stub
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		// copy filename into a const char * array;
		char cfname[maxDataFilenameLength];
		const char * ccfname = s2ccwp( fname, cfname,
					       maxDataFilenameLength );
		
		// remove const-ness
		int matid = const_matid; // const
		int ier = 0;
		
		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		
		gchgrids( ccfname, matid, nt, nrho, nhnu,
			  ngray, nmg, ier );

		return ier;
		
#else
                return 1; // requested file not found.
#endif
	    } // end of gchgrids
	
	//----------------------------------------//
	//                ggetgray                //
	//----------------------------------------//
	
	int wggetgray( const string &fname,   
		       const int &const_matid, const string skey,
		       vector<double> &temps,  const int &const_kt,
		       vector<double> &rhos,   const int &const_krho,
		       vector<double> &data,   const int &const_kgray )
	    {
#ifndef rtt_cdi_gandolf_stub
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		// copy filename into a const char * array;
		char cfname[maxDataFilenameLength];
		const char * ccfname = s2ccwp( fname, cfname,
					       maxDataFilenameLength );
		
		// copy skey into a const char * array;
		char key[ key_length ];                           
		const char * cckey = s2ccwp( skey, key, key_length );
		
		// remove const-ness
		int matid = const_matid; // const
		int kt    = const_kt;    // const
		int krho  = const_krho;  // const
		int kgray = const_kgray; // const
		
		int ier = 0;
		
		// Allocate memory for double arrays (temps,rhos,data).
		// These will be copied into vector<double> objects later.
		double *array_temps = new double [kt];
		double *array_rhos  = new double [krho];
		double *array_data  = new double [kgray];

		// temporaries
		// since we have already loaded the grid size by
		// calling wgchgrids() our values for kXXX should be
		// identical to nXXX returned by ggetgray().
		int nt, nrho, ngray;
		
		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		
		ggetgray( ccfname,     matid, cckey,
			  array_temps, kt,    nt, 
			  array_rhos,  krho,  nrho,
			  array_data,  kgray, ngray,
			  ier );

		// We should have already loaded the grid size in
		// wgchgrids().  These numbers had better match!
// 		if ( nt    != kt     ||
// 		     nrho  != krho   ||
// 		     ngray != kgray ) 
// 		    return -20;

		// ----------------------------------------
		// Copy the data back into C++ data types
		// ----------------------------------------
		
		if ( ier == 0 ) // If ggetgray() returns an error
		    // return the error without filling the arrays.
		    {
			temps.resize(nt);
			rhos.resize(nrho);
			data.resize(ngray);
			
			std::copy( array_temps, array_temps+nt,   temps.begin() );
			std::copy( array_rhos,  array_rhos+nrho,  rhos.begin()  );
			std::copy( array_data,  array_data+ngray, data.begin()  );
		    }

		delete [] array_temps;
		delete [] array_rhos;
		delete [] array_data;
		
		return ier;

#else // #ifndef rtt_cdi_gandolf_stub
                return 1; // requested file not found.
#endif
	    } // end of ggetgray
	
	//----------------------------------------//
	//                gintgrlog               //
	//----------------------------------------//
	
	double wgintgrlog( const vector<double> &temps, const int &const_nt,
			   const vector<double> &rhos,  const int &const_nrho,
			   const vector<double> &data,  const int &const_ngray,
			   const double &const_tlog, const double &const_rlog ) 
	    {
#ifndef rtt_cdi_gandolf_stub
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		// remove const-ness;
		int nt    = const_nt;   
		int nrho  = const_nrho; 
		int ngray = const_ngray;
		double tlog = const_tlog; 
		double rlog = const_rlog;
		
		// Allocate memory for double arrays (temps,rhos,data).
		// We copy vector objects into these arrays before calling 
		// gintgrlog().
		double *array_temps = new double [nt];
		double *array_rhos  = new double [nrho];
		double *array_data  = new double [ngray];
		
		std::copy( temps.begin(), temps.end(), array_temps );
		std::copy( rhos.begin(),  rhos.end(),  array_rhos );
		std::copy( data.begin(),  data.end(),  array_data );
		
		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		
		// the solution
		double ans;

		gintgrlog( array_temps, nt, array_rhos, nrho,
			   array_data, ngray, tlog, rlog, ans );
		
		// no error code is returned from this function.
		// we don't need to copy any data back into C++ data
		// types.  The only return value is "ans" and it is
		// already in the correct format.
		
		delete [] array_temps;
		delete [] array_rhos;
		delete [] array_data;

		return ans;
		
#else // #ifndef rtt_cdi_gandolf_stub
                double ans = 0.0;
                return ans; // dummy scalar.
#endif
	    } // end of ginggrlog
	
	//----------------------------------------//
	//                ggetmg                  //
	//----------------------------------------//
	
	// Read data grid (temp,density,energy_bounds) and mg opacity
	// data.  Retrieve both the size of the data and the actual data.
	
	int wggetmg( const string &fname,  
		     const int &const_matid, const string skey,
		     vector<double> &temps,  const int &const_kt,
		     vector<double> &rhos,   const int &const_krho,
		     vector<double> &hnus,   const int &const_khnu,
		     vector<double> &data,   const int &const_kdata )
	    {
#ifndef rtt_cdi_gandolf_stub
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		// copy filename into a const char * array;
		char cfname[maxDataFilenameLength];
		const char * ccfname = s2ccwp( fname, cfname,
					       maxDataFilenameLength );
		
		// copy skey into a const char * array;
		char key[ key_length ];                           
		const char * cckey = s2ccwp( skey, key, key_length );
		
		// remove const-ness
		int matid = const_matid; 
		int kt    = const_kt;    
		int krho  = const_krho;  
		int khnu  = const_khnu;  
		int kdata = const_kdata; 
		
		// Allocate memory for double arrays (temps,rhos,data).
		// These will be copied into vector<double> objects later.
		double *array_temps = new double [kt];
		double *array_rhos  = new double [krho];
		double *array_hnus  = new double [khnu];
		double *array_data  = new double [kdata];

		// temporaries
		// since we have already loaded the grid size by
		// calling wgchgrids() our values for kXXX should be
		// identical to nXXX returned by ggetmg().
		int nt, nrho, nhnu, ndata, ier = 0;
		
		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		
		ggetmg( ccfname, matid, cckey,
			array_temps, kt,    nt, 
			array_rhos,  krho,  nrho,
			array_hnus,  khnu,  nhnu,
			array_data,  kdata, ndata,
			ier );
		
		// ----------------------------------------
		// Copy the data back into C++ data types
		// ----------------------------------------
		
		// If ggetmg() returns an error code then we don't
		// fill these vectors.  We simply return the error code.
		if ( ier == 0 )
		    {
			// resize data found in the Opacity object.
			temps.resize(nt);
			rhos.resize(nrho);
			hnus.resize(nhnu);
			data.resize(ndata);
			
			std::copy( array_temps, array_temps+nt,   temps.begin() );
			std::copy( array_rhos,  array_rhos+nrho,  rhos.begin()  );
			std::copy( array_hnus,  array_hnus+nhnu,  hnus.begin()  );
			std::copy( array_data,  array_data+ndata, data.begin()  );
			
		    }
		
		// free up dynamically allocated memory
		
		delete [] array_temps;
		delete [] array_rhos;
		delete [] array_hnus;
		delete [] array_data;

		return ier;
		
#else // #ifndef rtt_cdi_gandolf_stub
                return 1; // requested file not found.
#endif
	    } // end of ggetgray
	
	//----------------------------------------//
	//                gintmglog               //
	//----------------------------------------//
	
	vector<double>
	    wgintmglog( const vector<double> &temps, const int &const_nt,
			const vector<double> &rhos,  const int &const_nrho,
			const int &const_nhnu,
			const vector<double> &data,  const int &const_ndata,
			const double &const_tlog, const double &const_rlog )
	    {
#ifndef rtt_cdi_gandolf_stub
		// ----------------------------------------
		// Create simple flat data types
		// ----------------------------------------
		
		const int ngroups = const_nhnu-1;
		
		// Remove const-ness.
		int nt    = const_nt; 
		int nrho  = const_nrho;
		int nhnu  = const_nhnu;
		int ndata = const_ndata;
		double tlog = const_tlog;
		double rlog = const_rlog;
		
		// Allocate memory for double arrays (temps,rhos,data).
		// We copy vector objects into these arrays before calling 
		// gintgrlog().
		double *array_temps = new double [nt];
		double *array_rhos  = new double [nrho];
		double *array_data  = new double [ndata];
		
		std::copy( temps.begin(), temps.end(), array_temps );
		std::copy( rhos.begin(),  rhos.end(),  array_rhos );
		std::copy( data.begin(),  data.end(),  array_data );
		
		// Allocate apace for the solution.
		double *array_ansmg = new double [ngroups];
		
		// --------------------------------------------------
		// call the Gandolf library function
		// --------------------------------------------------
		
		gintmglog( array_temps, nt, array_rhos, nrho,
			   nhnu, array_data, ndata, tlog, rlog,
			   array_ansmg ); 
		
		// ----------------------------------------
		// Copy the data back into C++ data types
		// ----------------------------------------
		
		// Create a vector<double> container for the solution;
		vector<double> ansmg( ngroups );
		// Copy the Multigroup Opacity into the new container.
		std::copy( array_ansmg, array_ansmg+ngroups, 
			   ansmg.begin() );
		
		// release space required by temps;
		delete [] array_temps;
		delete [] array_rhos;
		delete [] array_data;
		delete [] array_ansmg;

		return ansmg;
		
#else // rtt_cdi_gandolf_stub
        vector<double> ansmg( const_nhnu - 1, 0.0 );
        return ansmg; // dummy vector.
#endif
	    } // end of wgintmglog
	
    } // end namespace wrapper
    
} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
//                         end of cdi/GandolfWrapper.cc
//---------------------------------------------------------------------------//
