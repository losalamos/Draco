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

#include "ds++/Assert.hh"

#include <iostream>
#include <cstring>

namespace rtt_cdi_gandolf {

namespace wrapper {

using std::string;

 /*!
  * \brief Converts a const sring into a const char * that is padded
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
    
    void gmatids( const std::string &fname , vector<int> &matids, 
		  const int const_kmat, int &nmat, int &ier ) 
	{

	    // I could change this subroutine so that it identifies
	    // nmat=kmat by repeatedly calling gmatids_().

	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------

	    // copy filename into a const char * array;
    	    char cfname[maxDataFilenameLength];
 	    const char * ccfname = s2ccwp( fname, cfname,
					   maxDataFilenameLength );

	    // Remove constness from kmat.
	    int kmat = const_kmat; 

	    // we don't know the value of nmat until runtime so we
	    // must dynamically allocate a_matids.
	    int *a_matids = new int [ kmat ];

	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------

	    extc_gmatids( ccfname, a_matids, kmat, nmat, ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    // resize and update the vector matids fromt he array version.
	    matids.resize( nmat );
	    std::copy( a_matids, a_matids+nmat, matids.begin() );
	    
	    // Free up dynamic memory and return.
	    delete [] a_matids;

	    return;

	} // end of gmatids


    //----------------------------------------//
    //                gkeys                   //
    //----------------------------------------//
    
    void gkeys( const std::string &fname, const int &const_matid, 
		vector<string> &vkeys,
		const int const_kkeys, int &nkeys, int &ier)
	{
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

	    // we do not know the value of numKeys until after we call 
	    // gkeys() so we create the character array keys[][] to be 
	    // maxKeys long.  This array will later be copied into the
	    // vector vkeys that is returned to the calling program.

	    std::cout << "GandolfWrapper.cc::gkeys()  --> const char * ???" << std::endl;
	    
	    char keys[maxKeys][key_length];

	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------
	    
	    extc_gkeys( ccfname, matid, keys, kkeys, nkeys, ier );

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
		    strncpy( key, keys[i], key_length );
		    // kill trailing whitespace.
		    strtok( key, " " );
		    // store the keyword in a vector.
		    vkeys[i].assign( key, 0, strlen(key) );
		}
	    
	} // end of gkeys


    //----------------------------------------//
    //                gchgrids                //
    //----------------------------------------//
    
    void gchgrids( const std::string &fname, const int &matid,
		   int &nt, int &nrho, int &nhnu, int &ngray, int &nmg,
		   int &ier )
	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------
	    
 	    // copy filename into a const char * array;
	    char cfname[maxDataFilenameLength];
 	    const char * ccfname = s2ccwp( fname, cfname,
					   maxDataFilenameLength );
	    
	    // remove const-ness
	    int nc_matid = matid; // const

	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------
	    
	    extc_gchgrids( ccfname, nc_matid, nt, nrho, nhnu,
			   ngray, nmg, ier );

    } // end of gchgrids

    //----------------------------------------//
    //                ggetgray                //
    //----------------------------------------//
    
    void ggetgray( const string &fname,   const int &const_matid, const string skey,
		   vector<double> &temps, const int &const_kt,    int &nt, 
		   vector<double> &rhos,  const int &const_krho,  int &nrho,
		   vector<double> &data,  const int &const_kgray, int &ngray,
		   int &ier )
	{
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
	    
	    // Allocate memory for double arrays (temps,rhos,data).
	    // These will be copied into vector<double> objects later.
	    double *array_temps = new double [kt];
	    double *array_rhos  = new double [krho];
	    double *array_data  = new double [kgray];
	    
	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------

	    extc_ggetgray( ccfname,     matid, cckey,
			   array_temps, kt,    nt, 
			   array_rhos,  krho,  nrho,
			   array_data,  kgray, ngray,
			   ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    temps.resize(nt);
	    rhos.resize(nrho);
	    data.resize(ngray);

	    std::copy( array_temps, array_temps+nt,   temps.begin() );
	    std::copy( array_rhos,  array_rhos+nrho,  rhos.begin()  );
	    std::copy( array_data,  array_data+ngray, data.begin()  );

	    delete [] array_temps;
	    delete [] array_rhos;
	    delete [] array_data;

	} // end of ggetgray

    //----------------------------------------//
    //                gintgrlog               //
    //----------------------------------------//
    
    void gintgrlog( const vector<double> &temps, const int &const_nt,
		    const vector<double> &rhos,  const int &const_nrho,
		    const vector<double> &data,  const int &const_ngray,
		    const double &const_tlog, const double &const_rlog, 
		    double &ans )
 	{
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

	    extc_gintgrlog( array_temps, nt, array_rhos, nrho,
			    array_data, ngray, tlog, rlog, ans );
	    
	    // no error code is returned from this function.
	    // we don't need to copy any data back into C++ data
	    // types.  The only return value is "ans" and it is
	    // already in the correct format.

	    delete [] array_temps;
	    delete [] array_rhos;
	    delete [] array_data;

	} // end of ginggrlog

    //----------------------------------------//
    //                ggetmg                  //
    //----------------------------------------//
    
    // Read data grid (temp,density,energy_bounds) and mg opacity
    // data.  Retrieve both the size of the data and the actual data.

    void ggetmg( const string &fname,   const int &const_matid, const string skey,
		 vector<double> &temps, const int &const_kt,    int &nt, 
		 vector<double> &rhos,  const int &const_krho,  int &nrho,
		 vector<double> &hnus,  const int &const_khnu,  int &nhnu,
		 vector<double> &data,  const int &const_kdata, int &ndata,
		 int &ier )
	{
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
	    
	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------

	    extc_ggetmg( ccfname, matid, cckey,
			 array_temps, kt,    nt, 
			 array_rhos,  krho,  nrho,
			 array_hnus,  khnu,  nhnu,
			 array_data,  kdata, ndata,
			 ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    // resize data found in the Opacity object.
	    temps.resize(nt);
	    rhos.resize(nrho);
	    hnus.resize(nhnu);
	    data.resize(ndata);

	    std::copy( array_temps, array_temps+nt,   temps.begin() );
	    std::copy( array_rhos,  array_rhos+nrho,  rhos.begin()  );
	    std::copy( array_hnus,  array_hnus+nhnu,  hnus.begin()  );
	    std::copy( array_data,  array_data+ndata, data.begin()  );
	    
	    // free up dynamically allocated memory
	    
	    delete [] array_temps;
	    delete [] array_rhos;
	    delete [] array_hnus;
	    delete [] array_data;

	} // end of ggetgray

    //----------------------------------------//
    //                gintmglog               //
    //----------------------------------------//
    
    void gintmglog( const vector<double> &temps, const int &const_nt,
		    const vector<double> &rhos,  const int &const_nrho,
		    const int &const_nhnu,
		    const vector<double> &data,  const int &const_ndata,
		    const double &const_tlog, const double &const_rlog, 
		    vector<double> &ansmg )
 	{
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

	    extc_gintmglog( array_temps, nt, array_rhos, nrho,
			    nhnu, array_data, ndata, tlog, rlog,
			    array_ansmg ); 
	    
	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    std::copy( array_ansmg, array_ansmg+ngroups, 
		       ansmg.begin() );

	    // release space required by temps;
	    delete [] array_temps;
	    delete [] array_rhos;
	    delete [] array_data;
	    delete [] array_ansmg;

	} // end of ginggrlog

} // end namespace wrapper

} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
//                         end of cdi/GandolfWrapper.cc
//---------------------------------------------------------------------------//
