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

// void string2char( const string &source, char target[], 
// 		  int targetLength )
//     {
// 	// If the target c-string is longer than the source string then 
// 	// copy the string into the c-string and pad the extra cells
// 	// with blanks.	
	
// 	// If the target c-string is shorter than the source string
// 	// then copy targetLength characters from the source string
// 	// into the target string.

// 	if ( targetLength >= source.length() ) {
// 	    for ( int i=0; i < source.length(); ++i )
// 		target[i] = source[i];
// 	    for ( int i=source.length(); i < targetLength; ++i )
// 		target[i] = ' ';
// 	} else {
// 	    for ( int i=0; i < targetLength; ++i )
// 		target[i] = source[i];
// 	}
//     }

 const char *s2ccwp( const string &source, const int n )
     {
	 // remove const-ness
	 string s1(source);
	 // create a string of spaces
	 string padding(n-source.size(), ' ');
	 // append padding to the source string
	 s1.append(padding);
	 // create a c-string of the correct length;
	 char *c1 = new char [n];  // MEMORY LEAK ????????
	 // copy the data from s1 into c1
 	 for ( int i=0; i<n; ++i )
 	     c1[i] = s1[i];
	 // append a string terminator.
	 c1[79] = '\0';
	 // return the new c-string
	 return c1;	 
     }


    //----------------------------------------//
    //                gmatids                 //
    //----------------------------------------//
    
    void gmatids( const std::string &fname , vector<int> &matids, 
		  const int kmat, int &nmat, int &ier ) 
	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------

	    // copy filename into a char array;
	    // also remove constness.
//   	    char cfname[maxDataFilenameLength];
//  	    string2char( fname, cfname, maxDataFilenameLength );

	    // create "long int" versions of variables.
	    long int li_nmat = nmat;
	    long int li_ier  = ier;
	    // also remove constness from kmat.
	    long int li_kmat = kmat; 

	    // we don't know the value of nmat until runtime so we
	    // must dynamically allocate li_matids.
	    long int *li_matids = new long int [ kmat ];

	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------
	    
	    extc_gmatids( s2ccwp(fname,maxDataFilenameLength),
			  li_matids, li_kmat, li_nmat, li_ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    // update the function arguments from their "long int"
	    // counterparts.  We don't update kmat since it is const. 
	    ier  = li_ier;
	    nmat = li_nmat;
	    matids.resize( nmat );
	    std::copy( li_matids, li_matids+nmat, matids.begin() );
	    
	    // Free up dynamic memory and return.
	    delete [] li_matids;

	    return;

	} // end of gmatids


    //----------------------------------------//
    //                gkeys                   //
    //----------------------------------------//
    
    void gkeys( const std::string &fname, const int &matid, 
		vector<string> &vkeys,
		const int kkeys, int &nkeys, int &ier)
	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------
	    
	    // copy filename into a char array;
	    // also remove constness.
// 	    char cfname[maxDataFilenameLength];
// 	    string2char( fname, cfname, maxDataFilenameLength );

	    long int li_matid = matid; // const
	    long int li_kkeys = kkeys; // const
	    long int li_nkeys = nkeys;
	    long int li_ier   = ier;

	    // we do not know the value of numKeys until after we call 
	    // gkeys() so we create the character array keys[][] to be 
	    // maxKeys long.  This array will later be copied into the
	    // vector vkeys that is returned to the calling program.
	    char keys[maxKeys][key_length];

	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------
	    
	    extc_gkeys( s2ccwp(fname,maxDataFilenameLength),
			li_matid, keys, li_kkeys, li_nkeys,
			li_ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    // we don't modify matID or kkeys because these are const
	    // values.
	    nkeys = li_nkeys;	    
	    ier   = li_ier;

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
	    
// 	    // copy filename into a char array;
// 	    // also remove constness.
// 	    char cfname[maxDataFilenameLength];
// 	    string2char( fname, cfname, maxDataFilenameLength );
	    
	    long int li_matid = matid; // const
	    long int li_nt    = nt; 
	    long int li_nrho  = nrho;
	    long int li_nhnu  = nhnu;
	    long int li_ngray = ngray;
	    long int li_nmg   = nmg;
	    long int li_ier   = ier;

	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------
	    
	    extc_gchgrids( s2ccwp(fname,maxDataFilenameLength), 
			   li_matid, li_nt, li_nrho, li_nhnu,
			   li_ngray, li_nmg, li_ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    // copy data back into standard ojects.
	    // we don't modify matID because it is a const value.
	    nt    = li_nt;
	    nrho  = li_nrho;
	    nhnu  = li_nhnu;
	    ngray = li_ngray;
	    nmg   = li_nmg;
	    ier   = li_ier;

    } // end of gchgrids

    //----------------------------------------//
    //                ggetgray                //
    //----------------------------------------//
    
    void ggetgray( const string &fname,   const int &matid, const string skey,
		   vector<double> &temps, const int &kt,    int &nt, 
		   vector<double> &rhos,  const int &krho,  int &nrho,
		   vector<double> &data,  const int &kgray, int &ngray,
		   int &ier )
	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------

	    // copy filename into a char array;
	    // also remove constness.
// 	    char cfname2[maxDataFilenameLength];
// 	    string2char( fname, cfname, maxDataFilenameLength );

// 	    char key[ key_length ];                           
// 	    string2char( skey, key, key_length );
	    
	    // cast all integers as long integers before calling ggetgray.
	    long int li_matid = matid; // const
	    long int li_kt    = kt;    // const
	    long int li_nt    = nt;
	    long int li_krho  = krho;  // const
	    long int li_nrho  = nrho;
	    long int li_kgray = kgray; // const
	    long int li_ngray = ngray;
	    long int li_ier   = ier;
	    
	    // Allocate memory for double arrays (temps,rhos,data).
	    // These will be copied into vector<double> objects later.
	    double *array_temps = new double [kt];
	    double *array_rhos  = new double [krho];
	    double *array_data  = new double [kgray];
	    
	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------

	    extc_ggetgray( s2ccwp(fname,maxDataFilenameLength),
			   li_matid, 
			   s2ccwp(skey,key_length),
			   array_temps, li_kt,    li_nt, 
			   array_rhos,  li_krho,  li_nrho,
			   array_data,  li_kgray, li_ngray,
			   li_ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    nt    = li_nt;
	    nrho  = li_nrho;
	    ngray = li_ngray;
	    ier   = li_ier;

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
    
    void gintgrlog( const vector<double> &temps, const int &nt,
		    const vector<double> &rhos,  const int &nrho,
		    const vector<double> &data,  const int &ngray,
		    const double &const_tlog, const double &const_rlog, 
		    double &ans )
 	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------

	    // cast all integers as long integers before calling ggetgray.
	    long int li_nt    = nt;    // const
	    long int li_nrho  = nrho;  // const
	    long int li_ngray = ngray; // const
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

	    extc_gintgrlog( array_temps, li_nt, array_rhos, li_nrho,
			    array_data, li_ngray, tlog, rlog, ans );
	    
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

    void ggetmg( const string &fname,   const int &matid, const string skey,
		 vector<double> &temps, const int &kt,    int &nt, 
		 vector<double> &rhos,  const int &krho,  int &nrho,
		 vector<double> &hnus,  const int &khnu,  int &nhnu,
		 vector<double> &data,  const int &kdata, int &ndata,
		 int &ier )
	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------

// 	    // copy filename into a char array;
// 	    // also remove constness.
// 	    char cfname[maxDataFilenameLength];
// 	    string2char( fname, cfname, maxDataFilenameLength );

// 	    char key[ key_length ];                           
// 	    string2char( skey, key, key_length);

	    // cast all integers as long integers before calling ggetgray.
	    long int li_matid = matid; // const
	    long int li_kt    = kt;    // const
	    long int li_nt    = nt;
	    long int li_krho  = krho;  // const
	    long int li_nrho  = nrho;
	    long int li_khnu  = khnu;  // const
	    long int li_nhnu  = nhnu;
	    long int li_kdata = kdata; // const
	    long int li_ndata = ndata;
	    long int li_ier   = ier;
	    
	    // Allocate memory for double arrays (temps,rhos,data).
	    // These will be copied into vector<double> objects later.
	    double *array_temps = new double [kt];
	    double *array_rhos  = new double [krho];
	    double *array_hnus  = new double [khnu];
	    double *array_data  = new double [kdata];
	    
	    // --------------------------------------------------
	    // call the Gandolf library function
	    // --------------------------------------------------

	    extc_ggetmg( s2ccwp(fname,maxDataFilenameLength),
			 li_matid, 
			 s2ccwp(skey,key_length),
			 array_temps, li_kt,    li_nt, 
			 array_rhos,  li_krho,  li_nrho,
			 array_hnus,  li_khnu,  li_nhnu,
			 array_data,  li_kdata, li_ndata,
			 li_ier );

	    // ----------------------------------------
	    // Copy the data back into C++ data types
	    // ----------------------------------------

	    nt    = li_nt;
	    nrho  = li_nrho;
	    nhnu  = li_nhnu;
	    ndata = li_ndata;
	    ier   = li_ier;

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
    
    void gintmglog( const vector<double> &temps, const int &nt,
		    const vector<double> &rhos,  const int &nrho,
		    const int &nhnu,
		    const vector<double> &data,  const int &ndata,
		    const double &const_tlog, const double &const_rlog, 
		    vector<double> &ansmg )
 	{
	    // ----------------------------------------
	    // Create simple flat data types
	    // ----------------------------------------

	    int ngroups = nhnu-1;

	    // cast all integers as long integers before calling
	    // ggetgray.  Also remove "const-ness".	    
	    long int li_nt    = nt;    // const
	    long int li_nrho  = nrho;  // const
	    long int li_nhnu  = nhnu;  // const 
	    long int li_ndata = ndata; // const
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

	    extc_gintmglog( array_temps, li_nt, array_rhos, li_nrho,
			    li_nhnu, array_data, li_ndata, tlog, rlog,
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
