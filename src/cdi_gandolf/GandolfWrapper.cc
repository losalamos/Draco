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

#include <string>
#include <iostream>
#include <iomanip>

namespace rtt_cdi_gandolf {

    using std::cout;
    using std::endl;
    using std::string;

    //----------------------------------------//
    //                gmatids                 //
    //----------------------------------------//
    
    void gmatids( const std::string &fname , int matids[], 
		  const int kmat, int &nmat, int &ier ) 
	{
	    //	    cout << "In GandolfWrapper::gmatids()" << endl;

	    // Gandolf will not accept a filename that is longer that
	    // 80 chars in length.	  
	    int len = fname.length();
	    Require( len<80 ); 

	    // copy filename into a char array;
	    // also remove constness.
	    char cfname[80];
	    int i = 0;
	    for ( ; i < len; ++i )
		cfname[i] =fname[i];
	    for ( ; i < 80; ++i )
		cfname[i] = ' ';

	    // create "long int" versions of variables.
	    long int li_nmat = static_cast<long int>(nmat);
	    long int li_ier  = static_cast<long int>(ier);
	    // also remove constness from kmat.
	    long int li_kmat = static_cast<long int>(kmat); 

	    // we don't know the value of kmat until runtime so we
	    // must dynamically allocate li_matids.
	    long int *li_matids = new long int [ kmat ];
	    for ( int i=0; i<kmat; ++i )
		li_matids[i] = static_cast<long int>(matids[i]);

// 	    cout << "We have the following values before calling extc_gmatids:" << endl
// 		 << "   cfname = ";
// 	    for ( int i=0; i<len; ++i ) 
// 		cout << cfname[i];
// 	    cout << endl;
// 	    for ( int i=0; i<li_kmat; ++i)
// 		cout << "   li_matids[" << i << "] = " << li_matids[i] << endl;
// 	    cout << "   li_kmat = " << li_kmat << endl
// 		 << "   li_nmat = " << li_nmat << endl
// 		 << "   li_ier  = " << li_ier  << endl
// 		 << endl;



	    // call the Gandolf library function
	    extc_gmatids( cfname, li_matids, li_kmat, li_nmat, li_ier );



	    //	    cout << "     Back from extc_gmatids:" << endl;
// 	    cout << "   cfname = ";
// 	    for ( int i=0; i<len; ++i ) 
// 		cout << cfname[i];
// 	    cout << endl;
// 	    for ( int i=0; i<li_kmat; ++i)
// 		cout << "   li_matids[" << i << "] = " << li_matids[i] << endl;
// 	    cout << "   li_kmat = " << li_kmat << endl
// 		 << "   li_nmat = " << li_nmat << endl
// 		 << "   li_ier  = " << li_ier  << endl
// 		 << endl;

	    // Abort if Gandolf returns an error.

	    switch ( li_ier ) {
	    case 0: // no errors
		break;
	    case 1: // IPCRESS file not found.
		Insist( false, "The IPCRESS file was not found.");
		break;
	    case 2: // File is not IPCRESS.
		Insist( false, "The file does not appear to be in IPCRESS format");
		break;
	    case 3: // Problem reading file
		Insist( false, "Having trouble reading the IPCRESS file.");
		break;
	    case 4: // No material ID's found in file.
		Insist( false, "No material ID's were found in the IPCRESS data file.");
		break;
	    case 5: // too many matids found ( nmat > kmat )
		Insist( false, "Too many materials were found in the data file ( nmat > kmat ).");
		break;
	    default: // unknown error.
		Insist( false, "Unknown error returned from Gandolf::gmatids().");
		break;
	    }


	    // update the function arguments from their "long int"
	    // counterparts.
	    for ( int i=0; i<kmat; ++i )
		matids[i] = static_cast<int>(li_matids[i]);
	    ier  = static_cast<int>(li_ier);
	    nmat = static_cast<int>(li_nmat);
	    // we don't update kmat since it is const.
	    
	    delete [] li_matids;

	    return;

	} // end of gmatids


    //----------------------------------------//
    //                gmatids                 //
    //----------------------------------------//
    
    void gkeys( const std::string &fname, const int &matid, 
		char keys[][key_length],
		const int kkeys, int &nkeys, int &ier)
	{
	    //	    cout << "In GandolfWrapper::gkeys()." << endl;
	    
	    // Gandolf will not accept a filename that is longer that
	    // 80 chars in length.	  
	    int len = fname.length();
	    Require( len<80 ); 

	    // copy filename into a char array;
	    // also remove constness.
	    char cfname[80];
	    int i = 0;
	    for ( ; i < len; ++i )
		cfname[i] =fname[i];
	    for ( ; i < 80; ++i )
		cfname[i] = ' ';
	    
	    long int li_matid = static_cast<long int>(matid); // const
	    long int li_kkeys = static_cast<long int>(kkeys); // const
	    long int li_nkeys = static_cast<long int>(nkeys);
	    long int li_ier   = static_cast<long int>(ier);

	    // call the Gandolf library function
	    extc_gkeys( cfname, li_matid, keys, li_kkeys, li_nkeys,
			li_ier );

	    //	    cout << "     Back from extc_gkeys:" << endl;

	    // Abort if Gandolf returns an error.
	    switch ( li_ier ) {
	    case 0: // no errors
		break;
	    case 1: // IPCRESS file not found.
		Insist( false, "The IPCRESS file was not found.");
		break;
	    case 2: // File is not IPCRESS.
		Insist( false, "The file does not appear to be in IPCRESS format");
		break;
	    case 3: // Problem reading file
		Insist( false, "Having trouble reading the IPCRESS file.");
		break;
	    case 4: // No keys found for this material.
 		Insist( false, "No keys were found for this material");
		break;
	    case 5: // Too many keys found.
		Insist( false, "Too many keys for array ( nkeys > kkeys )." );
		break;
	    default: // unknown error.
		Insist( false, "Unknown error returned from Gandolf::gmatids().");
		break;
	    }

	    // copy data back into standard ojects.
	    // we don't modify matID or kkeys because these are const
	    // values.
	    nkeys = static_cast<int>(li_nkeys);	    
	    ier   = static_cast<int>(li_ier);

// 	    cout << "   cfname = ";
// 	    for ( int i=0; i<len; ++i ) 
// 		cout << cfname[i];
// 	    cout << endl;
// 	    cout << "   li_matid = " << li_matid << endl;
// 	    cout << "   li_kkeys = " << li_kkeys << endl
// 		 << "   li_nkeys = " << li_nkeys << endl
// 		 << "   li_ier   = " << li_ier  << endl
// 		 << endl;
// 	    for ( int i=0; i<nkeys; ++i ) {
// 		    cout << "   key[" << i << "] = ";
// 		    for ( int j=0; j<key_length; ++j)
// 			cout << keys[i][j];
// 		    cout << endl; }
// 	    cout << endl;


	} // end of gkeys


    //----------------------------------------//
    //                gchgrids                //
    //----------------------------------------//
    
    void gchgrids( const std::string &fname, const int &matid,
		   int &nt, int &nrho, int &nhnu, int &ngray, int &nmg,
		   int &ier )
	{
	    //	    cout << "In GandolfWrapper::gchgrids()." << endl;
	    
	    // Gandolf will not accept a filename that is longer that
	    // 80 chars in length.	  
	    int len = fname.length();
	    Require( len<80 ); 

	    // copy filename into a char array;
	    // also remove constness.
	    char cfname[80];
	    int i = 0;
	    for ( ; i < len; ++i )
		cfname[i] =fname[i];
	    for ( ; i < 80; ++i )
		cfname[i] = ' ';
	    
	    long int li_matid = static_cast<long int>(matid); // const
	    long int li_nt    = static_cast<long int>(nt); 
	    long int li_nrho  = static_cast<long int>(nrho);
	    long int li_nhnu  = static_cast<long int>(nhnu);
	    long int li_ngray = static_cast<long int>(ngray);
	    long int li_nmg   = static_cast<long int>(nmg);
	    long int li_ier   = static_cast<long int>(ier);

	    // call the Gandolf library function
	    extc_gchgrids( cfname, li_matid, li_nt, li_nrho, li_nhnu,
			   li_ngray, li_nmg, li_ier );

	    //	    cout << "     Back from extc_gchgrids:" << endl;

	    // Abort if Gandolf returns an error.
	    switch ( li_ier ) {
	    case 0: // no errors
		break;
	    case -1: // return with etas, not densities
		Insist( false, "IPCRESS file returned ETAs not Densities.");
		break;		
	    case 1: // IPCRESS file not found.
		Insist( false, "The IPCRESS file was not found.");
		break;
	    case 2: // File is not IPCRESS.
		Insist( false, "The file does not appear to be in IPCRESS format");
		break;
	    case 3: // Problem reading file
		Insist( false, "Having trouble reading the IPCRESS file.");
		break;
	    case 4: // Inconsistent gray grids, mg not checked
 		Insist( false, "Gray grid inconsistent with the temp/density grid.");
		break;
	    case 5: // ngray != nt*nrho, mg not checked
		Insist( false, "Wrong number of gray opacities found (ngray != nt*nrho)." );
		break;
	    case 6: // inconsistent mg grid.
		Insist( false, "MG grid inconsistent with the temp/density/hnu grid.");
		break;
	    case 7: //  nmg != nt*nrho*(nhnu-1).
		Insist( false, "Wrong number of MG opacities found (nmg != nt*nrho*(nhnu-1)).");
		break;
	    default: // unknown error.
		Insist( false, "Unknown error returned from Gandolf::gmatids().");
		break;
	    }

	    // copy data back into standard ojects.
	    // we don't modify matID it is a const value.
	    nt    = static_cast<int>(li_nt);
	    nrho  = static_cast<int>(li_nrho);
	    nhnu  = static_cast<int>(li_nhnu);
	    ngray = static_cast<int>(li_ngray);
	    nmg   = static_cast<int>(li_nmg);
	    ier   = static_cast<int>(li_ier);

// 	    cout << "   cfname = ";
// 	    for ( int i=0; i<len; ++i ) 
// 		cout << cfname[i];
// 	    cout << endl;
// 	    cout << "   li_matid = " << li_matid << endl;
// 	    cout << "   li_nt    = " << li_nt    << endl
// 		 << "   li_nrho  = " << li_nrho  << endl
// 		 << "   li_nhnu  = " << li_nhnu  << endl
// 		 << "   li_ngray = " << li_ngray << endl
// 		 << "   li_nmg   = " << li_nmg   << endl
// 		 << "   li_ier   = " << li_ier   << endl
// 		 << endl;
// 	    cout << endl;

    } // end of gchgrids

    //----------------------------------------//
    //                ggetgray                //
    //----------------------------------------//
    
    void ggetgray( const string &fname,   const int &matid, char *key, 
		   vector<double> &temps, const int &kt,    int &nt, 
		   vector<double> &rhos,  const int &krho,  int &nrho,
		   vector<double> &data,  const int &kgray, int &ngray,
		   int &ier )
	{
	    //	    cout << "In GandolfWrapper::ggetgray()." << endl;
	    
	    // Gandolf will not accept a filename that is longer that
	    // 80 chars in length.	  
	    int len = fname.length();
	    Require( len<80 ); 

	    // copy filename into a char array;
	    // also remove constness.
	    char cfname[80];
	    int i = 0;
	    for ( ; i < len; ++i )
		cfname[i] =fname[i];
	    for ( ; i < 80; ++i )
		cfname[i] = ' ';

	    // cast all integers as long integers before calling ggetgray.
	    long int li_matid = static_cast<long int>(matid); // const
	    long int li_kt    = static_cast<long int>(kt);    // const
	    long int li_nt    = static_cast<long int>(nt);
	    long int li_krho  = static_cast<long int>(krho);  // const
	    long int li_nrho  = static_cast<long int>(nrho);
	    long int li_kgray = static_cast<long int>(kgray); // const
	    long int li_ngray = static_cast<long int>(ngray);
	    long int li_ier   = static_cast<long int>(ier);
	    
	    // Allocate memory for double arrays (temps,rhos,data).
	    // These will be copied into vector<double> objects later.
	    double *array_temps = new double [kt];
	    double *array_rhos  = new double [krho];
	    double *array_data  = new double [kgray];
	    
	    // Call the routine from the Gandolf library.
	    extc_ggetgray( cfname,      li_matid, key, 
			   array_temps, li_kt,    li_nt, 
			   array_rhos,  li_krho,  li_nrho,
			   array_data,  li_kgray, li_ngray,
			   li_ier );
	    
	    //	    cout << "     Back from extc_ggetgray:" << endl;

	    // abort if Gandolf returned an error.
	    switch ( li_ier ) {
	    case 0: // no errors
		break;
	    case -1: // return with etas, not densities
		Insist( false, "IPCRESS file returned ETAs not Densities.");
		break;		
	    case 1: // IPCRESS file not found.
		Insist( false, "The IPCRESS file was not found.");
		break;
	    case 2: // File is not IPCRESS.
		Insist( false, "The file does not appear to be in IPCRESS format");
		break;
	    case 3: // Problem reading file
		Insist( false, "Having trouble reading the IPCRESS file.");
		break;
	    case 4: // Data not found
 		Insist( false, "Requested data not found.  Check nt, nrho, ngray.");
		break;
	    case 5: // Data larger than allocated arrays.
		Insist( false, "Data found is larger than allocated array size." );
		break;
	    case 6: // Data size not equal to nt*nrho
		Insist( false, "Data size not equal to expected size (ndata != nt*nrho)");
		break;
	    default: // unknown error.
		Insist( false, "Unknown error returned from Gandolf::gmatids().");
		break;
	    }

	    // copy flat data structures back into 
	    nt    = static_cast<int>(li_nt);
	    nrho  = static_cast<int>(li_nrho);
	    ngray = static_cast<int>(li_ngray);
	    ier   = static_cast<int>(li_ier);

	    temps.resize(nt);
	    rhos.resize(nrho);
	    data.resize(ngray);

	    std::copy( array_temps, array_temps+nt,   temps.begin() );
	    std::copy( array_rhos,  array_rhos+nrho,  rhos.begin()  );
	    std::copy( array_data,  array_data+ngray, data.begin()  );

// 	    cout << endl
// 		 << "Tabulated values from the IPCRESS file:" << endl;
// 	    cout << endl;
// 	    for( int i=0; i<nt; ++i)
// 		cout << "temps[" << i << "] = " 
// 		     << std::scientific << std::setprecision(2)
// 		     << temps[i] << " keV" << endl;
// 	    cout << endl;
// 	    for( int i=0; i<nrho; ++i)
// 		cout << "rhos[" << i << "] = " 
// 		     << std::setw(4)
// 		     << rhos[i] << " g/cm^3" << endl;
// 	    cout << endl;
// 	    for( int i=0; i<ngray; ++i)
// 		cout << "rss[" << i << "] = " 
// 		     << std::scientific << std::setprecision(8)
// 		     << data[i] << " cm^2/g" << endl;
// 	    cout << endl;

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
	    //	    cout << "In GandolfWrapper::gintgrlog()." << endl;
	    
	    // cast all integers as long integers before calling ggetgray.
	    long int li_nt    = static_cast<long int>(nt);    // const
	    long int li_nrho  = static_cast<long int>(nrho);  // const
	    long int li_ngray = static_cast<long int>(ngray); // const
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
	    
	    // Call the routine from the Gandolf library.
	    extc_gintgrlog( array_temps, li_nt, array_rhos, li_nrho,
			    array_data, li_ngray, tlog, rlog, ans );
	    
	    //	    cout << "     Back from extc_gintgrlog:" << endl;

	    // Gandolf does not return an error code for this function.

	    // copy flat data structures back into 
	    // (none??? -- maybe ans???)

	    delete [] array_temps;
	    delete [] array_rhos;
	    delete [] array_data;

	} // end of ginggrlog

    //----------------------------------------//
    //                ggetmg                  //
    //----------------------------------------//
    
    // Read data grid (temp,density,energy_bounds) and mg opacity
    // data.  Retrieve both the size of the data and the actual data.

    void ggetmg( const string &fname,   const int &matid, char *key, 
		 vector<double> &temps, const int &kt,    int &nt, 
		 vector<double> &rhos,  const int &krho,  int &nrho,
		 vector<double> &hnus,  const int &khnu,  int &nhnu,
		 vector<double> &data,  const int &kdata, int &ndata,
		 int &ier )
	{
	    //	    cout << "In GandolfWrapper::ggetmg()." << endl;
	    
	    // Gandolf will not accept a filename that is longer that
	    // 80 chars in length.	  
	    int len = fname.length();
	    Require( len<80 ); 

	    // copy filename into a char array;
	    // also remove constness.
	    char cfname[80];
	    int i = 0;
	    for ( ; i < len; ++i )
		cfname[i] =fname[i];
	    for ( ; i < 80; ++i )
		cfname[i] = ' ';

	    // cast all integers as long integers before calling ggetgray.
	    long int li_matid = static_cast<long int>(matid); // const
	    long int li_kt    = static_cast<long int>(kt);    // const
	    long int li_nt    = static_cast<long int>(nt);
	    long int li_krho  = static_cast<long int>(krho);  // const
	    long int li_nrho  = static_cast<long int>(nrho);
	    long int li_khnu  = static_cast<long int>(khnu);  // const
	    long int li_nhnu  = static_cast<long int>(nhnu);
	    long int li_kdata = static_cast<long int>(kdata); // const
	    long int li_ndata = static_cast<long int>(ndata);
	    long int li_ier   = static_cast<long int>(ier);
	    
	    // Allocate memory for double arrays (temps,rhos,data).
	    // These will be copied into vector<double> objects later.
	    double *array_temps = new double [kt];
	    double *array_rhos  = new double [krho];
	    double *array_hnus  = new double [khnu];
	    double *array_data  = new double [kdata];
	    
	    // Call the routine from the Gandolf library.
	    extc_ggetmg( cfname,      li_matid, key, 
			 array_temps, li_kt,    li_nt, 
			 array_rhos,  li_krho,  li_nrho,
			 array_hnus,  li_khnu,  li_nhnu,
			 array_data,  li_kdata, li_ndata,
			 li_ier );
	    
	    //	    cout << "     Back from extc_ggetmg:" << endl;

	    // abort if Gandolf returned an error.
	    switch ( li_ier ) {
	    case 0: // no errors
		break;
	    case -1: // return with etas, not densities
		Insist( false, "IPCRESS file returned ETAs not Densities.");
		break;		
	    case 1: // IPCRESS file not found.
		Insist( false, "The IPCRESS file was not found.");
		break;
	    case 2: // File is not IPCRESS.
		Insist( false, "The file does not appear to be in IPCRESS format");
		break;
	    case 3: // Problem reading file
		Insist( false, "Having trouble reading the IPCRESS file.");
		break;
	    case 4: // Data not found
 		Insist( false, "Requested data not found.  Check nt, nrho, nhnu and ndata.");
		break;
	    case 5: // Data larger than allocated arrays.
		Insist( false, "Data found is larger than allocated array size." );
		break;
	    case 6: // Data size not equal to nt*nrho
		Insist( false, "Data size not equal to expected size (ndata != nt*nrho*(nhnu-1))");
		break;
	    default: // unknown error.
		Insist( false, "Unknown error returned from Gandolf::gmatids().");
		break;
	    }

	    // copy flat data structures back into 
	    nt    = static_cast<int>(li_nt);
	    nrho  = static_cast<int>(li_nrho);
	    nhnu  = static_cast<int>(li_nhnu);
	    ndata = static_cast<int>(li_ndata);
	    ier   = static_cast<int>(li_ier);

	    // resize data found in the Opacity object.
	    temps.resize(nt);
	    rhos.resize(nrho);
	    hnus.resize(nhnu);
	    data.resize(ndata);

	    std::copy( array_temps, array_temps+nt,   temps.begin() );
	    std::copy( array_rhos,  array_rhos+nrho,  rhos.begin()  );
	    std::copy( array_hnus,  array_hnus+nhnu,  hnus.begin()  );
	    std::copy( array_data,  array_data+ndata, data.begin()  );

// 	    cout << endl
// 		 << "Tabulated values from the IPCRESS file:" << endl;
// 	    cout << endl;
// 	    for( int i=0; i<nt; ++i)
// 		cout << "temps[" << i << "] = " 
// 		     << std::scientific << std::setprecision(2)
// 		     << temps[i] << " keV" << endl;
// 	    cout << endl;
// 	    for( int i=0; i<nrho; ++i)
// 		cout << "rhos[" << i << "] = " 
// 		     << std::setw(4)
// 		     << rhos[i] << " g/cm^3" << endl;
// 	    cout << endl;
// 	    for( int i=0; i<ndata; ++i)
// 		cout << "rss[" << i << "] = " 
// 		     << std::scientific << std::setprecision(8)
// 		     << data[i] << " cm^2/g" << endl;
// 	    cout << endl;

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
	    //	    cout << "In GandolfWrapper::gintmglog()." << endl;
	    
	    int ngroups = nhnu-1;

	    // cast all integers as long integers before calling
	    // ggetgray.  Also remove "const-ness".	    
	    long int li_nt    = static_cast<long int>(nt);    // const
	    long int li_nrho  = static_cast<long int>(nrho);  // const
	    long int li_nhnu  = static_cast<long int>(nhnu);  // const 
	    long int li_ndata = static_cast<long int>(ndata); // const
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
	    
	    // Call the routine from the Gandolf library.
	    extc_gintmglog( array_temps, li_nt, array_rhos, li_nrho,
			    li_nhnu, array_data, li_ndata, tlog, rlog,
			    array_ansmg ); 
	    
	    //	    cout << "     Back from extc_gintmglog:" << endl;

	    // Gandolf does not return an error code for this function.

	    // copy flat data structures back into standard data
	    // types.
	    std::copy( array_ansmg, array_ansmg+ngroups, 
		       ansmg.begin() );

	    // release space required by temps;
	    delete [] array_temps;
	    delete [] array_rhos;
	    delete [] array_data;
	    delete [] array_ansmg;

	} // end of ginggrlog


} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
//                         end of cdi/GandolfWrapper.cc
//---------------------------------------------------------------------------//
