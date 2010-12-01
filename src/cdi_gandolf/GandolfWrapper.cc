//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfWrapper.cc
 * \author Kelly Thompson
 * \date   Thu Jul 13 15:31:56 2000
 * \brief  The ANSI-C functions wrap the Fortran functions found in
 *         libgandolf.
 * \note   Copyright (C) 2000-2010 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfWrapper.hh"
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
const char *s2ccwp( string const &source,
                    char         *c1,
                    int           n )
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
	
#ifndef rtt_cdi_gandolf_stub

//----------------------------------------//
//                gmatids                 //
//----------------------------------------//

int wgmatids( string const & fname,
              vector<int>  & matids, 
              int    const   const_kmat,
              int          & nmat ) 
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
} // end of gmatids
			
//----------------------------------------//
//                gkeys                   //
//----------------------------------------//
	
int wgkeys( string const   & fname,
            int    const   & const_matid, 
            vector<string> & vkeys,
            int    const   & const_kkeys,
            size_t         & nkeys_in )
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
    int ier = 0;

    int nkeys( static_cast<int>(nkeys_in) );
                
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

    gkeys( ccfname, matid, bjkeys, kkeys, nkeys, ier );
		
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
} // end of gkeys

	
//----------------------------------------//
//                gchgrids                //
//----------------------------------------//
	
int wgchgrids( string const & fname,
               int    const & const_matid, 
               int          & nt,
               int          & nrho,
               int          & nhnu,
               int          & ngray, 
               int          & nmg )
{
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
} // end of gchgrids


	
//----------------------------------------//
//                ggetgray                //
//----------------------------------------//
	
int wggetgray( string const   & fname,   
               int    const   & const_matid,
               string const     skey,
               vector<double> & temps,
               int    const   & const_kt,
               vector<double> & rhos,
               int    const   & const_krho,
               vector<double> & data,
               int    const   & const_kgray )
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
} // end of ggetgray

	
//----------------------------------------//
//                gintgrlog               //
//----------------------------------------//
	
double wgintgrlog( vector<double> const & temps,
                   int            const & const_nt,
                   vector<double> const & rhos,
                   int            const & const_nrho,
                   vector<double> const & data,
                   int            const & const_ngray,
                   double         const & const_tlog,
                   double         const & const_rlog ) 
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
} // end of ginggrlog
	

//----------------------------------------//
//                ggetmg                  //
//----------------------------------------//
	
// Read data grid (temp,density,energy_bounds) and mg opacity
// data.  Retrieve both the size of the data and the actual data.
	
int wggetmg( string const   & fname,  
             int    const   & const_matid,
             string const     skey,
             vector<double> & temps,
             int    const   & const_kt,
             vector<double> & rhos,
             int    const   & const_krho,
             vector<double> & hnus,
             int    const   & const_khnu,
             vector<double> & data,
             int    const   & const_kdata )
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
} // end of wggetmg
	

//----------------------------------------//
//                gintmglog               //
//----------------------------------------//
	
vector<double>
wgintmglog( vector<double> const & temps,
            int            const & const_nt,
            vector<double> const & rhos,
            int            const & const_nrho,
            int            const & const_nhnu,
            vector<double> const & data,
            int            const & const_ndata,
            double         const & const_tlog,
            double         const & const_rlog )
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
} // end of wgintmglog



// -------------------------------------------------------------------------- //    
// Stubbed out routines
// -------------------------------------------------------------------------- //    
#else // ifndef rtt_cdi_gandolf_stub

int wgmatids( string const &,
              vector<int>  &, 
              int    const,
              int          & ) { return 1; }
int wgkeys( string const   &,
            int    const   &, 
            vector<string> &,
            int    const   &,
            size_t         & ) { return 1; }
int wgchgrids( string const &,
               int    const &, 
               int          &,
               int          &,
               int          &,
               int          &, 
               int          & ) { return 1; }
int wggetgray( string const   &,   
               int    const   &,
               string const,
               vector<double> &,
               int    const   &,
               vector<double> &,
               int    const   &,
               vector<double> &,
               int    const   & ) { return 1; }
double wgintgrlog( vector<double> const &,
                   int            const &,
                   vector<double> const &,
                   int            const &,
                   vector<double> const &,
                   int            const &,
                   double         const &,
                   double         const & ) { return 0.0; }
int wggetmg( string   const &,  
             int      const &,
             string   const,
             vector<double> &,
             int      const &,
             vector<double> &,
             int      const &,
             vector<double> &,
             int      const &,
             vector<double> &,
             int      const & ) { return 1; }
vector<double>
wgintmglog( vector<double> const &,
            int            const &,
            vector<double> const &,
            int            const &,
            int            const & const_nhnu,
            vector<double> const &,
            int            const &,
            double         const &,
            double         const & )
{
    vector<double> ansmg( const_nhnu - 1, 0.0 );
    return ansmg; // dummy vector.
}

#endif // ifndef rtt_cdi_gandolf_stub
	
	
} // end namespace wrapper
} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
//                         end of cdi/GandolfWrapper.cc
//---------------------------------------------------------------------------//
