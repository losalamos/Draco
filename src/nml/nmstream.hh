//---------------------------------*-C++-*---------------------------------//
// nmstream.hh
// Geoffrey Furnish
// 3 June 1992
//-------------------------------------------------------------------------//
// @> Definitions of abstractions for stream based i/o.  Implement
// @> this interface on your favorite hw.
//-------------------------------------------------------------------------//

#ifndef __nml_nmstream_hh__
#define __nml_nmstream_hh__

#include "ds++/String.hh"

#include "nml/Group.hh"

#define ATT_io
#include <iostream.h>
#include <fstream.h>

class nmistream {

#ifdef ATT_io
    ifstream inf;
#endif

  public:
    nmistream( const char *name, char *mode ="r" );
    ~nmistream(void);
    void close(void);

    int get( char& c );
    int getline( char *buff, int len );
    void putback( char c );

    nmistream& operator>>( NML_Block& b );
    nmistream& operator>>( NML_Group& g );
};

class nmostream {

#ifdef ATT_io
    ofstream of;
#endif

  public:
    nmostream( const char *name, char *mode );
    ~nmostream();
    void close();

    nmostream& operator<<( const char& c );
    nmostream& operator<<( const char *s );

    nmostream& operator<<( const NML_Block& b );
    nmostream& operator<<( const NML_Group& g );

    nmostream& operator<<( const dsxx::String& s );
};

#endif				// __nml_nmstream_hh__

//-------------------------------------------------------------------------//
//                              end of nml/nmstream.hh
//-------------------------------------------------------------------------//
