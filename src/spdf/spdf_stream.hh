//----------------------------------*-C++-*----------------------------------//
// spdf_stream.hh
// Geoffrey Furnish
// Thu Aug 10 23:26:33 1995
//---------------------------------------------------------------------------//
// @> C++ implementation of Maurice's SPDF.
//---------------------------------------------------------------------------//

#ifndef __spdf_spdf_stream_hh__
#define __spdf_spdf_stream_hh__

#include <stdio.h>

#include <algorithm>

#include "spdf/config.hh"

/* Error numbers */

#define SPDF_ERROR		 1	/* Unknown error 		*/
#define SPDF_FNOPEN		 2	/* File not open 		*/
#define SPDF_FAOPEN		 3	/* File already open 		*/
#define SPDF_BADUN		 4	/* Bad unit number 		*/
#define SPDF_BADNBITS		 5	/* Invalid # of bits 		*/
#define SPDF_RDERR		 6	/* Read error 			*/
#define SPDF_WRERR		 7	/* Write error 			*/
#define SPDF_NOTSPDF		 8	/* Not a valid SPDF file 	*/
#define SPDF_BADLEN		 9	/* Record too long 		*/
#define SPDF_WRONGTYPE		10	/* Record of wrong type 	*/
#define SPDF_EOF		11	/* Unexpected end of file 	*/

//===========================================================================//
// class spdf_stream - Simple Portable Data File stream

// This class defines the basic capabilities of the simple portable data file
// facility.  However, the actual i/o methods are pure virtual.  You have to
// derive from this class and implement the low level i/o routines in order
// to implement spdf for the device of your choice.
//===========================================================================//

class spdf_stream : public spdf_types {

  protected:
    int bp;
    int debug;

  public:
    spdf_stream();
    virtual ~spdf_stream() {}

    void dbug_enter( char *m ) { m++; /* :-) */ }
    void dbug_on()  { debug = 1; }
    void dbug_off() { debug = 0; }

    void clear() { bp=0; }
    int  size() const { return bp; }

    int sizeof_int()   { return 4; }
    int sizeof_float() { return 4; }

// The next set of stuff comes mostly from PLplot's pdf.

    virtual int spdf_putc( int c )   =0;
    virtual int spdf_getc()          =0;
    virtual int spdf_ungetc( int c ) =0;
    virtual void wrx( const U_CHAR *x, int nitems ) =0;
    virtual void rdx( U_CHAR *x, int nitems ) =0;

    int  wr_header( char *header );
    int  rd_header( char *header );
    void wr_1byte( U_CHAR  s );
    void rd_1byte( U_CHAR& s );
    void wr_2bytes( U_SHORT  s );
    void rd_2bytes( U_SHORT& s );
    void wr_2nbytes( U_SHORT *s, int n );
    void rd_2nbytes( U_SHORT *s, int n );
    void wr_4bytes( u32  s );
    void rd_4bytes( u32& s );
    void wr_8bytes( u64  s );
    void rd_8bytes( u64& s );

    void wr_ieee_bytes( u32& f ) { wr_4bytes(f); }
    void wr_ieee_bytes( u64& d ) { wr_8bytes(d); }

    void rd_ieee_bytes( u32& f ) { rd_4bytes(f); }
    void rd_ieee_bytes( u64& d ) { rd_8bytes(d); }

    template<class T> void wr_ieee( T  d );
    template<class T> void rd_ieee( T& d );

    void print_ieeef( void *vx, void *vy );
    void print_ieeed( void *vx, void *vy );

// These things come basically from the old utils/spdf suite.

    U_CHAR getbid();

    int eod();

    void match_byte_id( int bid )
    {
	U_CHAR rbid = getbid();
	if (rbid != bid) {
	    if (rbid == 0)
		throw "spdf::end of file";
	    else {
		spdf_ungetc( rbid );
		throw "Wrong byte id";
	    }
	}
    }

    void wr_rep( int v )    { wr_4bytes(v); }
    void wr_rep( float v )  { wr_ieee(v); }
    void wr_rep( double v ) { wr_ieee(v); }

    void rd_rep( int& v )    { u32 x; rd_4bytes(x); v = static_cast<int>(x); }
    void rd_rep( float& v )  { rd_ieee( v ); }
    void rd_rep( double& v ) { rd_ieee( v ); }

    template<class T>
    void wr( T v )
    {
	spdf_putc( spdf_traits<T>::scalar_byte_id );
	wr_rep( v );
    }

    template<class T>
    void rd( T& v )
    {
	match_byte_id( spdf_traits<T>::scalar_byte_id );
	rd_rep( v );
    }

    template<class T>
    void wr( T *v, int len )
    {
	spdf_putc( spdf_traits<T>::record_byte_id );
	wr_rep( len );
	for( int i=0; i < len; i++ )
	    wr_rep( v[i] );
    }

    template<class T> void wrp( T *v, int len,
				int nbits = spdf_traits<T>::compress_50 );
    template<class T> void rdp( T *v, int& len, int max );

    template<class T>
    void rd( T *v, int& len, int max )
    {
	U_CHAR rbid = getbid();
	if (rbid == spdf_traits<T>::record_byte_id) {
	    u32 ul_len;
	    rd_4bytes( ul_len );
	    len = static_cast<int>(ul_len);
	    if (len > max)
		throw "bad length";
	    for( int i=0; i < len; i++ )
		rd_rep( v[i] );
	}
	else if (rbid == spdf_traits<T>::packed_record_byte_id) {
	    spdf_ungetc( rbid );
	    rdp( v, len, max );
	}
	else if (rbid == 0) {
	    throw "spdf::end of file";
	}
	else {
	    spdf_ungetc(rbid);
	    throw "Wrong byte id";
	}
    }

// This is a lot like the template above, except int's don't support
// packing. 

    void rd( int *v, int& len, int max )
    {
	u32 ul_len;
	match_byte_id( spdf_traits<int>::record_byte_id );
	rd_4bytes( ul_len );
	len = static_cast<int>(ul_len);
	if (len > max) throw "bad length";
	for( int i=0; i < len; i++ )
	    rd_rep( v[i] );
    }

  private:
    void packw( U_CHAR *p, int nbits, u32 *u, int ilim  );
    void unpackw( U_CHAR *p, int nbits, u32 *u, int ilim );

    template<class T> void minmax( T *v, int len, T& min, T& max )
    {
	T *pmin = std::min_element( v, v+len );
	T *pmax = std::max_element( v, v+len );
	min = *pmin;
	max = *pmax;
    }
};

//===========================================================================//
// class stdio_spdf_stream - An spdf_stream for i/o on FILE *

// This class implements the pure virtual methods of spdf_stream, threby
// providing an implementation over stdio.
//===========================================================================//

class stdio_spdf_stream : public spdf_stream
{
    FILE *fp;

  public:
    stdio_spdf_stream( const char *n, const char *a );
    ~stdio_spdf_stream();

    int spdf_putc( int c );
    int spdf_getc();
    int spdf_ungetc( int c );
    void wrx( const U_CHAR *x, int nitems );
    void rdx( U_CHAR *x, int nitems );
};

//===========================================================================//
// class mbuf_spdf_stream - An spdf_stream for i/o to a memory buffer

// This class implements the pure virtual methods of spdf_stream, threby
// providing an implementation of spdf using a memory buffer as the "storage
// device".  This is primarily useful for exploiting the SPDF data packing
// capabilities in contexts where the actual i/o is handled by some other
// facility. 
//===========================================================================//

class mbuf_spdf_stream : public spdf_stream
{
    U_CHAR *v;
    int bufmax, may_delete;

  public:
    mbuf_spdf_stream( int _bufmax =2048 );
    mbuf_spdf_stream( U_CHAR *_v, int _bufmax );
    ~mbuf_spdf_stream();

    void reset() { bp = 0; }
    U_CHAR *buf() { return v; }

    void set_min_size( int s );

    int spdf_putc( int c );
    int spdf_getc();
    int spdf_ungetc( int c );
    void wrx( const U_CHAR *x, int nitems );
    void rdx( U_CHAR *x, int nitems );
};

#include <iostream.h>
#include <fstream.h>

//===========================================================================//
// class spdf_ifstream - Do spdf input from an ifstream

// This class makes it easy to do spdf input from an existing ifstream.  What
// is significant about this is that you do not have to manage all aspects of
// the file i/o through spdf.  You can read spdf records right out of the
// middle of an ifstream.
//===========================================================================//

class spdf_ifstream : public spdf_stream
{
    ifstream& ifs;

  public:
    spdf_ifstream( ifstream& _ifs ) : ifs(_ifs) {}

    int spdf_putc( int c );
    int spdf_getc();
    int spdf_ungetc( int c );
    void wrx( const U_CHAR *x, int nitems );
    void rdx( U_CHAR *x, int nitems );
};

//===========================================================================//
// class spdf_ofstream - Do spdf output to an ofstream

// This class makes it easy to do spdf output to an existing ofstream.  What
// is significant about this is that you do not have to manage all aspects of
// the file i/o through spdf.  You can write spdf records right into the
// middle of an ofstream.  This makes it a lot easier to use intelligent file
// structures in C++, and still use spdf.
//===========================================================================//

class spdf_ofstream : public spdf_stream
{
    ofstream& ofs;

  public:
    spdf_ofstream( ofstream& _ofs ) : ofs(_ofs) {}

    int spdf_putc( int c );
    int spdf_getc();
    int spdf_ungetc( int c );
    void wrx( const U_CHAR *x, int nitems );
    void rdx( U_CHAR *x, int nitems );
};

#endif                          // __spdf_spdf_stream_hh__

//---------------------------------------------------------------------------//
//                              end of spdf/spdf_stream.hh
//---------------------------------------------------------------------------//
