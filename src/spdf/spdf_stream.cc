//----------------------------------*-C++-*----------------------------------//
// spdf_stream.cc
// Geoffrey Furnish
// Thu Aug 10 23:26:33 1995
//---------------------------------------------------------------------------//
// @> C++ implementation of Maurice's SPDF.
//---------------------------------------------------------------------------//

#include <stdio.h>
#include <math.h>

#include <memory>

#include "ds++/Assert.hh"

#include "spdf/spdf_stream.hh"

// Scale factors for bit packing.  Will eventually be eliminated.

const U32 scale32bits = (U32) 0xFFFFFFFF;	/* 2^32 - 1 */
const U32 scale16bits = (U32) 0xFFFF;		/* 2^16 - 1 */

static const int NBUFSIZ = 512;

//---------------------------------------------------------------------------//
// Set initial state.
//---------------------------------------------------------------------------//

spdf_stream::spdf_stream()
{
    bp = 0;
    debug = 0;
}

/*--------------------------------------------------------------------------*\
 * pdf_wr_header()
 *
 * Writes a header string.  Input string must be NULL-terminated.  The
 * written string is terminated by a new-line, not a NULL.  This is done
 * so you can type e.g. "% strings <file> | head" and get sensible output.
\*--------------------------------------------------------------------------*/

int spdf_stream::wr_header( char *header )
{
    int i;

    dbug_enter("pdf_wr_header");

    for (i = 0; i < 79; i++) {
	if (header[i] == '\0')
	    break;
	if (spdf_putc(header[i]) == EOF)
	    return SPDF_WRERR;
    }
    if (spdf_putc('\n') == EOF)
	return SPDF_WRERR;

    return 0;
}

/*--------------------------------------------------------------------------*\
 * int pdf_rd_header
 *
 * Reads a newline-terminated header string from PDFstrm *pdfs, and
 * converts to a usual NULL-terminated string.  80 chars maximum assumed.
\*--------------------------------------------------------------------------*/

int spdf_stream::rd_header( char *header )
{
    int i;

    dbug_enter("pdf_rd_header");

    for (i = 0; i < 79; i++) {
	int c = spdf_getc();
	if (c == EOF)
	    return SPDF_RDERR;
	header[i] = c;

	if (header[i] == '\n')
	    break;
    }
    header[i] = '\0';		/* NULL terminate */
    return 0;
}

/*--------------------------------------------------------------------------*\
 * int pdf_wr_1byte()
 *
 * Writes a U_CHAR as a single byte.
\*--------------------------------------------------------------------------*/

void spdf_stream::wr_1byte( U_CHAR s )
{
    U_CHAR x[1];

    x[0] = s;

    wrx(x,1);
}

/*--------------------------------------------------------------------------*\
 * int pdf_rd_1byte()
 *
 * Reads a single byte, storing into a U_CHAR.
\*--------------------------------------------------------------------------*/

void spdf_stream::rd_1byte( U_CHAR& s )
{
    U_CHAR x[1];

    rdx(x,1);

    s = x[0];
}

/*--------------------------------------------------------------------------*\
 * pdf_wr_2bytes()
 *
 * Writes a U_SHORT as two single bytes, low end first.
\*--------------------------------------------------------------------------*/

void spdf_stream::wr_2bytes( U_SHORT s )
{
    U_CHAR x[2];

    x[0] = (U_CHAR) ((u32) (s & (u32) 0x00FF));
    x[1] = (U_CHAR) ((u32) (s & (u32) 0xFF00) >> 8);

    wrx(x,2);
}

/*--------------------------------------------------------------------------*\
 * pdf_rd_2bytes()
 *
 * Reads a U_SHORT from two single bytes, low end first.
\*--------------------------------------------------------------------------*/

void spdf_stream::rd_2bytes( U_SHORT& s )
{
    U_CHAR x[2];

    rdx(x,2);

    s = 0;
    s |= (U_SHORT) x[0];
    s |= (U_SHORT) x[1] << 8;
}

/*--------------------------------------------------------------------------*\
 * pdf_wr_2nbytes()
 *
 * Writes n U_SHORT's as 2n single bytes, low end first.
\*--------------------------------------------------------------------------*/

void spdf_stream::wr_2nbytes( U_SHORT *s, int n )
{
    int i;
    U_CHAR x[2];

    for (i = 0; i < n; i++) {
	x[0] = (U_CHAR) ((u32) (s[i] & (u32) 0x00FF));
	x[1] = (U_CHAR) ((u32) (s[i] & (u32) 0xFF00) >> 8);

	wrx(x,2);
    }
}

/*--------------------------------------------------------------------------*\
 * pdf_rd_2nbytes()
 *
 * Reads n U_SHORT's from 2n single bytes, low end first.
\*--------------------------------------------------------------------------*/

void spdf_stream::rd_2nbytes( U_SHORT *s, int n )
{
    int i;
    U_CHAR x[2];

    for (i = 0; i < n; i++) {
	rdx(x,2);

	s[i] = 0;
	s[i] |= (U_SHORT) x[0];
	s[i] |= (U_SHORT) x[1] << 8;
    }
}

/*--------------------------------------------------------------------------*\
 * pdf_wr_4bytes()
 *
 * Writes an unsigned long as four single bytes, low end first.
\*--------------------------------------------------------------------------*/

void spdf_stream::wr_4bytes( u32 s )
{
    U_CHAR x[4];

    x[0] = static_cast<U_CHAR> ((s & 0x000000FF));
    x[1] = static_cast<U_CHAR> ((s & 0x0000FF00) >> 8);
    x[2] = static_cast<U_CHAR> ((s & 0x00FF0000) >> 16);
    x[3] = static_cast<U_CHAR> ((s & 0xFF000000) >> 24);

    wrx(x,4);
}

/*--------------------------------------------------------------------------*\
 * pdf_rd_4bytes()
 *
 * Reads an unsigned long from 4 single bytes, low end first.
\*--------------------------------------------------------------------------*/

void spdf_stream::rd_4bytes( u32& s )
{
    U_CHAR x[4];

    rdx(x,4);

    s = 0;
    s |= static_cast<u32>( x[0] );
    s |= static_cast<u32>( x[1] ) << 8;
    s |= static_cast<u32>( x[2] ) << 16;
    s |= static_cast<u32>( x[3] ) << 24;
}

//---------------------------------------------------------------------------//
// Writes an unsigned long as eight single bytes, low end first.
//---------------------------------------------------------------------------//

void spdf_stream::wr_8bytes( u64 s )
{
    U_CHAR x[8];

    x[0] = static_cast<U_CHAR> ((s & 0x00000000000000FF));
    x[1] = static_cast<U_CHAR> ((s & 0x000000000000FF00) >> 8);
    x[2] = static_cast<U_CHAR> ((s & 0x0000000000FF0000) >> 16);
    x[3] = static_cast<U_CHAR> ((s & 0x00000000FF000000) >> 24);
    x[4] = static_cast<U_CHAR> ((s & 0x000000FF00000000) >> 32);
    x[5] = static_cast<U_CHAR> ((s & 0x0000FF0000000000) >> 40);
    x[6] = static_cast<U_CHAR> ((s & 0x00FF000000000000) >> 48);
    x[7] = static_cast<U_CHAR> ((s & 0xFF00000000000000) >> 56);

    wrx( x, 8 );
}

//---------------------------------------------------------------------------//
// Reads an unsigned long from 8 single bytes, low end first.
//---------------------------------------------------------------------------//

void spdf_stream::rd_8bytes( u64& s )
{
    U_CHAR x[8];

    rdx( x, 8 );

    s = 0;
    s |= static_cast<u64>( x[0] );
    s |= static_cast<u64>( x[1] ) << 8;
    s |= static_cast<u64>( x[2] ) << 16;
    s |= static_cast<u64>( x[3] ) << 24;
    s |= static_cast<u64>( x[4] ) << 32;
    s |= static_cast<u64>( x[5] ) << 40;
    s |= static_cast<u64>( x[6] ) << 48;
    s |= static_cast<u64>( x[7] ) << 56;
}

/*--------------------------------------------------------------------------*\
 * Here is the IEEE floating point specification in both 32 bit and 64 bit
 * precisions, from page 9 of "IEEE Standard for Binary Floating-Point
 * Arithmetic", copyright 1985, IEEE Std 754-1985:
 * 
 * 
 *                             Single Format
 * 
 * msb means most significant bit
 * lsb means least significant bit
 * 
 *   1         8                                23
 * _____________________________________________________________________
 * |   |                |                                              |
 * | s |       e        |                        f                     |
 * |___|________________|______________________________________________|
 *      msb          lsb msb                                        lsb
 * 
 * 
 * 
 *                             Double Format
 * 
 * msb means most significant bit
 * lsb means least significant bit
 * 
 *   1        11                                52
 * _____________________________________________________________________
 * |   |                |                                              |
 * | s |       e        |                        f                     |
 * |___|________________|______________________________________________|
 *      msb          lsb msb                                        lsb
 * 
 * 
 * (Thanks to: Andy Mai (mai@ncar.ucar.edu))
 * 
 * 
 * According to "inmos: Transputer instruction set" the IEEE standard
 * specifies the floating format as:
 * 
 *      s exp frac
 * 
 * Where: s = sign bit  (1 bit)
 *      exp = exponent (8 bits for 32 bit float / 11 bits for 64 bit float)
 *      frac = fraction (23 bits for 32 bit float / 52 bits for 64 bit float)
 * 
 * value of (s exp frac) = (-1)^s * 1.frac * 2^(exp-bias) ; if exp not 0
 *                         (-1)^s * 0.frac * 2^(1-bias) ; if exp = 0
 * 
 * where bias = 127 for 32 bit float
 *       bias = 1023 for 64 bit float
 * 
 * (Thanks to: Tom Bjorkholm(TBJORKHOLM@abo.fi))
 * 
\*--------------------------------------------------------------------------*/

//---------------------------------------------------------------------------//
// Writes a floating point numeric type in the relevant IEEE format.
//---------------------------------------------------------------------------//

template<class T>
void spdf_stream::wr_ieee( T  d )
{
// There are probably ANSI C++ ways (traits) to do most of this, but till I
// see what the actual standard provides, I'll just stuff traits for these
// operations directly into spdf_traits<T>.

    union {
	T dval;
	spdf_traits<T>::bitrep bval;
    } v;

#ifdef IEEE_BIG_ENDIAN
    v.dval = d;
#else

    if ( d == T(0) ) {
	v.bval = 0;
    }
    else {

	T d_tmp;
	int exp, e_off, bias = spdf_traits<T>::bias;
	spdf_traits<T>::bitrep s_ieee, e_ieee, d_ieee;

	double d_dbl = d;
	double dmant = frexp( d_dbl, &exp );

	if (dmant < 0)
	    s_ieee = 1;
	else
	    s_ieee = 0;

	dmant = fabs(dmant);
	double d_new = 2 * dmant;
	int e_new = exp - 1;

	if (e_new < 1 - bias) {
	    e_off = e_new - (1 - bias);
	    e_ieee = 0;
	    d_tmp = d_new * pow( 2.0, double(e_off) );
	}
	else {
	    e_ieee = e_new + bias;
	    d_tmp = d_new - 1;
	}

	d_ieee = d_tmp * spdf_traits<T>::shift_factor();

	if (e_ieee > spdf_traits<T>::max_exp) {
	    if (debug)
		fprintf(stderr, "pdf_wr_ieeef: Warning -- overflow\n");
	    e_ieee = spdf_traits<T>::max_exp;
	}

	s_ieee = s_ieee << (spdf_traits<T>::bits_mantissa +
			    spdf_traits<T>::bits_exponent);
	e_ieee = e_ieee << spdf_traits<T>::bits_mantissa;

	v.bval = s_ieee | e_ieee | d_ieee;
    }
#endif
    wr_ieee_bytes( v.bval );

    if (debug) {
    // Need to fix this stuff to be T aware.

	fprintf(stderr, "Double value (written):      %lg\n", d );
// 	print_ieeed( &d, &v.bval );
    }
}

//---------------------------------------------------------------------------//
// Read a floating point numeric type in the relevant IEEE format.
//---------------------------------------------------------------------------//

template<class T>
void spdf_stream::rd_ieee( T& d )
{
    union {
	T dval;
	spdf_traits<T>::bitrep bval;
    } v;

    rd_ieee_bytes( v.bval );

#ifdef IEEE_BIG_ENDIAN
    d = v.dval;
#else
    double d_new;
    int exp, bias = spdf_traits<T>::bias;

    spdf_traits<T>::bitrep s_ieee = v.bval >> (spdf_traits<T>::bits_mantissa +
					       spdf_traits<T>::bits_exponent);
    spdf_traits<T>::bitrep e_ieee =
	v.bval >> spdf_traits<T>::bits_mantissa & spdf_traits<T>::exp_mask;
    spdf_traits<T>::bitrep d_ieee = v.bval & spdf_traits<T>::mantissa_mask;
	
    double d_tmp = double(d_ieee) / spdf_traits<T>::shift_factor();

    if (e_ieee == 0) {
	exp = 1 - bias;
	d_new = d_tmp;
    }
    else {
	exp = int(e_ieee) - bias;
	d_new = 1.0 + d_tmp;
    }

    d = d_new * pow( 2.0, double(exp) );
    if (s_ieee == 1)
	d = -d;
#endif

    if (debug) {
    // Make this stuff T aware.
	fprintf(stderr, "Float value (read):      %lg\n", d);
// 	print_ieeed(&d, &v.bval);
    }
}

/*--------------------------------------------------------------------------*\
 * print_ieeef()
 *
 * Prints binary representation for numbers pointed to by arguments.
 * The first argument is the original float, the second is the
 * IEEE representation.  They should be the same on any machine that
 * uses IEEE floats.
\*--------------------------------------------------------------------------*/

void spdf_stream::print_ieeef( void *vx, void *vy )
{
    int i;
    u32 f, *x = (u32 *) vx, *y = (u32 *) vy;
    char bitrep[33];

    bitrep[32] = '\0';

    f = *x;
    for (i = 0; i < 32; i++) {
	if (f & 1)
	    bitrep[32 - i - 1] = '1';
	else
	    bitrep[32 - i - 1] = '0';
	f = f >> 1;
    }
    fprintf(stderr, "Binary representation:      ");
    fprintf(stderr, "%s\n", bitrep);

    f = *y;
    for (i = 0; i < 32; i++) {
	if (f & 1)
	    bitrep[32 - i - 1] = '1';
	else
	    bitrep[32 - i - 1] = '0';
	f = f >> 1;
    }
    fprintf(stderr, "Converted representation:   ");
    fprintf(stderr, "%s\n\n", bitrep);

    return;
}

//---------------------------------------------------------------------------//
// Same as above, but for doubles.
//---------------------------------------------------------------------------//

void spdf_stream::print_ieeed( void *vx, void *vy )
{
    int i;
    u32 f, *x = (u32 *) vx, *y = (u32 *) vy;
    char bitrep[65];

    bitrep[64] = '\0';

    f = *x;
    for (i = 0; i < 64; i++) {
	if (f & 1)
	    bitrep[64 - i - 1] = '1';
	else
	    bitrep[64 - i - 1] = '0';
	f = f >> 1;
    }
    fprintf(stderr, "Binary representation:      ");
    fprintf(stderr, "%s\n", bitrep);

    f = *y;
    for (i = 0; i < 64; i++) {
	if (f & 1)
	    bitrep[64 - i - 1] = '1';
	else
	    bitrep[64 - i - 1] = '0';
	f = f >> 1;
    }
    fprintf(stderr, "Converted representation:   ");
    fprintf(stderr, "%s\n\n", bitrep);

    return;
}

//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Get the byte id tag from the file, and watch for end of file.
//---------------------------------------------------------------------------//

U_CHAR spdf_stream::getbid()
{
    int bid = spdf_getc();
    if (bid == EOF)
	return 0;

    return (U_CHAR) bid;
}

int spdf_stream::eod()
{
// Write end of data byte id tag.

    spdf_putc( SPDF_TYPE_EOD );

    return 0;
}

//---------------------------------------------------------------------------//
// SPDF write, record of packed single precision floating point data.
//---------------------------------------------------------------------------//

template<class T>
void spdf_stream::wrp( T *v, int len, int nbits )
{
    u32 ubuf[NBUFSIZ], scale;
    U_CHAR pbuf[NBUFSIZ * 4];

// Get scale factor.

    switch (nbits) {
    case 16:
	scale = scale16bits;
	break;
    case 32:
	scale = scale32bits;
	break;
    default:
// 	fprintf(stderr,
// 		"%s: %d bit writes not currently supported\n",
// 		function, nbits);
// 	  return SPDF_BADNBITS;
	  throw "spdf_badnbits";
    }

// Write byte id tag.

    spdf_putc( spdf_traits<T>::packed_record_byte_id );

// Write header info.

    wr_4bytes(len);

    T vmin, vmax;
    minmax( v, len, vmin, vmax );
    wr_ieee( vmin );
    wr_ieee( vmax );

    wr_1byte((U_CHAR) nbits);

    if (vmin == vmax) return;	// constant data optimization.

// Write data.
    for( int ib=0; ib < len; ib += NBUFSIZ ) {
	int ilim = // std::min( NBUFSIZ, len - ib );
	    ( NBUFSIZ < (len-ib) ? NBUFSIZ : len-ib );

	int nbytes = nbits / 8 * ilim;
	for( int i=0; i < ilim; i++ ) {
	    int j = ib + i;
	    ubuf[i] = scale * ( (v[j] - vmin) / (vmax - vmin) );
	}
	packw( pbuf, nbits, ubuf, ilim );
	wrx( pbuf, nbytes );
    }
}

/*--------------------------------------------------------------------------*\
* int spdf_rdfrp()
*
* SPDF read, record of floating point packed data.
* Returns the actual number of floats read in *p_len.
* If the number of floats in the record exceeds lmax, the function aborts.
\*--------------------------------------------------------------------------*/

template<class T>
void spdf_stream::rdp( T *v, int& len, int max )
{
    u32 ubuf[NBUFSIZ], scale;
    U_CHAR pbuf[NBUFSIZ * 4];

// Read & check byte id tag.

    match_byte_id( spdf_traits<T>::packed_record_byte_id );

// Read header info.

    u32 ul_len;
    rd_4bytes( ul_len );

    len = static_cast<int>( ul_len );
    if (len > max)
	throw "SPDF_BADLEN";

    T vmin, vmax;
    rd_ieee( vmin );
    rd_ieee( vmax );

    U_CHAR s;
    rd_1byte( s );
    int nbits = s;

    if (vmin == vmax) {
	for( int i=0; i < len; i++) {
	    v[i] = vmin;
	}
	return;
    }

// Get scale factor.

    switch (nbits) {
    case 16:
	scale = scale16bits;
	break;
    case 32:
	scale = scale32bits;
	break;
    default:
	throw "SPDF_BADNBITS";
    }

// Read data
// NOT RIGHT IF NBITS NOT EVENLY DIVISIBLE BY 8!

    for( int ib=0; ib < len; ib += NBUFSIZ ) {
	int ilim = // std::min( NBUFSIZ, len - ib );
	    ( NBUFSIZ < (len-ib) ? NBUFSIZ : len-ib );

	int nbytes = nbits / 8 * ilim;
	rdx( pbuf, nbytes );
	unpackw(pbuf, nbits, ubuf, ilim);
	for( int i=0; i < ilim; i++ ) {
	    int j = ib + i;
	    v[j] = vmin + ((ubuf[i] + 0.5)/ T(scale)) * (vmax - vmin);
	}
    }
}

//---------------------------------------------------------------------------//
// Copies stream buffer of longs to a stream of bytes, keeping only 'nbits'
// bits of precision (the lower bits are dropped off).
//
// This routine assumes there is at maximum 32 bits of information in each
// u[i] (so if a long is 64 bits, the upper 32 won't be used here).
// This is pretty easy right now since the number of bits is always a
// multiple of the byte size (8).  Implementing 12 or 24 bit packing will
// require a bit more work.
//---------------------------------------------------------------------------//

void spdf_stream::packw( U_CHAR *p, int nbits, u32 *u, int ilim  )
{
    int i, m = 0;

    if (nbits == 16) {

	for (i = 0; i < ilim; i++) {

	    p[m++] = (U_CHAR) ((u[i] & (u32) 0x000000FF));
	    p[m++] = (U_CHAR) ((u[i] & (u32) 0x0000FF00) >> 8);
	}
    }
    else if (nbits == 32) {

	for (i = 0; i < ilim; i++) {

	    p[m++] = (U_CHAR) ((u[i] & (u32) 0x000000FF));
	    p[m++] = (U_CHAR) ((u[i] & (u32) 0x0000FF00) >> 8);
	    p[m++] = (U_CHAR) ((u[i] & (u32) 0x00FF0000) >> 16);
	    p[m++] = (U_CHAR) ((u[i] & (u32) 0xFF000000) >> 24);
	}
    }
    return;
}

//---------------------------------------------------------------------------//
// Copies stream buffer of bytes to a stream of longs, storing only 'nbits'
// bits of precision.  Only supports nbits = 16 or nbits = 32 at present.
//---------------------------------------------------------------------------//

void spdf_stream::unpackw( U_CHAR *p, int nbits, u32 *u, int ilim )
{
    int i, m = 0;

    if (nbits == 16) {

	for( i = 0; i < ilim; i++ ) {
	    u[i] = 0;
	    u[i] |= (u32) p[m++];
	    u[i] |= (u32) p[m++] << 8;
	}
    }
    else if (nbits == 32) {

	for( i = 0; i < ilim; i++ ) {
	    u[i] = 0;
	    u[i] |= (u32) p[m++];
	    u[i] |= (u32) p[m++] << 8;
	    u[i] |= (u32) p[m++] << 16;
	    u[i] |= (u32) p[m++] << 24;
	}
    }
    return;
}

//---------------------------------------------------------------------------//
// Minimum and maximum elements of a 1-d data.
//---------------------------------------------------------------------------//
/*
void spdf_stream::fminmx( float *pf, int len, float& fmin, float& fmax )
{
    fmin = fmax = pf[0];

    for( int i=1; i < len; i++ ) {
	if (pf[i] < fmin) fmin = pf[i];
	if (pf[i] > fmax) fmax = pf[i];
    }
}

//---------------------------------------------------------------------------//
// Minimum and maximum elements of a 1-d data.  Double precision.
//---------------------------------------------------------------------------//

void spdf_stream::fminmx( double *pf, int len, double& fmin, double& fmax )
{
    fmin = fmax = pf[0];

    for( int i=1; i < len; i++ ) {
	if (pf[i] < fmin) fmin = pf[i];
	if (pf[i] > fmax) fmax = pf[i];
    }
}
*/
//---------------------------------------------------------------------------//
// Derived classes.  stdio_spdf_stream implements an spdf_stream for
// interfaceing with the file system through stdio.  mbuf_spdf_stream uses a
// memory buffer as its "storage medium".
//---------------------------------------------------------------------------//

stdio_spdf_stream::stdio_spdf_stream( const char *n, const char *a )
{
    fp = fopen( n, a );

    if (fp == NULL) {
	throw( "Unable to open file." );
    }
}

stdio_spdf_stream::~stdio_spdf_stream()
{
    fclose(fp);
}

int stdio_spdf_stream::spdf_putc( int c )
{
    int r = putc( c, fp );
    bp++;

    return r;
}

int stdio_spdf_stream::spdf_getc()
{
    int r = getc( fp );
    bp++;

    return r;
}

int stdio_spdf_stream::spdf_ungetc( int c )
{
    int r = ungetc( c, fp );
    if (bp) bp--;

    return r;
}

//---------------------------------------------------------------------------//
// Writes a record.
//---------------------------------------------------------------------------//

void stdio_spdf_stream::wrx( const U_CHAR *x, int nitems )
{
    int result = fwrite( x, 1, nitems, fp );

    if (result != nitems)
	throw( "stdio_spdf_stream::wrx, unable to write data" );

    bp += nitems;
}

//---------------------------------------------------------------------------//
// Reads a record.
//---------------------------------------------------------------------------//

void stdio_spdf_stream::rdx( U_CHAR *x, int nitems )
{
    int result = fread( x, 1, nitems, fp );

    if (result != nitems)
	throw( "stdio_spdf_stream::rdx, unable to read data" );

    bp += nitems;
}

//---------------------------------------------------------------------------//

mbuf_spdf_stream::mbuf_spdf_stream( int _bufmax /*=2048*/ )
    : bufmax(_bufmax)
{
    v = new U_CHAR[ bufmax ];
    may_delete = 1;
}

mbuf_spdf_stream::mbuf_spdf_stream( U_CHAR *_v, int _bufmax )
    : v(_v), bufmax(_bufmax)
{
    may_delete = 0;
}

mbuf_spdf_stream::~mbuf_spdf_stream()
{
    if (may_delete) delete[] v;
}

//---------------------------------------------------------------------------//
// Ensure that the memory buffer is at least a certain size.  This is so that
// we can do a direct binary write into the buffer without fearing overflow.
//---------------------------------------------------------------------------//

void mbuf_spdf_stream::set_min_size( int s )
{
    if (s <= bufmax )
	return;

    Insist( may_delete && bp == 0,
	    "Can't expand buffer in mbuf_spdf_stream::set_min_size." );

    delete[] v;
    bufmax = s;
    v = new U_CHAR[ bufmax ];
}

int mbuf_spdf_stream::spdf_putc( int c )
{
    Assert( v );

    if (bp >= bufmax) {
	Insist( may_delete,
		"Can't expand buffer in mbuf_spdf_stream::spdf_putc." );

	int nbufmax = bufmax + 512;
	U_CHAR *nv = new U_CHAR[ nbufmax ];

	std::copy( v, v+bp, nv );

	delete[] v;

	bufmax = nbufmax;
	v = nv;
    }

    v[bp++] = c;

    return c;
}

int mbuf_spdf_stream::spdf_getc()
{
    int r = EOF;

    if (bp < bufmax)
	return v[bp++];

    return r;
}

int mbuf_spdf_stream::spdf_ungetc( int c )
{
    if (bp)
	v[--bp] = c;

    return c;
}

//---------------------------------------------------------------------------//
// Writes a record.
//---------------------------------------------------------------------------//

void mbuf_spdf_stream::wrx( const U_CHAR *x, int nitems )
{
    for( int i = 0; i < nitems; i++) {
	if ( bp >= bufmax ) {
	    if (may_delete) {
		int nbufmax = bufmax + 512;
		U_CHAR *nv = new U_CHAR[ nbufmax ];

		std::copy( v, v+bp, nv );

		delete[] v;

		bufmax = nbufmax;
		v = nv;
	    } else {
		throw( "Out of space, but not allowed to realloc!" );
	    }
	}
	v[bp++] = x[i];
    }
}

//---------------------------------------------------------------------------//
// Reads a record.
//---------------------------------------------------------------------------//

void mbuf_spdf_stream::rdx( U_CHAR *x, int nitems )
{
    int i, result = 0;

    for (i = 0; i < nitems; i++) {
	if (bp > bufmax)
	    break;
	x[i] = v[bp++];
    }
    result = i;

    if (result != nitems)
	throw( "mbuf_spdf_stream::rdx, not enough data." );
}

//---------------------------------------------------------------------------//
// Putchar.  Illegal for this class.
//---------------------------------------------------------------------------//

int spdf_ifstream::spdf_putc( int c )
{
    throw( "Cannot write to an spdf_ifstream." );
}

//---------------------------------------------------------------------------//
// Fetch a single char from the managed ifstream.
//---------------------------------------------------------------------------//

int spdf_ifstream::spdf_getc()
{
//     char r = EOF;

//     if (!ifs.eof())
// 	ifs.get(r);

//     return r;

// Try it this way for a while, see if it works out okay...  (gmf, 8/13/96).

    if (ifs.eof())
	return EOF;

    char r;
    ifs.get(r);
    return r;
}

//---------------------------------------------------------------------------//
// PUsh a character back onto the stream.
//---------------------------------------------------------------------------//

int spdf_ifstream::spdf_ungetc( int c )
{
    char ch = c;
    ifs.putback(ch);

    return 1;
}

//---------------------------------------------------------------------------//
// Output to an ifstream is a no no.
//---------------------------------------------------------------------------//

void spdf_ifstream::wrx( const U_CHAR *x, int nitems )
{
    throw( "Cannot write to an spdf_ifstream." );
}

//---------------------------------------------------------------------------//
// Binary read from the managed ifstream.
//---------------------------------------------------------------------------//

void spdf_ifstream::rdx( U_CHAR *x, int nitems )
{
    ifs.read( (char *) x, nitems );
}


//---------------------------------------------------------------------------//
// Put a character out to the ofstream.
//---------------------------------------------------------------------------//

int spdf_ofstream::spdf_putc( int c )
{
    char ch = c;

    ofs << ch;

    return 1;
}

//---------------------------------------------------------------------------//
// Cannot fetch from an ofstream.
//---------------------------------------------------------------------------//

int spdf_ofstream::spdf_getc()
{
    throw( "Cannot fetch from an ofstream." );
}

//---------------------------------------------------------------------------//
// Also illegal.
//---------------------------------------------------------------------------//

int spdf_ofstream::spdf_ungetc( int c )
{
    throw( "Should not be ungetting, b/c cannot get." );
}

//---------------------------------------------------------------------------//
// Binary output to the ofstream.
//---------------------------------------------------------------------------//

void spdf_ofstream::wrx( const U_CHAR *x, int nitems )
{
    ofs.write( (const char *) x, nitems );
}

//---------------------------------------------------------------------------//
// Also illegal.
//---------------------------------------------------------------------------//

void spdf_ofstream::rdx( U_CHAR *x, int nitems )
{
    throw( "CAnnot read from an ofstream." );
}

//---------------------------------------------------------------------------//
//                              end of spdf_stream.cc
//---------------------------------------------------------------------------//
