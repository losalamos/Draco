//----------------------------------*-C++-*----------------------------------//
// config.hh
// Geoffrey Furnish
// Mon Nov 25 17:07:59 1996
//---------------------------------------------------------------------------//
// @> Configure SPDF for machine architecture.
//---------------------------------------------------------------------------//

#ifndef __spdf_config_hh__
#define __spdf_config_hh__

#if defined(__hpux)
#define IEEE_BIG_ENDIAN
#endif

// Define architecture specializations here.

// DEC Alphaservers running OSF, are wacky beyond belief.  float and double
// appear to have bits ordered in the opposite direction, etc.

#if defined(__alpha) && defined(__osf__)

#define U32 unsigned int
#define U64 unsigned long

#endif

#if defined(__sgi)
#define U32 unsigned int
#endif

// Define default values for some unsigned types here.  Use ifndef since they
// may have been specialized above, for particular machine architectures.

#ifndef U_CHAR
#define U_CHAR unsigned char
#endif

#ifndef U_SHORT
#define U_SHORT unsigned short
#endif

#ifndef U_INT
#define U_INT unsigned int
#endif

#ifndef U32
#define U32 unsigned long
#endif

#ifndef U64
#define U64 unsigned long long
#endif


//---------------------------------------------------------------------------//
// Record types.  Each write starts out with a 1-byte type flag.

enum {
    SPDF_TYPE_EOD	= 0x10,	// End of data
    SPDF_TYPE_S	= 0x20,		// String
    SPDF_TYPE_I	= 0x30,		// Integer
    SPDF_TYPE_IR	= 0x31,	// Integer record
    SPDF_TYPE_F	= 0x40,		// Float
    SPDF_TYPE_FR	= 0x41,	// Float record
    SPDF_TYPE_FRP	= 0x42,	// Float record, packed
    SPDF_TYPE_D	= 0x50,		// Double
    SPDF_TYPE_DR	= 0x51,	// Double record
    SPDF_TYPE_DRP	= 0x52	// Double record, packed
};

//---------------------------------------------------------------------------//
// Now lets state what the various types we will be using will be called.

class spdf_types {
  public:
    typedef unsigned char uchar;
    typedef U32 u32;
    typedef U64 u64;
};

template<class T> class spdf_traits {};

template<>
class spdf_traits<int> {
  public:
    enum { size = 4 };
    static const int scalar_byte_id = SPDF_TYPE_I;
    static const int record_byte_id = SPDF_TYPE_IR;
};

template<>
class spdf_traits<float> {
  public:
    enum { size = 4 };
    static const char *name() { return "float"; }
    static const int scalar_byte_id = SPDF_TYPE_F;
    static const int record_byte_id = SPDF_TYPE_FR;
    static const int packed_record_byte_id = SPDF_TYPE_FRP;
    typedef spdf_types::u32 bitrep;
    static const int bias = 127;
    static float shift_factor() { return 8388608; }// 2^23
    static const int max_exp = 255;
    static const int bits_mantissa = 23;
    static const int bits_exponent = 8;
    static const bitrep exp_mask = 0xFF;
    static const bitrep mantissa_mask = 0x007FFFFF;
    static const int compress_50 = 16;
};

template<>
class spdf_traits<double> {
  public:
    enum { size = 8 };
    static const char *name() { return "double"; }
    static const int scalar_byte_id = SPDF_TYPE_D;
    static const int record_byte_id = SPDF_TYPE_DR;
    static const int packed_record_byte_id = SPDF_TYPE_DRP;
    typedef spdf_types::u64 bitrep;
    static const int bias = 1023;
    static double shift_factor() { return 8388608. * 8388608. * 64.; } // 2^52
    static const int max_exp = 2047;
    static const int bits_mantissa = 52;
    static const int bits_exponent = 11;
    static const bitrep exp_mask = 0x7FF;
    static const bitrep mantissa_mask = 0x0fffffffffffffL;
    static const int compress_50 = 32;
};

#endif                          // __spdf_config_hh__

//---------------------------------------------------------------------------//
//                              end of spdf/config.hh
//---------------------------------------------------------------------------//
