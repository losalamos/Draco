//----------------------------------*-C++-*----------------------------------//
// persist.hh
// Geoffrey Furnish
// 16 May 1994
//---------------------------------------------------------------------------//
// @> A persistence class to help work with the Perstalyzer.
//---------------------------------------------------------------------------//

#ifndef __util_persist_hh__
#define __util_persist_hh__

#include <iostream.h>
#include <fstream.h>

#include "ds++/String.hh"
#include "ds++/Mat.hh"
#include "ds++/DynArray.hh"

using dsxx::Mat1;

template<class T>
class prst_read_xcpt {
  public:
    ifstream& ifs;
    const char *name;
    T& ref;

    prst_read_xcpt( ifstream& _ifs, const char *_name, T& _ref )
	: ifs(_ifs), name(_name), ref(_ref) {}
};

class persist {

    char *buf;
    char *ss;
    enum { buflen = 200 };

  public:
    persist();
    persist( const persist& p );
    ~persist();
    persist& operator=(const persist& p );

    int p_printf( ofstream& ifs, const char *s, ... );

// Unfortunately, we can't have a p_scanf() because too many systems
// don't have a vsscanf(), and I can't find a truly portable version.

    void p_require( ifstream& ifs, const char *s );

    void prst_write( ofstream& ofs, const char *s, int  v );
    void prst_read(  ifstream& ifs, const char *s, int& v );

    void prst_write( ofstream& ofs, const char *s, String  v );
    void prst_read(  ifstream& ifs, const char *s, String& v );

    void prst_write( ofstream& ofs, const char *s, DynArray<int>& v );
    void prst_read(  ifstream& ifs, const char *s, DynArray<int>& v );

    void prst_write( ofstream& ofs, const char *s, DynArray<float>& v );
    void prst_read(  ifstream& ifs, const char *s, DynArray<float>& v );

    void prst_write( ofstream& ofs, const char *s, DynArray<double>& v );
    void prst_read(  ifstream& ifs, const char *s, DynArray<double>& v );

    void prst_write( ofstream& ofs, const char *s, float  v );
    void prst_read(  ifstream& ifs, const char *s, float& v );

    void prst_write( ofstream& ofs, const char *s, double  v );
    void prst_read(  ifstream& ifs, const char *s, double& v );

    void prst_write( ofstream& ofs, const char *s, Mat1<int>& v );
    void prst_read(  ifstream& ifs, const char *s, Mat1<int>& v );

    void prst_write( ofstream& ofs, const char *s, Mat1<float>& v );
    void prst_read(  ifstream& ifs, const char *s, Mat1<float>& v );

    void prst_write( ofstream& ofs, const char *s, Mat1<double>& v );
    void prst_read(  ifstream& ifs, const char *s, Mat1<double>& v );
};

#endif                          // __util_persist_hh__

//---------------------------------------------------------------------------//
//                              end of util/persist.hh
//---------------------------------------------------------------------------//
