//----------------------------------*-C++-*----------------------------------//
// persist.hh
// Geoffrey Furnish
// 16 May 1994
//---------------------------------------------------------------------------//
// @> A persistence class to help work with the Perstalyzer.
//---------------------------------------------------------------------------//

#ifndef __util_persist_hh__
#define __util_persist_hh__

#include <iostream>
#include <fstream>

#include "ds++/String.hh"
#include "ds++/Mat.hh"
#include "ds++/DynArray.hh"

template<class T>
class prst_read_xcpt {
  public:
    std::ifstream& ifs;
    const char *name;
    T& ref;

    prst_read_xcpt( std::ifstream& _ifs, const char *_name, T& _ref )
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

    int p_printf( std::ofstream& ifs, const char *s, ... );

// Unfortunately, we can't have a p_scanf() because too many systems
// don't have a vsscanf(), and I can't find a truly portable version.

    void p_require( std::ifstream& ifs, const char *s );

    void prst_write( std::ofstream& ofs, const char *s, int  v );
    void prst_read(  std::ifstream& ifs, const char *s, int& v );

    void prst_write( std::ofstream& ofs, const char *s, dsxx::String  v );
    void prst_read(  std::ifstream& ifs, const char *s, dsxx::String& v );

    void prst_write( std::ofstream& ofs, const char *s,
		     dsxx::DynArray<int>& v );
    void prst_read(  std::ifstream& ifs, const char *s,
		     dsxx::DynArray<int>& v );

    void prst_write( std::ofstream& ofs, const char *s,
		     dsxx::DynArray<float>& v );
    void prst_read(  std::ifstream& ifs, const char *s,
		     dsxx::DynArray<float>& v );

    void prst_write( std::ofstream& ofs, const char *s,
		     dsxx::DynArray<double>& v );
    void prst_read(  std::ifstream& ifs, const char *s,
		     dsxx::DynArray<double>& v );

    void prst_write( std::ofstream& ofs, const char *s, float  v );
    void prst_read(  std::ifstream& ifs, const char *s, float& v );

    void prst_write( std::ofstream& ofs, const char *s, double  v );
    void prst_read(  std::ifstream& ifs, const char *s, double& v );

    void prst_write( std::ofstream& ofs, const char *s, dsxx::Mat1<int>& v );
    void prst_read(  std::ifstream& ifs, const char *s, dsxx::Mat1<int>& v );

    void prst_write( std::ofstream& ofs, const char *s, dsxx::Mat1<float>& v );
    void prst_read(  std::ifstream& ifs, const char *s, dsxx::Mat1<float>& v );

    void prst_write( std::ofstream& ofs, const char *s, dsxx::Mat1<double>& v );
    void prst_read(  std::ifstream& ifs, const char *s, dsxx::Mat1<double>& v );
};

#endif                          // __util_persist_hh__

//---------------------------------------------------------------------------//
//                              end of util/persist.hh
//---------------------------------------------------------------------------//
