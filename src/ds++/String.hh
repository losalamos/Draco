//----------------------------------*-C++-*----------------------------------//
// String.hh
// Geoffrey Furnish
// 28 January 1994
//---------------------------------------------------------------------------//
// @> A basic String class.
//---------------------------------------------------------------------------//

#ifndef __ds_String_hh__
#define __ds_String_hh__

#include <iosfwd>

#include <string.h>

class String {

    struct srep {
	char *s;
	int n;			// reference count
	int len;		// strlen(s)
	int sz;			// size of s

	srep( int _sz ) : sz(_sz), n(1), len(0) {
	    s = new char[ sz ];
	    s[0] = '\0';
	}

	~srep() { delete[] s; }
    };
    srep *p;

    char bogus;			// return this instead of throw(range);

  public:
    String( int startlen =1 );
    String( const String& s );
    String( const char *s );
    String( const char c );

    ~String();

    String& operator=( const char c );
    String& operator=( const char *s );
    String& operator=( const String& s );

    operator const char *() const { return p->s; }
    int len() const { return p->len; }
    int length() const { return p->len; } // SC/libg++ compatibility.
    char& operator[]( int i );
    char  operator[]( int i ) const;

  friend std::ostream& operator<<( std::ostream& s, const String& x );
  friend std::istream& operator>>( std::istream& s, String& x );

  friend inline int operator==( const String& x, const String& y )
    { return strcmp( x.p->s, y.p->s ) == 0; }
  friend inline int operator==( const String& x, const char *s )
    { return strcmp( x.p->s, s ) == 0; }
    
  friend inline int operator!=( const String& x, const String& y )
    { return strcmp( x.p->s, y.p->s ) != 0; }
  friend inline int operator!=( const String& x, const char *s )
    { return strcmp( x.p->s, s ) != 0; }
    
// These relational operators work fine on ASCII machines, but if you
// have EBCDIC or some such, these probably don't do what you want.
    
  friend inline int operator<( const String& x, const String& y )
    { return strcmp( x.p->s, y.p->s ) < 0; }
  friend inline int operator<( const String& x, const char *s )
    { return strcmp( x.p->s, s ) < 0; }
    
  friend inline int operator<=( const String& x, const String& y )
    { return strcmp( x.p->s, y.p->s ) <= 0; }
  friend inline int operator<=( const String& x, const char *s )
    { return strcmp( x.p->s, s ) <= 0; }
    
  friend inline int operator>( const String& x, const String& y )
    { return strcmp( x.p->s, y.p->s ) > 0; }
  friend inline int operator>( const String& x, const char *s )
    { return strcmp( x.p->s, s ) > 0; }
    
  friend inline int operator>=( const String& x, const String& y )
    { return strcmp( x.p->s, y.p->s ) >= 0; }
  friend inline int operator>=( const String& x, const char *s )
    { return strcmp( x.p->s, s ) >= 0; }
    
// Now some operators for doing concatenation and such.
    
    String& operator+=( const String& s );
    String& operator+=( const char *s );
    String& operator+=( const char c );

  friend String operator+( const String& x, const String& y );
  friend String operator+( const String& x, const char *s );
  friend String operator+( const char *s, const String& x );

  friend String operator+( const String& s, const char c );
  friend String operator+( const char c, const String& s );

};

#endif				// __ds_String_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/String.hh
//---------------------------------------------------------------------------//
