//----------------------------------*-C++-*----------------------------------//
// String.cc
// Geoffrey Furnish
// 28 Januar 1994
//---------------------------------------------------------------------------//
// @> A basic String class.
//---------------------------------------------------------------------------//

#include <iostream>

#include "String.hh"
#include <string.h>

#include "Assert.hh"

NAMESPACE_DS_BEG

//===========================================================================//
// class String - char * made useful

// Everybody needs a String class.  The basic purpose is to put an end to all
// the tedium of working with char *'s.  
//
// This String class is loosely based on Stroustrup's String class presented
// on pp 248 in the 2nd edition.  Originally it was going to be almost
// exactly his class, but as I was typing it in, I realized it would be
// possible to do significantly better in terms of runtime performance, for
// only a small amount of additional effort.  Not intended as criticism of
// his class, since he always stresses language features over algorithms.
//
// The proscription for a performance oriented String class stems from the
// age-old tradeoff of size versus speed.  Figuring that memory space is not
// such a problem nowadays as it was back in the early 80's, this String
// class saves some extra state so that operations like strlen() aren't
// required more than once.  Furthermore, and undoubdtedly more important,
// calls to new and delete are limited as much as possible.  Storage is only
// freed when the String is actually being destroyed.  The char * buffer used
// by this String class is then of monotonically increasing size, until such
// time as the String is destructed.
//
// A variety of convenience functions are provided for working with Strings,
// but the class is not cluttered with things which might negatively impact
// portability.  For example, no support is provided for Regexp's, a common
// part of String class culture, since the author uses a variety of embedded
// systems for which Unix libc is not available.
//===========================================================================//

//---------------------------------------------------------------------------//
// Instantiate a String with a specified initial size.  This is not a hard
// limit though, it will still grow on demand.
//---------------------------------------------------------------------------//

String::String( int startlen /*=1*/ )
{
    p = new srep(startlen);
    p->s[0] = '\0';
}

//---------------------------------------------------------------------------//
// Instantiate a String from another.  This class uses reference counting for
// efficiency. 
//---------------------------------------------------------------------------//

String::String( const String& s )
{
    s.p->n++;
    p = s.p;
}

//---------------------------------------------------------------------------//
// Construct a String from a char *.
//---------------------------------------------------------------------------//

String::String( const char *s )
{
    int len = strlen(s);
    p = new srep(len+1);
    p->len = len;
    strcpy( p->s, s );
}

//---------------------------------------------------------------------------//
// Construct a String from a char.
//---------------------------------------------------------------------------//

String::String( const char c )
{
    p = new srep(5);
    p->len = 1;
    p->s[0] = c;
    p->s[1] = '\0';
}

//---------------------------------------------------------------------------//
// Detach ourselves from the managed String data.  If we were the last one
// using it, free it too.
//---------------------------------------------------------------------------//

String::~String()
{
    if ( --p->n == 0 )
      delete p;
}

//---------------------------------------------------------------------------//
// Safe assignment.  Use reference counting and copy on write semantics.
//---------------------------------------------------------------------------//

String& String::operator=( const char c )
{
    int newsz  = 2;

    if ( p->n > 1 ) {		// disconnect self
	p->n--;
	p = new srep(newsz+3);
    } else			// No one else is using this srep, so
				// let's make sure its big enough.
      if ( p->sz < newsz ) {
	  delete p;
	  p = new srep(newsz+3);
      }

    p->len = 1;
    p->s[0] = c;
    p->s[1] = '\0';
    return *this;
}

//---------------------------------------------------------------------------//
// Safe assignment.  Use reference counting and copy on write semantics.
//---------------------------------------------------------------------------//

String& String::operator=( const char *s )
{
    int newlen = strlen(s);
    int newsz  = newlen +1;

    if ( p->n > 1 ) {		// disconnect self
	p->n--;
	p = new srep(newsz);
    } else			// No one else is using this srep, so
				// let's make sure its big enough.
      if ( p->sz < newsz ) {
	  delete p;
	  p = new srep(newsz);
      }

    p->len = newlen;
    strcpy( p->s, s );
    return *this;
}

//---------------------------------------------------------------------------//
// Alias to another String.  Exploit reference counting, and copy on write as
// needed. 
//---------------------------------------------------------------------------//

String& String::operator=( const String& s )
{
    s.p->n++;			// protect against st = st;
    if ( --p->n == 0 )
      delete p;

    p = s.p;
    return *this;
}

// Safe indexing.

// NOTE: for these two [] funcs, we allow them to access with index == len,
// in order to allow detection of the terminating '\0'.

//---------------------------------------------------------------------------//
// Export a reference to a char in a String.  Note that since this reference
// is writeable, in order to preserve the integrity of any const String's
// which might be aliased to this one, we have to regard this as a "write"
// operation, and detach from any aliased entities. 
//---------------------------------------------------------------------------//

char& String::operator[]( int i )
{
    Insist( i >= 0 && i <= p->len, "Subscript out of range." );

    if (p->n > 1) {		// Then we must disconnect.
	p->n--;
	srep *np = new srep(p->len + 1);
	strcpy( np->s, p->s );
	np->len = p->len;
	p = np;
    }

    return p->s[i];
}

//---------------------------------------------------------------------------//
// Export the requested char from the const String.  This one does not expose
// the internal state of the String object.
//---------------------------------------------------------------------------//

char String::operator[]( int i ) const
{
    Insist( i >= 0 && i <= p->len, "Subscript out of range." );

    return p->s[i];
}

//---------------------------------------------------------------------------//
// Output to i/o streams.
//---------------------------------------------------------------------------//

std::ostream& operator<<( std::ostream& s, const String& x )
{
    return s << x.p->s;
}

//---------------------------------------------------------------------------//
// Input a String.  Have to look at this sometime to make sure we avoid buffer
// overflow. 
//---------------------------------------------------------------------------//

std::istream& operator>>( std::istream& s, String& x )
{
    char buf[256];

    s >> buf;			// Could overflow.

    x = buf;

    return s;
}

// Comparison functions.

//---------------------------------------------------------------------------//
// Concatenation functions.

// Okay, the logic is:
// If my reference count is > 1, then decrement it, get myself a new srep,
// and make sure it's big enough for the new String.  Be careful to get the
// copying done before losing the old srep.

// If reference count == 1, then this is a solo.  Check to see if
// concatenating to self.  If so, 
// ...
//---------------------------------------------------------------------------//

// This whole function should really be fixed to avoid use of strcat, since
// strcpy ought to be faster since we already know the end of hte existing
// string is p->s + p-> len.  Will try to make this optimization later.

//---------------------------------------------------------------------------//
// Concatenate a String to this one.
//---------------------------------------------------------------------------//

String& String::operator+=( const String& s )
{
    int slen = s.len();

    if ( p->n > 1 ) {		// Shared on entry.
	p->n--;			// Disconnect self.

// Get new srep, and initialize with right values

	srep *ps = new srep( p->len + slen + 1 );
	strcpy( ps->s, p->s );
	strcpy( ps->s + p->len, s.p->s ); // faster than strcat.
	ps->len = p->len + slen;

// Set new srep and go home.

	p = ps;
	return *this;
    }

// Not sharing with anyone.

    if (this == &s) {		// Watch out for st += st;

// Working with self, so don't adjust reference count.

	if ( p->sz > 2*p->len ) { // There's already enough space.

	    strncpy( p->s + p->len, p->s, p->len );
	    p->len *= 2;
	    p->s[p->len] = '\0';
	    return *this;

	} else {		// Need to allocate more space, carefully!

	    srep *ps = new srep( 2*p->len + 1 );
	    ps->len = 2*p->len;
	    strcpy( ps->s, p->s );
	    strcat( ps->s, p->s );
	    delete p;
	    p = ps;
	    return *this;
	}
    }

// Okay, they're different strings, good.

// Have to be sure that p->len + slen + 1 ('\0') fit into p->sz.

    if ( p->len + slen >= p->sz ) { // Then we need more space.

	srep *ps = new srep( p->len + slen + 1 );
	ps->len = p->len + slen;
	strcpy( ps->s, p->s );
	strcat( ps->s, s.p->s );
	delete p;		// since no one will be using it.
	p = ps; 
	return *this;

    } else {			// We already have enough space.

	p->len += slen;
	strcat( p->s, s.p->s );
	return *this;
    }
}

//---------------------------------------------------------------------------//
// Concatenate a char * to this String.
//---------------------------------------------------------------------------//

String& String::operator+=( const char *s )
{
    int slen = strlen(s);

    if ( p->n > 1 ) {		// Shared on entry.
	p->n--;			// Disconnect self.

	srep *ps = new srep( p->len + slen + 1 );
	ps->len = p->len + slen;
	strcpy( ps->s, p->s );
	strcpy( ps->s + p->len, s );	// faster than strcat.
	p = ps;
	return *this;
    }

    if ( p->len + slen >= p->sz ) { // then we need more space.
	srep *ps = new srep( p->len + slen + 1 );
	ps->len = p->len + slen;
	strcpy( ps->s, p->s );
	strcpy( ps->s + p->len, s ); // faster than strcat.
	delete p;
	p = ps;
	return *this;

    } else {

	strcpy( p->s + p->len, s );
	p->len += slen;
	return *this;
    }
}
    
//---------------------------------------------------------------------------//
// Postpend a character.
//---------------------------------------------------------------------------//

String& String::operator+=( const char c )
{
    if ( p->n > 1 ) {		// Shared on entry.
	p->n--;			// Disconnect self.

	srep *ps = new srep( p->len + 5 + 1 );
	ps->len = p->len + 1;
	strcpy( ps->s, p->s );
	ps->s[p->len] = c;
	ps->s[p->len + 1] = '\0';
	p = ps;
	return *this;
    }

    if ( p->len + 1 >= p->sz ) { // then we need more space.
	srep *ps = new srep( p->len + 5 + 1 );
	ps->len = p->len + 1;
	strcpy( ps->s, p->s );
	ps->s[p->len] = c;
	ps->s[p->len + 1] = '\0';
	delete p;
	p = ps;
	return *this;

    } else {

	p->s[p->len] = c;
	p->s[p->len + 1] = '\0';
	p->len += 1;
	return *this;
    }
}

//---------------------------------------------------------------------------//
// Produce new String which is concatenation of two others.
//---------------------------------------------------------------------------//

String operator+( const String& x, const String& y )
{
    String s = x;		// Hmm, maybe we should find a way to
    return s += y;		// preallocate space in s for x+y.
}

//---------------------------------------------------------------------------//
// Concatenate a String and a char *.
//---------------------------------------------------------------------------//

String operator+( const String& x, const char *s )
{
    String ss = x;		// ditto.
    return ss += s;
}

//---------------------------------------------------------------------------//
// Concatenate a char * and a String.
//---------------------------------------------------------------------------//

String operator+( const char *s, const String& x )
{
    String ss(s);
    return ss += x;
}

//---------------------------------------------------------------------------//
// Concatenate a String and a char.
//---------------------------------------------------------------------------//

String operator+( const String& s, const char c )
{
    String ss = s;
    return ss += c;
}

//---------------------------------------------------------------------------//
// Concatenate a char and a String.
//---------------------------------------------------------------------------//

String operator+( const char c, const String& s )
{
    String ss(c);
    return ss += s;
}

NAMESPACE_DS_END

//---------------------------------------------------------------------------//
//                              end of String.cc
//---------------------------------------------------------------------------//
