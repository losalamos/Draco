//----------------------------------*-C++-*----------------------------------//
// Assert.cc
// Geoffrey Furnish
// Fri Jul 25 08:41:38 1997
//---------------------------------------------------------------------------//
// @> Helper functions for the Assert facility.
//---------------------------------------------------------------------------//

#include "Assert.hh"

#include <string.h>
#include <stdio.h>

#include <iostream.h>

NAMESPACE_DS_BEG

assertion::assertion( const char *cond, const char *file, int line )
{
    msg = new char[ strlen(cond) + strlen(file) + 500 ];
    sprintf( msg, "Assertion: %s, failed in %s, line %8d.",
	     cond, file, line );
}

assertion::assertion( const char *m )
{
    msg = new char[ strlen(m)+1 ];
    strcpy( msg, m );
}

assertion::assertion( const assertion& a )
{
    msg = new char[ strlen(a.msg)+1 ];
    strcpy( msg, a.msg );
}

assertion& assertion::operator=( const assertion& a )
{
    msg = new char[ strlen(a.msg)+1 ];
    strcpy( msg, a.msg );
    return *this;
}

//---------------------------------------------------------------------------//
// Function to perform the task of actually throwing an assertion.
//---------------------------------------------------------------------------//

void toss_cookies( const char *cond, const char *file, int line )
{
// cerr << ...
    throw assertion( cond, file, line );
}

//---------------------------------------------------------------------------//
// Function to perform the task of actually throwing an isistion.
//---------------------------------------------------------------------------//

void insist( const char *cond, const char *msg, const char *file, int line )
{
// cerr << "Insisting that ..." << endl;
    throw assertion( msg );
}

NAMESPACE_DS_END

//---------------------------------------------------------------------------//
//                              end of Assert.cc
//---------------------------------------------------------------------------//
