//----------------------------------*-C++-*----------------------------------//
// Algorithms.hh
// Geoffrey Furnish
// Mon Jan 27 12:11:50 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __ds_Algorithms_hh__
#define __ds_Algorithms_hh__

//---------------------------------------------------------------------------//
// Search through two ranges of memory, checking to see if the data they hold
// are equal.
//---------------------------------------------------------------------------//

template<class FwdItr>
bool range_equal( FwdItr s1, FwdItr e1, FwdItr s2 )
{
    for( ; s1 != e1; )
	if (*s1++ != *s2++) return false;

    return true;
}

template<class T>
T norm( T *begin, T *end )
{
    T r = T(0.);
    for( T *x = begin; x != end; x++ ) r += (*x) * (*x);
    return sqrt( r / T(end - begin) );
}

#endif                          // __ds_Algorithms_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Algorithms.hh
//---------------------------------------------------------------------------//
