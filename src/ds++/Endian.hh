//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Endian.hh
 * \author Mike Buksas
 * \date   Tue Oct 23 14:15:55 2007
 * \brief  Function declarations for endian conversions
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Endian_hh
#define dsxx_Endian_hh

#include <algorithm>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/**
 * \brief Elemetary byte-swapping routine.
 *
 * \arg The data to byte-swap, represented as character data.
 * \arg The size of the data array.
 *
 * This is a core routine used by other functions with a friendlier
 * interface.
 *
 * Note that we provide two versions for signed and unsigned character
 * data. Internally, we use unsigned. Certain applications use signed char
 * data, and the second form is provided if they need to manipulate the
 * character data directly, instead of using one of the byte_swap functions. 
 * 
 */
inline void char_byte_swap(unsigned char *data, int n)
{
    unsigned char *end = data+n-1;
    while (data < end) std::swap(*data++, *end--);
}

inline void char_byte_swap(char *data, int n)
{
    char* end = data+n-1;
    while (data < end) std::swap(*data++, *end--);
}

//---------------------------------------------------------------------------//
/**
 * \brief General byte-swapping routine
 *
 * This function operates in place on it's argument.
 * 
 */
template <typename T>
void byte_swap(T& value)
{
    char_byte_swap((unsigned char*)(&value), sizeof(T));
}


//---------------------------------------------------------------------------//
/**
 * \brief General byte-swapping routine.
 *
 * This function returns a bite-swapped copy of the argument.
 *
 */
template <typename T>
T byte_swap_copy(T value) 
{
    byte_swap(value);
    return value;
}

} // end namespace rtt_dsxx

#endif // dsxx_Endian_hh

//---------------------------------------------------------------------------//
//              end of ds++/Endian.hh
//---------------------------------------------------------------------------//
