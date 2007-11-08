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
 * \arg The data to byte-swap, represented as unsigned character data.
 * \arg The size of the data array.
 *
 * This is a core routine used by other functions with a friendlier
 * interface. 
 * 
 */
void char_byte_swap(unsigned char *data, int n)
{

    unsigned char *end = data+n-1;
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

} // end namespace rtt_dsxx

#endif // dsxx_Endian_hh

//---------------------------------------------------------------------------//
//              end of ds++/Endian.hh
//---------------------------------------------------------------------------//
