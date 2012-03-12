//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Endian.hh
 * \author Mike Buksas
 * \date   Tue Oct 23 14:15:55 2007
 * \brief  Function declarations for endian conversions
 * \note   Copyright (C) 2007-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Endian_hh
#define dsxx_Endian_hh

#include <ds++/config.h>
#include <algorithm>
#include <cstring>
#include <stdint.h>
#include <iomanip>
#include <iostream>

//---------------------------------------------------------------------------//
/*!
 * Endian conversion functions.
 *
 * The endian nature of a data representation describes the order in which the
 * constitutent bytes of a multi-byte data value are ordered. We are concerned
 * with converting between big and little endian orderings on platforms where
 * the char data type is one byte in size.
 *
 * If there are other endians out there, I seriously do not want to know about
 * them.
 *
 * To convert between big and little endian data we intrepret the data to be
 * converted as a character array by casting a pointer to the data to
 * (char*). We then manipulate the order, but not the contents, of the
 * character data.
 *
 * Note that we are implicitly assuming that the size of char on each platform
 * is one byte.
 *
 * In order for these functions to work on floating point data, we are
 * assuming that the floating point representations are identical on the two
 * architectures _except_ for the difference in endianness. Also, the sign and
 * exponent information of the floating point representation must fit within a
 * single byte of data, so that it does not require extra steps at the
 * bit-level for conversion.
 *
 */
//---------------------------------------------------------------------------//


namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*!
 * \brief Elemetary byte-swapping routine.
 *
 * \arg The data to byte-swap, represented as character data.
 * \arg The size of the data array.
 *
 * This is a core routine used by other functions to convert data between
 * endian representations.
 *
 * It swaps the elements of a character array of length n. Element 0 is
 * swapped with element n, 1 with n-1 etc... The contents of the individual
 * elements are not changed, only their order.
 *
 * For example, consider the unsigned integer value: \c 0xDEADBEEF.  (\c 0x
 * means this is a hexidecimal value) Two hexidecimal digits is a single byte
 * (16^2 = 2^8) so the layout of the value in big endian style is:
 * \verbatim
 *       0        1        2        3
 *     D  E     A  D     B  E     E  F 
 *  |--------|--------|--------|--------|
 *       ^        ^        ^        ^
 *       |        +--------+        |
 *       +--------------------------+
 *                 swapped
 * \endverbatim
 * The conversion to little endian involves the swap operations pictured in
 * the diagram above. The resulting value (if still interpreted as big-endian)
 * is \c 0xEFBEADDE.
 *
 * We provide two versions for signed and unsigned character data. Internally,
 * we use unsigned. Certain applications use signed char data, and the second
 * form is provided if they need to manipulate the character data directly,
 * instead of using one of the byte_swap functions.
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
/*!
 * \brief General byte-swapping routine
 *
 * This function operates in place on its argument.
 * 
 */
template <typename T>
void byte_swap(T& value)
{
    char_byte_swap((unsigned char*)(&value), sizeof(T));
}


//---------------------------------------------------------------------------//
/*!
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

//---------------------------------------------------------------------------//
/*!
 * \brief Does this platform use big or little endianness
 *
 * \return true if platform uses big endian format
 */
DLL_PUBLIC bool is_big_endian(void);

//---------------------------------------------------------------------------//
/*!
 * \brief Does this platform support IEEE float representation?
 *
 * Some older Cray machines did not support the IEEE float representation.
 * This simple test will identify machines that are IEEE compliant.
 * 
 * \return true if we support IEEE float representation.
 */
DLL_PUBLIC bool has_ieee_float_representation(void);

} // end namespace rtt_dsxx

#endif // dsxx_Endian_hh

//---------------------------------------------------------------------------//
//              end of ds++/Endian.hh
//---------------------------------------------------------------------------//
