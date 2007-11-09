//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Endian.cc
 * \author Mike Buksas
 * \date   Friday, November  9, 2007
 * \brief  Function definitions for Endian
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Endian.hh"

namespace rtt_dsxx
{

void char_byte_swap(unsigned char *data, int n)
{
    unsigned char *end = data+n-1;
    while (data < end) std::swap(*data++, *end--);
}

void char_byte_swap(char *data, int n)
{
    char* end = data+n-1;
    while (data < end) std::swap(*data++, *end--);
}


} // end namespace <namespace>

//---------------------------------------------------------------------------//
//                 end of Endian.cc
//---------------------------------------------------------------------------//
