// swap.hh
// T. M. Kelley
// Jun 04, 2007
// Header for swap
// (c) Copyright 2007 LANSLLC all rights reserved.

#ifndef SWAP_H
#define SWAP_H

#include <stdint.h>

namespace rtt_utils
{
    // Do byte-swapping one of two ways: either use GNU extended asm
    // for x86 (really 486+), or use the "poor people's" method of 
    // digging out one byte at a time and moving it to the right place.
    // The poor method is really not that bad: GCC for example seems to 
    // optimize it down to about 5 instructions in the 32 bit case, 
    // including loads, versus 2 instructions (including load) for the 
    // inline asm case.
    inline uint32_t swap(uint32_t const input)
    {
#ifdef __use_x86_gnu_asm
        uint32_t output = input;
        asm("bswap %0" : "+g" (output) : );
#else
        uint32_t byte, output;
        byte = input & 255U;
        output = (byte << 24);

        byte = input & 65280U;      // 255 << 8
        output = output | (byte << 8);
        
        byte = input & 16711680U;   // 255 << 16
        output = output | (byte >> 8);
        
        byte = input & 4278190080U; // 255 << 24
        output = output | (byte >> 24); // look out--algebraic shift r.
        
#endif // __use_x86_asm
        return output;
    } // int32_t swap( int32_t)


    inline double swap( double const input)
    {
#ifdef __use_x86_gnu_asm
        double output = input;
        asm("bswap %0" : "+g" (output) : );   
#else
        union
        {
            double d;
            uint64_t u;
        } b64;

        uint64_t byte, tmp, uinput;

        // change meaning of input bits to uint64_t:
        b64.d = input;
        uinput = b64.u;

        // 1
        byte = uinput & 255;
        tmp = (byte << 56);
        // 2
        byte = uinput & 65280;      // 255 << 8
        tmp = tmp | (byte << 40);
        // 3
        byte = uinput & 16711680;   // 255 << 16
        tmp = tmp | (byte << 24);
        // 4
        byte = uinput & 4278190080U; // 255 << 24
        tmp = tmp | (byte << 8); 
        // 5
        byte = uinput & 1095216660480ULL; // 255 << 32
        // byte = uinput & 1095216660480; // 255 << 32
        tmp = tmp | (byte >> 8);
        // 6 
        byte = uinput & 280375465082880ULL; // 255 << 40
        tmp = tmp | (byte >> 24);
        // 7
        byte = uinput & 71776119061217280ULL; // 255 << 48
        tmp = tmp | (byte >> 40);
        // 8
        byte = uinput & 18374686479671623680ULL; // 255 << 56
        // byte = uinput & 0xff00000000000000; // 255 << 56
        tmp = tmp | (byte >> 56); 
    
        // change meaning of bits in b64.
        b64.u = tmp;
        double output = b64.d;
#endif // __use_x86_gnu_asm
        
        return output;
    } // double swap( double)


} // rtt_utils::

#endif // include guard


// version
// $Id$

// End of file
