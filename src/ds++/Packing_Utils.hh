//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Packing_Utils.hh
 * \author Thomas M. Evans
 * \date   Thu Jul 19 11:27:46 2001
 * \brief  Packing Utilities, classes for packing stuff.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ds_Packing_Utils_hh__
#define __ds_Packing_Utils_hh__

#include "Assert.hh"

#include <cstring>

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \file ds++/Packing_Utils.hh

 * This file contains classes and utilities that are used to "pack" data into
 * byte-streams.  The byte-streams are represented by the char* type.  The
 * following classes are:

 * \arg \b Packer packing class
 * \arg \b Unpacker unpacking class

 */

/*!
 * \example ds++/test/tstPacking_Utils.cc

 * Test the Packer and Unpacker classes.

 */
// revision history:
// -----------------
// 0) original
//===========================================================================//
 
//===========================================================================//
/*!
 * \class Packer
 
 * \brief Pack data into a byte stream.

 * This class allows clients to "register" a char* stream and then load it
 * with data of any type.  This assumes that the sizeof(T) operator works and
 * has meaning for the type.  Under the hood it uses std::memcpy to perform
 * the loading.  This class is easily understood by checking the examples.

 * No memory allocation is performed by the Packer. 

 * The benefit of using the Packer class is that byte copies are isolated
 * into this one section of code, thus obviating the need for
 * reinterpret_cast statements in client code.  In fact, this functionality
 * conforms exactly to the ANSI C++ standard for copying byte-streams of data
 * (sec. 3.9). Additionally, bounds checking is performed on all stream
 * packing operations.  This bounds checking is always on.

 * This class returns real char * pointers through its query functions.  We
 * do not use the STL iterator notation, even though that is how the pointers
 * are used, so as not to confuse the fact that these char * streams are \i
 * continuous \i data byte-streams.  The pointers that are used to "iterate"
 * through the streams are real pointers, not an abstract iterator class.  So
 * one could think of these as iterators (they act like iterators) but they
 * are real pointers into a continguous memory char * stream.

 * Data can be unpacked using the Unpacker class.
 
 */
//===========================================================================//

class Packer
{
  public:
    // Typedefs.
    typedef char *       pointer;
    typedef const char * const_pointer;

  private:
    // Size of packed stream.
    unsigned int stream_size;
    
    // Pointer (mutable) into data stream.
    pointer ptr;

    // Pointers to begin and end of buffers.
    pointer begin_ptr;
    pointer end_ptr;

  public:
    //! Constructor.
    Packer() : stream_size(0), ptr(0), begin_ptr(0), end_ptr(0) {/*...*/}

    // Set the buffer.
    inline void set_buffer(unsigned int, pointer);

    // Pack values into the buffer.
    template<class T> inline void pack(const T&);

    // >>> ACCESSORS

    //! Get a pointer to the current position of the data stream.
    const_pointer get_ptr() const { return ptr; }

    //! Get a pointer to the beginning position of the data stream.
    const_pointer begin() const { return begin_ptr; }

    //! Get a pointer to the ending position of the data stream.
    const_pointer end() const { return end_ptr; }

    //! Get the size of the data stream.
    unsigned int size() const { return stream_size; }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Set an allocated buffer to write data into.

 * This function accepts an allocated char* buffer.  It assigns begin and end
 * pointers and a mutable position pointer that acts like an iterator.  The
 * Packer will write POD (Plain Old Data) data into this buffer starting at
 * the beginning address of the buffer.  This function must be called before
 * any Packer::pack calls can be made.

 * Once Packer::set_buffer is called, all subsequent calls to Packer::pack
 * will write data incrementally into the buffer set by set_buffer.  To write
 * data into a different buffer, call Packer::set_buffer again; at this point
 * the Packer no longer has any knowledge about the old buffer.

 * Note, the buffer must be allocated large enough to hold all the data that
 * the client intends to load into it.  There is no memory allocation
 * performed by the Packer class; thus, the buffer cannot be increased in
 * size if a value is written past the end of the buffer.  See the
 * Packer::pack function for more details.

 * \param size_in size of the buffer
 * \param buffer pointer to the char * buffer

 */
void Packer::set_buffer(unsigned int size_in, pointer buffer)
{
    Require (buffer);
    
    // set the size, begin and end pointers, and iterator
    stream_size  = size_in;
    ptr          = &buffer[0];
    begin_ptr    = &buffer[0];
    end_ptr      = begin_ptr + stream_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack data into a buffer.

 * This function packs a piece of data (single datum) into the buffer set by
 * Packer::set_buffer.  It advances the pointer (iterator) location to the
 * next location automatically.  It uses the sizeof(T) operator to get the
 * size of the data; thus, only data where sizeof() has meaning will be
 * properly written to the buffer.

 * Packer::pack() does bounds checking to ensure that the buffer and buffer
 * size defined by Packer::set_buffer are consistent.  This bounds-checking
 * is always on as the Packer is not normally used in compute-intensive
 * calculations.

 * \param value data of type T to pack into the buffer; the data size must be
 * accessible using the sizeof() operator

 */
template<class T>
void Packer::pack(const T &value)
{
    Require (begin_ptr);
    Insist  (ptr >= begin_ptr, "Bounds error in packer!");
    Insist  (ptr + sizeof(T) <= end_ptr, "Bounds error in packer!");

    // copy value into the buffer
    std::memcpy(ptr, &value, sizeof(T));

    // advance the iterator pointer to the next location
    ptr += sizeof(T);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream out (<<) operator for packing data.

 * The overloaded stream out operator can be used to pack data into streams
 * (Packer p; p.set_buffer(i,b); p << data;).  It simply calls the
 * Packer::pack function.  It returns a reference to the Packer object so
 * that stream out operations can be strung together.

 */
template<class T>
inline Packer& operator<<(Packer &p, const T &value)
{
    // pack the value
    p.pack(value);

    // return the packer object
    return p;
}

//===========================================================================//
/*!
 * \class Unpacker
 
 * \brief Unpack data from a byte stream.

 * This class allows clients to "register" a char* stream and then unload
 * data from it.  This assumes that the sizeof(T) operator works and has
 * meaning for the type.  Under the hood it uses std::memcpy to perform the
 * unloading.  This class is easily understood by checking the examples.

 * No memory allocation is performed by the Unpacker. 

 * The benefit of using the Unpacker class is that byte copies are isolated
 * into this one section of code, thus obviating the need for
 * reinterpret_cast statements in client code.  In fact, this functionality
 * conforms exactly to the ANSI C++ standard for copying byte-streams of data
 * (sec. 3.9). Additionally, bounds checking is performed on all stream
 * packing operations.  This bounds checking is always on.

 * This class returns real char * pointers through its query functions.  We
 * do not use the STL iterator notation, even though that is how the pointers
 * are used, so as not to confuse the fact that these char * streams are \i
 * continuous \i data byte-streams.  The pointers that are used to "iterate"
 * through the streams are real pointers, not an abstract iterator class.  So
 * one could think of these as iterators (they act like iterators) but they
 * are real pointers into a continguous memory char * stream.

 * This class is the complement to the Packer class.
 
 */
//===========================================================================//

class Unpacker
{
  public:
    // Typedefs.
    typedef char *       pointer;
    typedef const char * const_pointer;

  private:
    // Size of packed stream.
    unsigned int stream_size;
    
    // Pointer (mutable) into data stream.
    const_pointer ptr;

    // Pointers to begin and end of buffers.
    const_pointer begin_ptr;
    const_pointer end_ptr;

  public:
    //! Constructor.
    Unpacker() : stream_size(0), ptr(0), begin_ptr(0), end_ptr(0) {/*...*/} 

    // Set the buffer.
    inline void set_buffer(unsigned int, const_pointer);

    // Unpack value from buffer.
    template<class T> inline void unpack(T &);
    
    // >>> ACCESSORS

    //! Get a pointer to the current position of the data stream.
    const_pointer get_ptr() const { return ptr; }

    //! Get a pointer to the beginning position of the data stream.
    const_pointer begin() const { return begin_ptr; }

    //! Get a pointer to the ending position of the data stream.
    const_pointer end() const { return end_ptr; }

    //! Get the size of the data stream.
    unsigned int size() const { return stream_size; }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Set an allocated buffer to read data from.

 * This function accepts an allocated char* buffer.  It assigns begin and end
 * pointers and a mutable position pointer that acts like an iterator.  The
 * Unpacker will read POD data from this buffer starting at the beginning
 * address of the buffer.  This function must be called before any
 * Unpacker::unpack calls can be made.

 * Once Unpacker::set_buffer is called, all subsequent calls to
 * Unpacker::unpack will read data incrementally from the buffer set by
 * set_buffer.  To read data from a different buffer, call
 * Unpacker::set_buffer again; at this point the Unpacker no longer has any
 * knowledge about the old buffer.

 * Note, there is no memory allocation performed by the Unacker class.  Also,
 * the client must know how much data to read from the stream (of course
 * checks can be made telling where the end of the stream is located using
 * the Unpacker::get_ptr, Unpacker::begin, and Unpacker::end functions).

 * \param size_in size of the buffer
 * \param buffer const_pointer to the char * buffer

 */
void Unpacker::set_buffer(unsigned int size_in, const_pointer buffer)
{
    Require (buffer);
    
    // set the size, begin and end pointers, and iterator
    stream_size  = size_in;
    ptr          = &buffer[0];
    begin_ptr    = &buffer[0];
    end_ptr      = begin_ptr + stream_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack data from the buffer.

 * This function unpacks a piece of data (single datum) from the buffer set
 * by Unpacker::set_buffer.  It advances the pointer (iterator) location to
 * the next location automatically.  It uses the sizeof(T) operator to get
 * the size of the data; thus, only data where sizeof() has meaning will be
 * properly read from the buffer.POLYNOMIAL_Specific_Heat_ANALYTIC_EoS_MODEL

 * Unpacker::unpack() does bounds checking to ensure that the buffer and
 * buffer size defined by Unpacker::set_buffer are consistent.  This
 * bounds-checking is always on as this should not be used in computation
 * intensive parts of the code.

 * \param value data of type T to unpack from the buffer; the data size must
 * be accessible using the sizeof() operator

 */
template<class T>
void Unpacker::unpack(T &value)
{
    Require (begin_ptr);
    Insist  (ptr >= begin_ptr, "Bounds error in unpacker!");
    Insist  (ptr + sizeof(T) <= end_ptr, "Bounds error in unpacker!");

    // copy data into the value reference
    std::memcpy(&value, ptr, sizeof(T));

    // advance the iterator pointer to the next location
    ptr += sizeof(T);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream in (>>) operator for unpacking data.

 * The overloaded stream in operator can be used to unpack data from streams
 * (Unpacker u; u.set_buffer(i,b); u >> data;).  It simply calls the
 * Unpacker::unpack function.  It returns a reference to the Unpacker object
 * so that stream in operations can be strung together.

 */
template<class T>
inline Unpacker& operator>>(Unpacker &u, T &value)
{
    // unpack the value
    u.unpack(value);

    // return the unpacker object
    return u;
}

} // end namespace rtt_dsxx

#endif                          // __ds_Packing_Utils_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Packing_Utils.hh
//---------------------------------------------------------------------------//
