//----------------------------------*-C++-*----------------------------------//
// Shadow_Opaque_Pointers.hh
// Mark Gray (original) / B.T. Adams (modified to use smart pointers)
// 1 Sep 99
/*! 
 * \file   ds++/Shadow_Opaque_Pointers.hh
 * \author Mark Gray/B.T. Adams
 * \date   Wed 1 Sep 10:33:26 1999
 * \brief  Header file for opaque_pointers library.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//
#ifndef OPAQUE_POINTERS_HH
#define OPAQUE_POINTERS_HH

#include <map>

/*!
 * \brief RTT shadow object interface namespace.
 *
 * Provides namespace protection for the Draco RTT shadow object interface
 * opaque pointer utilities that are used to interface C++ objects with 
 * external codes that use other programming languages. An example code is 
 * also provide to illustrate the usage of all of the shadow object interface
 * functions to the amr_mesh package from a Fortran 90 code.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
namespace rtt_dsxx
{


/*! 
 * \brief  Defines the C++ data type that will be used to represent the opaque
 *         pointers (and will thus be exchanged with the external code in all 
 *         opaque pointer references).
 */
typedef long opaque_pointer_type;

/*! 
 * \brief  The templated opaque_pointers class provides utilities that can be 
 *         used to manage (i.e., create, destroy, execute, and access) C++ 
 *         objects from a Fortran 90 (or other programming language) external
 *         code. The opaque pointers that are used herein are a C++ data type
 *         (e.g., long) that can be exchanged with a F90 code in calling 
 *         arguments and are used internally as a "key" to reference smart
 *         pointers to the C++ objects of the specified type.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
template <class T>
class opaque_pointers
{
  public:

	// add t to list, return opaque pointer to it
/*! 
 * \brief  Stores the smart pointer to the specified C++ object in an 
 *         opaque_pointer class object of type T and returns the new opaque 
 *         pointer to the object. This also updates the next available opaque 
 *         pointer number for objects of this data type.
 * \param t Smart pointer to a C++ object of data type T.
 * \return  An opaque_pointer_type object that references the smart pointer 
 *          to the C++ object.
 */
    static opaque_pointer_type insert(SP<T> t);

	// is there no more room?
/*! 
 * \brief  Checks to insure that an opaque pointer can be added to the
 *         opaque_pointer class object of this type without exceeding
 *         storage limitations.
 * \return Status of existing storage < storage limit.
 */
    static bool is_full();

	// convert opaque pointer to real pointer
/*! 
 * \brief  Returns the smart pointer to the C++ object that corresponds to 
 *         the specified opaque pointer.
 * \param i Opaque pointer number.
 * \return Smart pointer to the C++ object.
 */
    static SP<T> item(opaque_pointer_type i);

	// is i associated?
/*! 
 * \brief  Checks to insure that the specified opaque pointer is assigned 
 *         to a C++ object.
 * \param i Opaque pointer number.
 * \return Status of the opaque pointer assignment to a C++ object.
 */
    static bool has(opaque_pointer_type i);

	// remove pointer referenced by opaque pointer
/*! 
 * \brief  Removes the smart pointer to the referenced C++ object from the 
 *         opaque_pointer class object container and eliminates the associated
 *         opaque pointer "key".
 * \param i Opaque pointer number.
 */
    static void erase(opaque_pointer_type i);

  private:

/*! 
 * \brief  Defines the C++ container type that will be used to store the 
 *         opaque pointer "keys" and the associated smart pointers to the 
 *         C++ objects.
 */
    typedef std::map<opaque_pointer_type, SP<T> > ptr_map;

/*! 
 * \brief  Defines a C++ data structure that is used to store the 
 *         opaque_pointers container and the next available opaque pointer 
 *         number for C++ objects of a specific type.
 */
    struct rep
    {
/*! 
 * \brief  The next available opaque pointer number for C++ objects of a 
 *         specific type.
 */
	opaque_pointer_type next_avail;
/*! 
 * \brief  The opaque_pointers container for the opaque pointers and the 
 *         associated smart pointers to the C++ objects.
 */
	ptr_map object_pointers;
    };

/*! 
 * \brief  Returns the C++ data structure that is used to store the 
 *         opaque_pointers container and the next available opaque pointer 
 *         number for C++ objects of a specific type.
 */
    static rep & get_rep();

/*! 
 * \brief  Returns the C++ container from a type rep structure.
 * \return The opaque_pointers container for C++ objects of this type.
 */
    static ptr_map & get_object_pointers()
    {
	return get_rep().object_pointers;
    }

/*! 
 * \brief  Returns the next available opaque pointer number from a type rep 
 *         structure.
 * \return The next available opaque pointer number for C++ objects of this 
 *         type.
 */
    static opaque_pointer_type & get_next_avail()
    {
	return get_rep().next_avail;
    }

};

// include implementation

#include "opaquePointers.t.hh"

} // end namespace rtt_dsxx

#endif
