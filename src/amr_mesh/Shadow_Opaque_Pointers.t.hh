//----------------------------------*-C++-*----------------------------------//
// Shadow_Opaque_Pointers.t.hh
// Mark Gray (original) / B.T. Adams (modified to use smart pointers)
// 1 Sep 99
/*! 
 * \file   amr_mesh/Shadow_Opaque_Pointers.t.hh
 * \author Mark Gray/B.T. Adams
 * \date   Wed 1 Sep 10:33:26 1999
 * \brief  Header file for Shadow_Opaque_Pointers library.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Shadow_Opaque_Pointers.hh"
#include <limits>
#include "assert.h"

template <class T>
opaque_pointers<T>::rep & opaque_pointers<T>::get_rep()
{
    static rep r;
    return r;
}

template <class T>
opaque_pointer_type opaque_pointers<T>::insert(SP<T> t)
{
    // add t to list, return opaque pointer to it

    assert(!opaque_pointers<T>::is_full()); // REQUIRE

    opaque_pointer_type i = get_next_avail();
    
    ptr_map::value_type value(i,t);
    bool succeeded = get_object_pointers().insert(value).second;

    assert(succeeded);
    
    get_next_avail()++;
    
    assert(t == opaque_pointers<T>::item(i)); // ENSURE

    return i;
}

template <class T>
bool opaque_pointers<T>::is_full()
{
    // is there no more room?
    return get_next_avail() >= std::numeric_limits<opaque_pointer_type>::max()
	|| get_object_pointers().size() >= get_object_pointers().max_size();
}

template <class T>
SP<T> opaque_pointers<T>::item(opaque_pointer_type i)
{
    // convert opaque pointer to real pointer
    assert(opaque_pointers<T>::has(i)); // ENSURE

    return get_object_pointers().find(i)->second;
}

template <class T> 
bool opaque_pointers<T>::has(opaque_pointer_type i)
{
    // is i associated?

    ptr_map::size_type nfound = get_object_pointers().count(i);
    return nfound != 0;
}

template <class T>
void opaque_pointers<T>::erase(opaque_pointer_type i)
{
    // remove pointer referenced by opaque pointer
    
    assert(opaque_pointers<T>::has(i)); // REQUIRE

    get_object_pointers().erase(i);
    
    assert(!opaque_pointers<T>::has(i)); // ENSURE
}
