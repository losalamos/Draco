#ifndef OPAQUE_POINTERS_HH
#define OPAQUE_POINTERS_HH

#include <map>

using dsxx::SP;

typedef long opaque_pointer_type;

template <class T>
class opaque_pointers
{
  public:

    typedef T value_type;
    typedef T *pointer_type;
    typedef T &reference;
    
    static opaque_pointer_type insert(SP<T>);
	// add t to list, return opaque pointer to it

    static bool is_full();
	// is there no more room?

    static SP<T> item(opaque_pointer_type i);
	// convert opaque pointer to real pointer

    static bool has(opaque_pointer_type i);
	// is i associated?

    static void erase(opaque_pointer_type i);
	// remove pointer referenced by opaque pointer

  private:

    typedef std::map<opaque_pointer_type, SP<T> > ptr_map;

    struct rep
    {
	opaque_pointer_type next_avail;
	ptr_map object_pointers;
    };

    static rep &get_rep();

    static ptr_map &get_object_pointers()
    {
	return get_rep().object_pointers;
    }

    static opaque_pointer_type &get_next_avail()
    {
	return get_rep().next_avail;
    }

};

// include implementation

#include "Shadow_Opaque_Pointers.t.hh"

#endif
