//----------------------------------*-C++-*----------------------------------//
// NestMap.hh
// Shawn Pautz
// Mon Mar 29 16:04:51 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifdef SDP

#ifndef __ds_NestMap_hh__
#define __ds_NestMap_hh__

#include <cstddef>
#include <map>
#include <stack>
#include <utility>

using std::map;
using std::stack;
using std::pair;

namespace rtt_dsxx
{
 
//===========================================================================//
// class NestMap - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Key, class T>
class NestMap 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef pair<const Key, T> value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const pointer const_pointer;
    typedef ptrdiff_t difference_type;
    typedef size_t size_type;

    class Node;

  private:

    typedef typename map<Key, Node>::iterator mapIterator;

  public:

    class iterator
    {
	// NESTED CLASSES AND TYPEDEFS

	// DATA

	stack<mapIterator> parentNodes;
	typename map<Key, Node>::iterator currentNode;

      public:

	// CREATORS

	iterator(const stack<typename map<Key, Node>::iterator>& parentNodes_);

	~iterator();

	// MANIPULATORS

	iterator& operator++();
	reference operator*() const;
	pointer operator->() const;

	// ACCESSORS

	bool operator!=(const iterator& iter) const;

      private:

	// IMPLEMENTATION
    };

    class Node
    {
	friend class NestMap;
	friend class NestMap::iterator;

	// NESTED CLASSES AND TYPEDEFS

	// DATA

	value_type data;
	int level;
	map<Key, Node> nestedNodes;

      public:

	// CREATORS

	Node(const Key& key_, const T& x_, int level_);

	~Node();

	// MANIPULATORS

	// ACCESSORS

	bool empty();
        typename map<Key, Node>::iterator begin();
        typename map<Key, Node>::iterator end();

      private:

	// IMPLEMENTATION
    };

  private:

    // DATA

    stack<typename map<Key, Node>::iterator> activeNodes;
    map<Key, Node> nestedNodes;
    size_type size;
    
  public:

    // CREATORS
    
    NestMap();
    NestMap(const NestMap &rhs);
    ~NestMap();

    // MANIPULATORS
    
    NestMap& operator=(const NestMap &rhs);

    iterator open(const Key& key, T t = T());
    void close();

    // ACCESSORS

    iterator begin();
    iterator end();
    iterator current();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_dsxx

#endif                          // __ds_NestMap_hh__

#endif

//---------------------------------------------------------------------------//
//                              end of NestMap.hh
//---------------------------------------------------------------------------//
