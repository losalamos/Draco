//----------------------------------*-C++-*----------------------------------//
// NestMap.t.hh
// Shawn Pautz
// Mon Mar 29 16:04:51 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __ds_NestMap_hh__
#define __ds_NestMap_hh__

#include <cstddef>
#include <map>
#include <stack>
#include <utility>

using std::map;
using std::stack;
using std::pair;

namespace rtt_ds
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
	mapIterator currentNode;

      public:

	// CREATORS

	iterator(const stack<mapIterator>& parentNodes_)
	    : parentNodes(parentNodes_), currentNode(parentNodes.top())
	{ parentNodes.pop(); }

	~iterator() {}

	// MANIPULATORS

	iterator& operator++();
	reference operator*() const { return currentNode->second.data; }
	pointer operator->() const { return &(currentNode->second.data); }

	// ACCESSORS

	bool operator!=(const iterator& iter) const
	{ return (currentNode != iter.currentNode); }

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

	Node(const Key& key_, const T& x_, int level_)
	    : data(pair<value_type::first_type,
                        value_type::second_type>(key_, x_)),
	      level(level_), nestedNodes() {}

	~Node() {}

	// MANIPULATORS

	// ACCESSORS

	bool empty() { return nestedNodes.empty(); }
        mapIterator begin() { return nestedNodes.begin(); }
        mapIterator end() { return nestedNodes.end(); }

      private:

	// IMPLEMENTATION
    };

  private:

    // DATA

    stack<mapIterator> activeNodes;
    map<Key, Node> nestedNodes;
    size_type size;
    
  public:

    // CREATORS
    
    NestMap() : activeNodes(), nestedNodes(), size(0) {}
    NestMap(const NestMap &rhs) : activeNodes(rhs.activeNodes),
	nestedNodes(rhs.nestedNodes), size(rhs.size) {}
    ~NestMap() {}

    // MANIPULATORS
    
    NestMap& operator=(const NestMap &rhs)
    {
	activeNodes = rhs.activeNodes;
	nestedNodes = rhs.nestedNodes;
	size = rhs.size;
    }

    iterator open(const Key& key, T t = T());
    void close() { activeNodes.pop(); }

    // ACCESSORS

    iterator begin();
    iterator end();
    iterator current() { return iterator(activeNodes); }

  private:
    
    // IMPLEMENTATION
};

template<class Key, class T>
NestMap<Key, T>::iterator& NestMap<Key, T>::iterator::operator++()
{
    if (currentNode->second.empty())
    {
	while(!parentNodes.empty() &&
              ++currentNode == (parentNodes.top())->second.end())
	{
	    currentNode = parentNodes.top();
	    parentNodes.pop();
	}
	if (parentNodes.empty())
	    ++currentNode;
    }
    else
    {
	parentNodes.push(currentNode);
	currentNode = currentNode->second.begin();
    }

    return *this;
}

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::open(const Key& key, T t = T())
{
    map<Key, Node>* pnodemap;
    mapIterator newiter;

    if (activeNodes.empty())
	pnodemap = &nestedNodes;
    else
	pnodemap = &((*(activeNodes.top())).second.nestedNodes);
    if (pnodemap->count(key) == 0)
    {
	Node newnode(key, t, activeNodes.size() + 1);
	map<Key, Node>::value_type newpair(key, newnode);
	newiter = (pnodemap->insert(newpair)).first;
	++size;
    }
    else
	newiter = pnodemap->find(key);
    activeNodes.push(newiter);

    return iterator(activeNodes);
}

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::begin()
{
    stack<mapIterator> iterStack;
    iterStack.push(nestedNodes.begin());

    return iterator(iterStack);
}

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::end()
{
    stack<mapIterator> iterStack;
    iterStack.push(nestedNodes.end());

    return iterator(iterStack);
}

} // end namespace rtt_ds

#endif                          // __ds_NestMap_hh__

//---------------------------------------------------------------------------//
//                              end of NestMap.t.hh
//---------------------------------------------------------------------------//
