//----------------------------------*-C++-*----------------------------------//
// NestMap.t.hh
// Shawn Pautz
// Mon Mar 29 16:04:51 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "NestMap.hh"

#ifdef SDP

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
NestMap<Key, T>::iterator::iterator
(const stack<typename map<Key, NestMap<Key, T>::Node>::iterator>& parentNodes_)
    : parentNodes(parentNodes_), currentNode(parentNodes.top())
{ parentNodes.pop(); }

template<class Key, class T>
NestMap<Key, T>::iterator::~iterator() {}

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
NestMap<Key, T>::reference NestMap<Key, T>::iterator::operator*() const
{ return currentNode->second.data; }

template<class Key, class T>
NestMap<Key, T>::pointer NestMap<Key, T>::iterator::operator->() const
{ return &(currentNode->second.data); }

template<class Key, class T>
bool NestMap<Key, T>::iterator::operator!=(const iterator& iter) const
{ return (currentNode != iter.currentNode); }

template<class Key, class T>
NestMap<Key, T>::Node::Node(const Key& key_, const T& x_, int level_)
    : data(pair<value_type::first_type,
                value_type::second_type>(key_, x_)),
      level(level_), nestedNodes() {}

template<class Key, class T>
NestMap<Key, T>::Node::~Node() {}

template<class Key, class T>
bool NestMap<Key, T>::Node::empty()
{ return nestedNodes.empty(); }

template<class Key, class T>
typename map<Key, typename NestMap<Key, T>::Node>::iterator
NestMap<Key, T>::Node::begin()
{ return nestedNodes.begin(); }

template<class Key, class T>
typename map<Key, typename NestMap<Key, T>::Node>::iterator
NestMap<Key, T>::Node::end()
{ return nestedNodes.end(); }

template<class Key, class T>
NestMap<Key, T>::NestMap()
    : activeNodes(), nestedNodes(), size(0) {}

template<class Key, class T>
NestMap<Key, T>::NestMap(const NestMap &rhs)
    : activeNodes(rhs.activeNodes),
      nestedNodes(rhs.nestedNodes), size(rhs.size) {}

template<class Key, class T>
NestMap<Key, T>::~NestMap() {}

template<class Key, class T>
NestMap<Key, T>& NestMap<Key, T>::operator=(const NestMap &rhs)
{
    activeNodes = rhs.activeNodes;
    nestedNodes = rhs.nestedNodes;
    size = rhs.size;

    return *this;
}

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::open(const Key& key, T t = T())
{
    map<Key, NestMap<Key, T>::Node>* pnodemap;
    map<Key, NestMap<Key, T>::Node>::iterator newiter;

    if (activeNodes.empty())
	pnodemap = &nestedNodes;
    else
	pnodemap = &((*(activeNodes.top())).second.nestedNodes);
    if (pnodemap->count(key) == 0)
    {
	NestMap<Key, T>::Node newnode(key, t, activeNodes.size() + 1);
	map<Key, NestMap<Key, T>::Node>::value_type newpair(key, newnode);
	newiter = (pnodemap->insert(newpair)).first;
	++size;
    }
    else
	newiter = pnodemap->find(key);
    activeNodes.push(newiter);

    return NestMap<Key, T>::iterator(activeNodes);
}

template<class Key, class T>
void NestMap<Key, T>::close()
{ activeNodes.pop(); }

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::begin()
{
    stack<map<Key, NestMap<Key, T>::Node>::iterator> iterStack;
    iterStack.push(nestedNodes.begin());

    return NestMap<Key, T>::iterator(iterStack);
}

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::end()
{
    stack<map<Key, NestMap<Key, T>::Node>::iterator> iterStack;
    iterStack.push(nestedNodes.end());

    return NestMap<Key, T>::iterator(iterStack);
}

template<class Key, class T>
NestMap<Key, T>::iterator NestMap<Key, T>::current()
{ return NestMap<Key, T>::iterator(activeNodes); }

} // end namespace rtt_dsxx

#endif

//---------------------------------------------------------------------------//
//                              end of NestMap.t.hh
//---------------------------------------------------------------------------//
