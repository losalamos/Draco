//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/StructuredMeshCFAVTest.t.hh
 * \author Randy M. Roberts
 * \date   Mon Sep 13 11:00:32 1999
 * \brief  Implementation file for the StructuredMeshCFAVTest class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "StructuredMeshCFAVTest.hh"
#include "compose.hh"

#include <sstream>
#include <algorithm>
#include <stdexcept>

//*****************************************
// Constructor: StructuredMeshCFAVTest(...)
//*****************************************

namespace rtt_meshTest
{
template<class MTFactory>
StructuredMeshCFAVTest<MTFactory>::
StructuredMeshCFAVTest(Tester &parent_in, MTFactory &meshFactory_in)
    : Tester("StructuredMeshCFAVTest", parent_in.os()),
      parent_m(parent_in),
      meshFactory_m(meshFactory_in),
      meshProduct_m(meshFactory_m.create()),
      mesh_m(meshProduct_m.mesh()),
      fCtor_m(meshProduct_m.fieldConstructor()),
      vindices_m(fCtor_m),
      findices_m(fCtor_m)
{
    // Set up the vindices_m field to contain the vertex
    // numbers, numbered from 0 (on each process).
    
    int ivrtx = 0;
    for (MT::vcif::iterator vit = vindices().begin();
	 vit != vindices().end();
	 vit++)
    {
	*vit = ivrtx++;
    }

    // Set up the findices_m field to contain the face
    // numbers, numbered from 0 (on each process).
    
    int iface = 0;
    for (MT::fcdif::iterator fit = findices().begin();
	 fit != findices().end();
	 fit++)
    {
	*fit = iface++;
    }
    
}

//*******
// run()
//*******

template<class MTFactory>
void StructuredMeshCFAVTest<MTFactory>::run()
{
    setPassed(true);
    
    os() << "Running t1()." << std::endl;
    
    t1();

    typedef typename MT::ConnFacesAroundVertices<typename MT::fcdif> ConnFcdif;

    os() << "Running t2<non-const>()." << std::endl;
    
    t2<ConnFcdif, ConnFcdif::iterator,
	ConnFcdif::value_type::iterator>("non-const");
    
    os() << "Running t2<const>()." << std::endl;
    
    t2<const ConnFcdif, ConnFcdif::const_iterator,
	ConnFcdif::value_type::const_iterator>("const");
    
    os() << "Running t3<non-const>()." << std::endl;
    
    t3<ConnFcdif, ConnFcdif::iterator,
	ConnFcdif::value_type::iterator>("non-const");
    
    os() << "Running t3<const>()." << std::endl;
    
    t3<const ConnFcdif, ConnFcdif::const_iterator,
	ConnFcdif::value_type::const_iterator>("const");
}

//*******
// t1()
//*******

template<class MTFactory>
void StructuredMeshCFAVTest<MTFactory>::t1()
{
    // Test if the sizes of the vertex-centered field is the size we
    // think it should be.

    const std::string prefix("t1: ");
    
    if (vindices().size() != NumVerticesPerCell() * get_ncells())
    {
	std::ostringstream ost;
	ost << "Size of MT::vcif != " << NumVerticesPerCell()
	    << " * ncells";
	testassert(false, prefix+ost.str(), __FILE__, __LINE__);
    }
    
    // Test if the sizes of the face-centered-disc field is the size we
    // think it should be.

    if (findices().size() != NumFacesPerCell() * get_ncells())
    {
	std::ostringstream ost;
	ost << "Size of MT::fcdif != " << NumFacesPerCell()
	    << " * ncells";
	testassert(false, prefix+ost.str(), __FILE__, __LINE__);
    }
    
}

//********
// t2(...)
//********

template<class MTFactory>
template<class CT, class CIT, class CVIT>
void StructuredMeshCFAVTest<MTFactory>::t2(const std::string &constness)
{
    const std::string prefix(std::string("t2<") + constness + ">: ");
    
    typename MT::fcdif ffield(fCtor_m);
    
    CT connFfield(ffield);

    // The outer iteration of the connFcdif should be the same size as
    // a vertex field.
    
    bool samesize = std::distance(connFfield.begin(), connFfield.end())
	== vindices().size();
    if (!samesize)
    {
	testassert(false,
		   prefix + "Size of ConnFacesAroundVertices<fcdif> field " +
		   " != size of vcif", __FILE__, __LINE__);
    }

    // Let's start the outter iteration, to check the inner iterations.
    
    for (CIT cfi = connFfield.begin(); cfi != connFfield.end(); cfi++)
    {
	// The inner iteration of the connFcdif should be number of
	// faces arround the interior of a vertex.
	
	bool correctSize = std::distance(cfi->begin(), cfi->end())
	    == NumFacesPerVertex();
	if (!correctSize)
	{
	    testassert(false,
		       prefix + "Size of ConnFacesAroundVertices<fcdif> "
		       + "inner iteration != number of faces arround a vertex",
		       __FILE__, __LINE__);
	    return;
	}
    }
}

//********
// t3(...)
//********

template<class MTFactory>
template<class CT, class CIT, class CVIT>
void StructuredMeshCFAVTest<MTFactory>::t3(const std::string &constness)
{
    const std::string prefix(std::string("t3<") + constness + ">: ");

    // Construct two data structures from the ConnFacesArroundVertices
    // construct.
    //
    // The first data structure, facesArroundVertices, contains
    // the set of faces surrounding each vertex of the mesh.
    // The second data structure, verticesArroundFaces, contains
    // the set of vertices contained by each face of the mesh.
    
    VecSetInt facesArroundVertices(vindices().size());
    VecSetInt verticesArroundFaces(findices().size());

    CT connFindices(findices());
    
    typename MT::vcif::const_iterator vi = vindices().begin();
    
    for (CIT cfi = connFindices.begin(); cfi != connFindices.end(); cfi++, vi++)
    {
	const int vindex = *vi;
	
	for (CVIT cfivi = cfi->begin(); cfivi != cfi->end(); cfivi++)
	{
	    const int findex = *cfivi;

	    // Insert the face into this vertex's face set,
	    // and check to see that this faces has not already been
	    // seen by this vertex.

	    if (facesArroundVertices[vindex].insert(findex).second == false)
	    {
		std::ostringstream ost;
		ost << "Duplicate face, "
		    << findex << ", found for vertex " << vindex;
		
		testassert(false, prefix+ost.str(), __FILE__, __LINE__);
		return;
	    }
	    
	    // Insert the vertex into this face's vertex set,
	    // and check to see that this face has not already been
	    // seen by this vertex

	    if (verticesArroundFaces[findex].insert(vindex).second == false)
	    {
		std::ostringstream ost;
		ost << "Duplicate vertex, "
		    << vindex << ", found for face " << findex;
		
		testassert(false, prefix+ost.str(), __FILE__, __LINE__);
		return;
	    }
	}
    }

    // Make sure that the outer iteration of Conn...
    // agrees with the iteration of the vertex centered field.
    
    if (vi != vindices().end())
	testassert(false, prefix + "Conn outer iteration did " +
		   "not agree with vindice's iteration.", __FILE__, __LINE__);

    // Check if this data structure signify the correct number of
    // faces per vertex.

    checkSize(facesArroundVertices, NumFacesPerVertex(),
	      "number of faces per vertex");

    // Check if this data structure signify the correct number of
    // vertices per face.

    checkSize(verticesArroundFaces, NumVerticesPerFace(),
	      "number of vertices per face");

    
    // Check to make sure these data structures make sense
    // as a structured mesh.
    
    checkConnectivity(facesArroundVertices, verticesArroundFaces);
}

//***********************
// checkConnectivity(...)
//***********************

template<class MTFactory>
void StructuredMeshCFAVTest<MTFactory>::
checkConnectivity(const VecSetInt &facesArroundVertices,
		  const VecSetInt &verticesArroundFaces)
{
    const std::string prefix("checkConnectivity");
    
    // From the sets of vertices surrounding each face
    // derive the set of edges in the mesh.

    std::set<Cell> edgeSet = getEdgeSet(verticesArroundFaces);

    // Make sure we have the correct number of edges in the mesh.

    if (edgeSet.size() != NumEdgesPerCell() * get_ncells())
    {
	std::ostringstream ost;
	ost << "Unexpected number of edges in edge set."
	    << " edge set size: " << edgeSet.size()
	    << ", expected number of edges: "
	    << NumEdgesPerCell() * get_ncells();
	
	testassert(false, prefix+ost.str(), __FILE__, __LINE__);
	return;
    }

    // Check that each vertex has the correct number of edges
    // sticking out of it.
    
    checkEdgeSetforVertexCount(edgeSet);
}

//***********************
// getEdgeSet(...)
//***********************

template<class MTFactory>
std::set<Cell> StructuredMeshCFAVTest<MTFactory>::
getEdgeSet(const VecSetInt &verticesArroundFaces)
{
    const std::string prefix("getEdgeSet");
    
    // We shall loop over all of the vertex sets for each pair of faces
    // to determine the intersections of vertices for each pair.
    // The intersecting vertices for a pair of faces determine an edge.
    // We shall return this set of edges.
    
    std::set<Cell> edgeSet;

    for (int iface1 = 0; iface1 < verticesArroundFaces.size(); iface1++)
    {
	for (int iface2 = iface1 + 1;
	     iface2 < verticesArroundFaces.size();
	     iface2++)
	{
	    std::vector<int> vint =
		findIntersection(verticesArroundFaces[iface1],
				 verticesArroundFaces[iface2]);

	    // If there is an intersection of vertices for two faces,
	    // Add the resulting edge to the edge set.
	    
	    if (!vint.empty())
	    {
		try
		{
		    // An edge is a Cell with dimension - 2 from what we
		    // consider a cell.

		    Cell edge(Dimension()-2, vint);
	    
		    if (edgeSet.insert(edge).second == false)
		    {
			std::ostringstream ost;
			ost << "Edge: (";
			std::vector<int>::iterator last = vint.begin();
			std::advance(last, vint.size()-1);
			std::copy(vint.begin(), last,
				  std::ostream_iterator<int>(ost, ","));
			ost << *last << ")"
			    << " already inserted into edge set.";
			testassert(false, prefix+ost.str(), __FILE__, __LINE__);
			return edgeSet;
		    }
		}
		catch (const std::invalid_argument &exc)
		{
		    testassert(false, prefix+exc.what(), __FILE__, __LINE__);
		    return edgeSet;
		}
	    }
	}
    }

    return edgeSet;
}

//***************************
// findIntersection(...)
//***************************

template<class MTFactory>
std::vector<int>
StructuredMeshCFAVTest<MTFactory>::findIntersection(const std::set<int> &set1,
						    const std::set<int> &set2)
{
    std::vector<int> vint;

    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
			  std::back_inserter(vint));

    return vint;
}

//********************************
// checkEdgeSetforVertexCount(...)
//********************************

template<class MTFactory>
void StructuredMeshCFAVTest<MTFactory>::
checkEdgeSetforVertexCount(const std::set<Cell> &edgeSet)
{
    const std::string prefix("checkEdgeSetforVertexCount");
    
    // Check to make sure that each vertex shows up in the edge set
    // the correct number of times.
    // The correct number of times is equal to the dimension of the hex,
    // i.e. a three-dimensional mesh should have each vertex subtend
    // three edges.
    
    std::vector<int> vrtxCount(vindices().size());

    for (std::set<Cell>::const_iterator eit = edgeSet.begin();
	 eit != edgeSet.end();
	 eit++)
    {
	const Cell::VertexList &vlist = eit->vertices();
	for (Cell::VertexList::const_iterator vit = vlist.begin();
	     vit != vlist.end();
	     vit++)
	{
	    if (*vit < 0 || *vit >= vindices().size())
	    {
		testassert(false,
			   prefix+"Found edge with illegal vertex number.",
			   __FILE__, __LINE__);
		return;
	    }
	    vrtxCount[*vit]++;
	}
    }

    // Find the vertex that does not have the correct number of edges.
    
    std::vector<int>::const_iterator badVrtx =
	std::find_if(vrtxCount.begin(), vrtxCount.end(),
		     std::bind2nd(std::not_equal_to<int>(), Dimension()));
    
    if (badVrtx	!= vrtxCount.end())
    {
	int badVrtxNum = badVrtx - vrtxCount.begin();
	std::ostringstream ost;
	ost << "Found vertex(" << badVrtxNum
	    << ") with wrong edge count. "
	    << "Was " << *badVrtx
	    << ", should have been " << Dimension() << ".";
	testassert(false, prefix+ost.str(), __FILE__, __LINE__);
    }
}

//***************
// checkSize(...)
//***************

template<class MTFactory>
void StructuredMeshCFAVTest<MTFactory>::checkSize(const VecSetInt &vecset,
						  int size,
						  const std::string &msg)
{
    const std::string prefix("checkSize");
    
    // Check to make sure that each set, in this vector of sets, is the
    // expected size.


    // The predicate of the following find_if algorithm does the following...
    //
    // returns true if a reference from the iteration's size() method
    // is not equal to the size input variable.
    //
    // This is done through composition and bondage:
    //
    //          compose(bind(operator!=(), size),
    //                  mem_fun(set<int>::size()))
    //
    //     which is equivalent to
    //
    //           iter->size() != size
    //
    // where iter is the current iterate in the vector,
    // and points to a set of ints.

    using std::find_if;
    using rtt_meshTest::compose1;
    using std::bind2nd;
    using std::not_equal_to;
    using std::mem_fun_ref;
    using std::set;
    
    VecSetInt::const_iterator badLoc =
	find_if(vecset.begin(), vecset.end(),
		compose1(bind2nd(not_equal_to<int>(), size),
			 mem_fun_ref(&set<int>::size)));

    if (badLoc != vecset.end())
    {
	int badNum = badLoc - vecset.begin();
	std::ostringstream ost;
	ost << "Found " << msg
	    << " at " << badNum
	    << " with wrong size.  "
	    << "Was " << badLoc->size()
	    << ", should have been " << size << ".";
	testassert(false, prefix+ost.str(), __FILE__, __LINE__);
    }
}

} // end namespace rtt_meshTest


//---------------------------------------------------------------------------//
//                              end of StructuredMeshCFAVTest.t.hh
//---------------------------------------------------------------------------//
