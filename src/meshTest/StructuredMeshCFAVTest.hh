//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/StructuredMeshCFAVTest.hh
 * \author Randy M. Roberts
 * \date   Wed Sep 22 13:23:07 1999
 * \brief  Header file for the base class StructuredMeshCFAVTest
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_StructuredMeshCFAVTest_hh__
#define __meshTest_StructuredMeshCFAVTest_hh__

#include "Tester.hh"
#include "Cell.hh"

#include <vector>
#include <set>
#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class StructuredMeshCFAVTest
 * \brief Class used by TestMTConnFacesArroundVrtx to test
 * the structured MT::ConnFacesAroundVertices service.
 *
 * The StructuredMeshCFAVTest class is used by TestMTConnFacesArroundVrtx
 * to run tests for a
 * structured Mesh MT.  These tests excersize the
 * MT::ConnFacesAroundVertices service.
 *
 */
// 
//===========================================================================//
template<class MTFactory>
class StructuredMeshCFAVTest : public Tester
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;
    
    typedef std::vector< std::set<int> > VecSetInt;

  private:
    
    // DATA

    Tester &parent_m;
    
    // The order of the following members are important to
    // the constructor.
    
    MTFactory &meshFactory_m;
    MTFactoryProduct meshProduct_m;
    MT &mesh_m;
    FieldConstructor &fCtor_m;
    typename MT::vcif vindices_m;
    typename MT::fcdif findices_m;

  public:

    // STATIC METHODS

    static int Dimension()
    {
	return MTFactory::Dimension();
    }

    static int NumVerticesPerCell()
    {
	// 2 ^ Dim
	return 1 << Dimension();
    }

    static int NumFacesPerCell()
    {
	// 2 * Dim
	return 2 * Dimension();
    }

    static int NumEdgesPerCell()
    {
	// Dim * 2 ^ (Dim -1)
	return Dimension() * (1 << (Dimension()-1));
    }
    
    static int NumFacesPerVertex()
    {
	return Dimension();
    }

    static int NumVerticesPerFace()
    {
	return 1 << (Dimension()-1);
    }

    // CREATORS
    
    //! Constructor

    StructuredMeshCFAVTest(Tester &parent_in, MTFactory &meshFactory_in);
    
    //! Destructor
    
    ~StructuredMeshCFAVTest() { /* empty */ }

    // MANIPULATORS
    
    //! The run method is the primary interface into the class.

    void run();

    // ACCESSORS

  protected:

    // PROTECTED MANIPULATORS

    virtual void testassert(bool passed, const std::string &msg)
    {
	parent_m.testassert(passed, Name()+": "+msg);
	if (!passed)
	    setPassed(false);
    }
    
    virtual void testassert(bool passed, const std::string &msg,
			    const std::string &file,
			    const std::string &line)
    {
	Tester::testassert(passed, msg, file, line);
    }
    
    virtual void testassert(bool passed, const std::string &msg,
			    const std::string &file, int line)
    {
	Tester::testassert(passed, msg, file,line);
    }
    
  private:

    // DISSALLOWED CREATORS

    StructuredMeshCFAVTest(const StructuredMeshCFAVTest &rhs);

    // DISSALLOWED MANIPULATORS
    
    StructuredMeshCFAVTest& operator=(const StructuredMeshCFAVTest &rhs);

    // IMPLEMENTATION

    int get_ncells() const { return mesh_m.get_ncells(); }

    void t1();

    template<class CT, class CIT, class CVIT>
    void t2(const std::string &constness);

    template<class CT, class CIT, class CVIT>
    void t3(const std::string &constness);

    template<class CT, class CIT, class CVIT>
    void t4(const std::string &constness);

    typename MT::vcif &vindices() { return vindices_m; }

    typename MT::fcdif &findices() { return findices_m; }
    
    void checkConnectivity(const VecSetInt &facesArroundVertices,
			   const VecSetInt &verticesArroundFaces);

    std::set<Cell> getEdgeSet(const VecSetInt &verticesArroundFaces);

    std::vector<int> findIntersection(const std::set<int> &set1,
				   const std::set<int> &set2);

    void checkEdgeSetforVertexCount(const std::set<Cell> &edgeSet);

    void checkSize(const VecSetInt &vecset, int size, const std::string &msg);
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_StructuredMeshCFAVTest_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/StructuredMeshCFAVTest.hh
//---------------------------------------------------------------------------//
