//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTConnFacesArroundVrtx.hh
 * \author Randy M. Roberts
 * \date   Fri Sep 10 12:55:42 1999
 * \brief  Header file for the TestMTConnFacesArroundVrtx class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTConnFacesArroundVrtx_hh__
#define __meshTest_TestMTConnFacesArroundVrtx_hh__

#include <string>
#include <iosfwd>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class TestMTConnFacesArroundVrtx
 *
 * \brief Main interface to testing the MT::ConnFacesArroundVerices service.
 *
 * The TestMTConnFacesArroundVrtx class is templated on an MTFactory.
 * The MTFactory is used to create coupled instances of meshes and field
 * constructors, along with typedefs that determine whether the mesh is
 * structured or unstructured.
 *
 * Required MTFactory services include:
 *
 * MTFactory::MT -- The mesh type
 *
 * MTFactory::FieldConstructor -- The field constructor type.
 *
 * MTFactory::Product -- A class that is returned by the MTFactroy::create()
 * method.  This class is responsible for returning references to the mesh
 * and field constructor via the mesh() and fieldConstructor() methods,
 * respectively.
 *
 * MTFactory::Structured -- A tag that, if used, implies that the mesh
 * is a structured mesh.
 *
 * MTFactory::UnStructured -- A tag that, if used, implies that the mesh
 * is an unstructured mesh.
 *
 * MTFactory::Structuring -- A typedef to either, Structured, or UnStructured,
 * used to determine whether the mesh is actually one or the other.
 *
 * MTFactory::Product MTFactory::create() -- The method that returns a new
 * pair of meshes and field constructors.
 */
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMTConnFacesArroundVrtx 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

    typedef typename MTFactory::Structured Structured;
    typedef typename MTFactory::UnStructured UnStructured;
    typedef typename MTFactory::Structuring Structuring;
    
    // DATA
    
    MTFactory &meshFactory_m;
    
    std::ostream &os_m;
    
    bool passed_m;

  public:

    // CREATORS

    //! Constructor
    
    TestMTConnFacesArroundVrtx(MTFactory &meshFactory_in, std::ostream &os_in)
	: meshFactory_m(meshFactory_in), os_m(os_in),
	  passed_m(false)
    {
	/* empty */
    }

    //! Destructor
    
    ~TestMTConnFacesArroundVrtx() { /* empty */ }

    // MANIPULATORS

    //! Main interface to testing class.
    
    void run();

    // ACCESSORS

    //! Returns success of previously ran run() method.
    
    bool passed() const { return passed_m; }
    
  private:
    
    // DISSALLOWED CREATORS

    TestMTConnFacesArroundVrtx(const TestMTConnFacesArroundVrtx &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTConnFacesArroundVrtx& operator=(const TestMTConnFacesArroundVrtx &rhs);

    // IMPLEMENTATION

    void error(bool &passed, const std::string &msg);

    void run(Structured);

    void run(UnStructured);

};

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTConnFacesArroundVrtx_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTConnFacesArroundVrtx.hh
//---------------------------------------------------------------------------//
