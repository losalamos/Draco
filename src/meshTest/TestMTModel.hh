//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTModel.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Fri Aug 20 09:11:42 1999
 * \brief  Header file for the TestMTModel class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTModel_hh__
#define __meshTest_TestMTModel_hh__

#include "Tester.hh"
#include <iosfwd>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class TestMTModel
 *
 * \brief  This class tests whether a given MT is a model of the
 *            "Solon" MT concept.
 *
 * The TestMTModel class is templated on an MTFactory.
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
class TestMTModel : public Tester
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

    // DATA

    MTFactory &meshFactory_m;
    
  public:

    // CREATORS

    //! Constructor
    
    TestMTModel(MTFactory &meshFactory_in, std::ostream &os_in)
	: Tester("TestMTModel", os_in), meshFactory_m(meshFactory_in)
    {
	/* empty */
    }

    //! Destructor
 
    ~TestMTModel() { /* empty */ }

    // MANIPULATORS

    //! Main interface to testing class.
    
    void run();

    // ACCESSORS
    
  private:
    
    // DISSALLOWED CREATORS

    TestMTModel(const TestMTModel &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTModel& operator=(const TestMTModel &rhs);

    // IMPLEMENTATION

    void t1();
    void t2();
    
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTModel_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTModel.hh
//---------------------------------------------------------------------------//
