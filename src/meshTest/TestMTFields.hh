//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTFields.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Mon Aug 23 16:49:16 1999
 * \brief  Header file for the TestMTFields class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTFields_hh__
#define __meshTest_TestMTFields_hh__

#include "FieldGroup.hh"
#include "DoubleContainer.hh"

#include <iosfwd>
#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class TestMTFields
 * \brief Main interface to testing the nested field class services.
 *
 * The TestMTFields class is templated on an MTFactory.
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
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMTFields 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

  public:
    
    enum FieldGroupID {ebstf, ecctf, efcdtf, enctf, evctf};

    inline static const char *Name(int id);
    
  public:

    //! This typedef is in lieu of compiler support for template template
    //! arguments.
    
    typedef FieldGroup<typename MT::template bstf<double>,
	typename MT::template bstf<int>, typename MT::template bstf<long>,
	typename MT::template bstf <DoubleContainer>, ebstf > BSTF;
    
    //! This typedef is in lieu of compiler support for template template
    //! arguments.
    
    typedef FieldGroup<typename MT::template cctf<double>,
	typename MT::template cctf<int>, typename MT::template cctf<long>,
	typename MT::template cctf<DoubleContainer>, ecctf > CCTF;
    
    //! This typedef is in lieu of compiler support for template template
    //! arguments.
    
    typedef FieldGroup<typename MT::template fcdtf<double>,
	typename MT::template fcdtf<int>, typename MT::template fcdtf<long>,
	typename MT::template fcdtf<DoubleContainer>, efcdtf > FCDTF;
    
    //! This typedef is in lieu of compiler support for template template
    //! arguments.
    
    typedef FieldGroup<typename MT::template nctf<double>,
	typename MT::template nctf<int>, typename MT::template nctf<long>,
	typename MT::template nctf<DoubleContainer>, enctf > NCTF;
    
    //! This typedef is in lieu of compiler support for template template
    //! arguments.
    
    typedef FieldGroup<typename MT::template vctf<double>,
	typename MT::template vctf<int>, typename MT::template vctf<long>,
	typename MT::template vctf<DoubleContainer>, evctf > VCTF;
    
    // DATA

  private:

    MTFactory &meshFactory_m;
    
    std::ostream &os_m;
    
    bool passed_m;

  public:

    // CREATORS
    
    //! Constructor
    
    TestMTFields(MTFactory &meshFactory_in, std::ostream &os_in)
	: meshFactory_m(meshFactory_in), os_m(os_in),
	  passed_m(false)
    {
	/* empty */
    }

    //! Destructor
    
    ~TestMTFields() { /* empty */ }

    // MANIPULATORS

    //! Primary interface into this class.
    
    template<class FGRP>
    void run();
    
    // ACCESSORS

    //! Returns success of previously ran run() method.
    
    bool passed() const { return passed_m; }

  private:
    
    // DISSALLOWED CREATORS

    TestMTFields(const TestMTFields &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTFields& operator=(const TestMTFields &rhs);

    // IMPLEMENTATION
};

template<class MTFactory>
const char *TestMTFields<MTFactory>::Name(int id)
{
    switch (id)
    {
    case ebstf:
	return "bstf";
    case ecctf:
	return "cctf";
    case efcdtf:
	return "fcdtf";
    case enctf:
	return "nctf";
    case evctf:
	return "vctf";
    }
    return "unknown type";
}

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTFields_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTFields.hh
//---------------------------------------------------------------------------//
