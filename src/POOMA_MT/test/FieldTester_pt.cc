//----------------------------------*-C++-*----------------------------------//
// FieldTester_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTFields.hh"
#include "meshTest/FieldTester.t.hh"

using rtt_meshTest::TestMTFields;
using rtt_meshTest::FieldTester;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

typedef TestMTFields<PoomaMesh_XYZFactory> Tester;

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::run();

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::t1();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::t2();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::t3();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::t4();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::t5();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::t6();

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::t1();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::t2();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::t3();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::t4();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::t5();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::t6();

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::t1();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::t2();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::t3();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::t4();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::t5();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::t6();

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::t1();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::t2();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::t3();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::t4();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::t5();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::t6();

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::t1();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::t2();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::t3();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::t4();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::t5();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::t6();

using rtt_meshTest::FieldTesters::f1;
using rtt_meshTest::FieldTesters::f2;
using rtt_meshTest::FieldTesters::f3;
using rtt_meshTest::FieldTesters::f4;
using rtt_meshTest::FieldTesters::f5;

#define F1(T1) \
template bool f1<T1>(const T1 &, const T1 &);

#define F2(T1) \
template bool f2<T1>(const T1::iterator &);

#define F3(T1) \
template bool f3<T1>(const T1::iterator &, const T1::iterator &);

#define F4(T1) \
template bool f4<T1>(const T1::const_iterator &);

#define F5(T1) \
template bool f5<T1>(const T1::const_iterator &, const T1::const_iterator &);

#define FNField(T1) \
F1(T1); \
F2(T1); \
F3(T1); \
F4(T1); \
F5(T1);

FNField(Tester::BSTF::XD);
FNField(Tester::CCTF::XD);
FNField(Tester::FCDTF::XD);
FNField(Tester::NCTF::XD);
FNField(Tester::VCTF::XD);

//---------------------------------------------------------------------------//
//                              end of FieldTester_pt.cc
//---------------------------------------------------------------------------//
