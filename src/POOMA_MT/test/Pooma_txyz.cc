//----------------------------------*-C++-*----------------------------------//
// Pooma_txyz.cc
// Julian C. Cummings
// Mon Oct 05 09:53:44 1998
//---------------------------------------------------------------------------//
// @>
//---------------------------------------------------------------------------//

// include files

#include "Utility/PoomaInfo.h"
#include "Utility/Inform.h"
#include "PETE/PoomaExpressions.h"
#include "POOMA_MT/PoomaMesh_XYZ.hh"

// main routine

int main(int argc, char* argv[])
{

  // initialize Pooma
  Pooma pooma(argc, argv);

  // create some Inform objects for tty output
  Inform msg0("txyz");
  Inform msgall("txyz",INFORM_ALL_NODES);

  msg0 << "Pooma initialized." << endl << endl;

  // hard-wire mesh type and parameters for now
  // UniformCartesian mesh: 4^3, parallel in each dimension
  int NumCells[3] = {4,4,4};
  double CellWidth[3] = {1.0,1.0,1.0};
  e_dim_tag Decomposition[3] = {PARALLEL,PARALLEL,PARALLEL};
  typedef UniformCartesian<3> PoomaMesh_t;
  typedef PoomaMesh_XYZ<PoomaMesh_t> MT_t;

  // create MT mesh object
  dsxx::SP<MT_t> spm = new MT_t(NumCells,CellWidth,Decomposition);

  // print out mesh sizes
  msg0 << "PoomaMesh constructed." << endl;
  msg0 << "ncx = " << spm->get_ncx() << ", ncy = "
       << spm->get_ncy() << ", ncz = " << spm->get_ncz() << endl << endl;

  // create some ccsf objects
  MT_t::ccsf x(spm), y(spm), z(spm);

  // print out field sizes
  msg0 << "MT::ccsf objects constructed." << endl << endl;
  msgall << "Size of MT::ccsf = " << x.size() << endl;

  // try some simple math with ccsf objects
  x = 1;
  y = 2;
  z = x + y;

  // iterate over elements and check correctness
  bool passed = true;
  MT_t::ccsf::const_iterator zit, zend = z.end();
  for (zit = z.begin(); zit != zend; ++zit)
    if (*zit != 3) passed = false;
  if (passed)
    msgall << "Simple arithmetic with MT::ccsf objects passed!" << endl;
  else
    msgall << "Simple arithmetic with MT::ccsf objects failed!" << endl;

  // test accessors for mesh info
  spm->get_xloc(x);
  spm->get_yloc(y);
  spm->get_zloc(z);

  // construct some fcdsf objects
  MT_t::fcdsf xf(spm), yf(spm), zf(spm);

  // try some simple math
  xf = 1;
  yf = xf;
  zf = xf - yf;

  // iterate over elements and check correctness
  passed = true;
  MT_t::fcdsf::const_iterator zfit, zfend = zf.end();
  for (zfit = zf.begin(); zfit != zfend; ++zfit)
    if (*zfit != 0) passed = false;
  if (passed)
    msgall << "Simple arithmetic with MT::fcdsf objects passed!" << endl;
  else
    msgall << "Simple arithmetic with MT::fcdsf objects failed!" << endl;

  // try a scatter assign
  MT_t::scatter(x,xf,OpAssign());

  msg0 << "Scatter to ccsf from fcdsf completed!" << endl << endl;

  // construct a bssf object
  MT_t::bssf b1(spm);

  // iterate through faces and assign
  b1.Uncompress();  // must uncompress before writing with iterator!
  MT_t::bssf::iterator b1it, b1end = b1.end();
  for (b1it = b1.begin(); b1it != b1end; ++b1it)
    *b1it = 0.0;

  msg0 << "Assignment to bssf completed!" << endl << endl;

  // try some gather and scatter operations
  MT_t::ccsf oneCC(spm);
  oneCC = 1.0;
  MT_t::fcdsf twoFC(spm);
  twoFC = 0.0;
  msg0 << "Before scatter-add to fcdsf from ccsf:" << endl;
  msgall << "oneCC = " << oneCC << endl;
  msgall << "twoFC = " << twoFC << endl;
  MT_t::scatter(twoFC, oneCC, OpAddAssign());
  msg0 << "After scatter-add to fcdsf from ccsf:" << endl;
  msgall << "oneCC = " << oneCC << endl;
  msgall << "twoFC = " << twoFC << endl;

  MT_t::ccsf threeCC(spm);
  threeCC = 3.0;
  MT_t::fcdsf nineFC(spm);
  nineFC = 1.0;
  msg0 << "Before scatter-multiply to fcdsf from ccsf:" << endl;
  msgall << "threeCC = " << threeCC << endl;
  msgall << "nineFC = " << nineFC << endl;
  MT_t::scatter(nineFC, threeCC, OpMultiplyAssign());
  msg0 << "After scatter-multiply to fcdsf from ccsf:" << endl;
  msgall << "threeCC = " << threeCC << endl;
  msgall << "nineFC = " << nineFC << endl;

  // apply boundary conditions using gather
  MT_t::gather(nineFC, b1, OpAssign());
  msg0 << "After applying zero boundary conditions:" << endl;
  msgall << "nineFC = " << nineFC << endl;

  // test some of the reduction functions
  double threeSum = MT_t::sum(threeCC);
  double nineMax = MT_t::max(nineFC);
  msg0 << "Sum of threeCC = " << threeSum << endl;
  msg0 << "Max of nineFC = " << nineMax << endl;

  Vektor<double,3> v1(3.0), v2(4.0);

  // double val = rtt_traits::vector_traits<Vektor<double,3> >::dot(v1, v2);
  double val = dot(v1, v2);

  msg0 << "dot of (3,3,3) and (4,4,4) is " << val << endl;

  return 0;
}

//---------------------------------------------------------------------------//
//                              end of Pooma_txyz.cc
//---------------------------------------------------------------------------//
