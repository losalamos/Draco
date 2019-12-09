//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstDracoStrings.cc
 * \author Kelly G. Thompson <kgt@lanl.gov>
 * \date   Wednesday, Aug 23, 2017, 13:30 pm
 * \brief  Test functions defined in ds++/DracoStrings.hh
 * \note   Copyright (C) 2017-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_trim(UnitTest &ut) {

  cout << "Begin test_trim checks...\n";

  string const case1("   This is a string.     ");
  FAIL_IF_NOT(trim(case1) == string("This is a string."));

  string const case2("This is a string.");
  FAIL_IF_NOT(trim(case2) == string("This is a string."));

  string const case3("#  This is a string.  ");
  FAIL_IF_NOT(trim(case3, string("# ")) == string("This is a string."));

  if (ut.numFails == 0)
    PASSMSG("test_trim: All tests pass.");
  else
    FAILMSG("test_trim: FAILED");

  return;
}

//---------------------------------------------------------------------------//
void test_prune(UnitTest &ut) {

  cout << "\nBegin test_prune checks...\n";
  unsigned const nf = ut.numFails;

  string const case2("This is a string.");
  FAIL_IF_NOT(prune(case2, string(" ")) == string("Thisisastring."));

  string const numbers("123456789");
  string const case3("123 hello 456 world 789!");
  FAIL_IF_NOT(trim(prune(case3, numbers)) == string("hello  world !"));

  if (ut.numFails == nf)
    PASSMSG("test_prune: All tests pass.");
  else
    FAILMSG("test_prune: FAILED");

  return;
}

//---------------------------------------------------------------------------//
void test_tokenize(UnitTest &ut) {

  cout << "\nBegin test_tokenize checks...\n";
  unsigned const nf = ut.numFails;

  string const case1("1, 2, 3, 4, 5, 6");
  vector<string> case1words = tokenize(case1, ",");
  FAIL_IF_NOT(case1words.size() == 6);
  FAIL_IF_NOT(case1words[3] == string("4"));

  string const case2("This is a string.");
  vector<string> case2words = tokenize(case2);
  FAIL_IF_NOT(case2words.size() == 4);
  FAIL_IF_NOT(case2words[0] == string("This"));

  if (ut.numFails == nf)
    PASSMSG("test_tokenize: All tests pass.");
  else
    FAILMSG("test_tokenize: FAILED");

  return;
}

//---------------------------------------------------------------------------//
void test_parse_number(UnitTest &ut) {

  cout << "\nBegin test_parse_number checks...\n";
  unsigned const nf = ut.numFails;

  string const case1("1");
  string const case3("3.1415");
  string const case4("four");

  double const deps = numeric_limits<double>::epsilon();
  float const feps = numeric_limits<float>::epsilon();

  FAIL_IF_NOT(parse_number<int>(case1) == 1);
  FAIL_IF_NOT(parse_number<long>(case1) == 1l);
  FAIL_IF_NOT(parse_number<unsigned long>(case1) == 1ul);
  FAIL_IF_NOT(parse_number<int32_t>(case1) == 1);
  FAIL_IF_NOT(parse_number<int64_t>(case1) == 1l);
  FAIL_IF_NOT(parse_number<uint32_t>(case1) == 1u);
  FAIL_IF_NOT(parse_number<uint64_t>(case1) == 1ul);
  FAIL_IF_NOT(soft_equiv(parse_number<float>(case1), 1.0f, feps));
  FAIL_IF_NOT(soft_equiv(parse_number<double>(case1), 1.0, deps));

  FAIL_IF_NOT(parse_number<int>(case3) == 3);
  FAIL_IF_NOT(soft_equiv(parse_number<float>(case3), 3.1415f, feps));
  FAIL_IF_NOT(soft_equiv(parse_number<double>(case3), 3.1415, deps));

  try {
    parse_number<int>(case4);
    FAILMSG("Failed to capture invalid_argument given to stoi.");
  } catch (std::invalid_argument & /*error*/) {
    PASSMSG("invalid_argument givent to stox(str).");
  } catch (...) {
    FAILMSG("Failed to capture invalid_argument given to stoi.");
  }

  try {
    parse_number<int>(string("5000000000000"));
    FAILMSG("Failed to capture out_of_range given to stoi.");
  } catch (std::out_of_range & /*error*/) {
    PASSMSG("invalid_argument givent to stox(str).");
  } catch (...) {
    FAILMSG("Failed to capture out_of_range given to stoi.");
  }

  if (ut.numFails == nf)
    PASSMSG("test_parse_number: All tests pass.");
  else
    FAILMSG("test_parse_number: FAILED");

  return;
}

//---------------------------------------------------------------------------//
void test_string_to_numvec(UnitTest &ut) {

  cout << "\nBegin test_string_to_numvec checks...\n";
  unsigned const nf = ut.numFails;

  string const case1("{ 1, 2, 3 }");
  string const case2("{ 1.1, 2.2, 3.3 }");
  string const case3("[ 1.1 2.2 3.3 )");

  vector<int> const case1v = string_to_numvec<int>(case1);
  vector<double> const case2v = string_to_numvec<double>(case2);
  vector<double> const case3v =
      string_to_numvec<double>(case3, string("[)"), string(" "));

  vector<int> const case1ref = {1, 2, 3};
  vector<double> const case2ref = {1.1, 2.2, 3.3};

  FAIL_IF_NOT(case1v == case1ref);
  FAIL_IF_NOT(soft_equiv(case2v, case2ref));
  FAIL_IF_NOT(soft_equiv(case3v, case2ref));

  if (ut.numFails == nf)
    PASSMSG("test_string_to_numvec: All tests pass.");
  else
    FAILMSG("test_string_to_numvec: FAILED");

  return;
}

//----------------------------------------------------------------------------//
void test_tostring(UnitTest &ut) {

  cout << "\nBegin test_tostring checks...\n";
  unsigned const nf = ut.numFails;

  string const pi = rtt_dsxx::to_string(M_PI, 20);

  if (soft_equiv(M_PI, atof(pi.c_str())))
    PASSMSG("pi correctly written/read");
  else
    FAILMSG("pi NOT correctly written/read");

  double const foo(2.11111111);
  unsigned int const p(23);
  // Must be careful to use rtt_dsxx::to_string and avoid std::to_string --
  // especially after 'using namespace std.'
  string s1(rtt_dsxx::to_string(foo, p));
  string s2(rtt_dsxx::to_string(foo));
  if (s1 == s2)
    PASSMSG("double printed using default formatting.");
  else
    FAILMSG("double printed with wrong format!");

  if (ut.numFails == nf)
    PASSMSG("test_tostring: All tests pass.");
  else
    FAILMSG("test_tostring: FAILED");

  return;
}

//----------------------------------------------------------------------------//
void test_upper_lower(UnitTest &ut) {

  cout << "\nBegin test_upper_lower checks...\n";
  unsigned const nf = ut.numFails;

  std::string const mixedCase("This StRiNg HAS mixed CAse.");
  std::string const upperCase = rtt_dsxx::string_toupper(mixedCase);
  std::string const lowerCase = rtt_dsxx::string_tolower(mixedCase);

  FAIL_IF_NOT(upperCase == "THIS STRING HAS MIXED CASE.");
  FAIL_IF_NOT(lowerCase == "this string has mixed case.");

  if (ut.numFails == nf)
    PASSMSG("test_upper_lower: All tests pass.");
  else
    FAILMSG("test_upper_lower: FAILED");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_trim(ut);
    test_prune(ut);
    test_tokenize(ut);
    test_parse_number(ut);
    test_string_to_numvec(ut);
    test_tostring(ut);
    test_upper_lower(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstDracoStrings.cc
//---------------------------------------------------------------------------//
