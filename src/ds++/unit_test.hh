//----------------------------------*-C++-*----------------------------------//
/*!
  \file   unit_test.hh
  \author Rob Lowrie
  \date   2001-05-20
  \brief  Testing utility for unit tests.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_unit_test_hh
#define rtt_dsxx_unit_test_hh

#include <iostream>

// A useful macro

#define UNIT_TEST(x) rtt_dsxx::unit_test(x, __LINE__, __FILE__)

namespace rtt_dsxx {

void unit_test(const bool pass, int line, char *file);

//---------------------------------------------------------------------------//
/*!
  \brief  Prints passed or failed, suitable for parsing
          by tools/test_filter.py.

  \param pass   If true, print pass. Otherwise print fail.
  \param line   Line number of calling routine.  Typically, pass
                __LINE__.
  \param file   File name containing calling routine.  Typically,
                pass __FILE__.
*/
//---------------------------------------------------------------------------//
void unit_test(const bool pass, int line, char *file)
{
    if ( pass ) {
	std::cout << " test: pass\n";
    }
    else {
	std::cout << " test: fail"
		  << "  (line " << line << " of " << file << ")\n";
    }
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_unit_test_hh

//---------------------------------------------------------------------------//
// end of ds++/unit_test.hh
//---------------------------------------------------------------------------//
