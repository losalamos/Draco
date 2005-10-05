//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstToken.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "parser_test.hh"
#include "../Release.hh"
#include "../Token.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

using namespace std;

using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void token_test()
{
    Token end_token  (END,               "parser_test 1");
    Token other_token('$',               "parser_test 2");
    Token real_token (REAL,  "+1.56e-3", "parser_test 3");

    if (end_token.Type()!=END) ITFAILS;
    if (other_token.Type()!=OTHER) ITFAILS;
    if (real_token.Type()!=REAL) ITFAILS;

    if (other_token.Text()!="$") ITFAILS;
    if (real_token.Text()!="+1.56e-3") ITFAILS;

    if (end_token.Location()!="parser_test 1") ITFAILS;
    if (other_token.Location()!="parser_test 2") ITFAILS;
    if (real_token.Location()!="parser_test 3") ITFAILS;

    if (!Is_Integer_Text("057133")) ITFAILS;
    if (Is_Integer_Text("08223")) ITFAILS;
    if (!Is_Integer_Text("663323")) ITFAILS;
    if (!Is_Integer_Text("0x33a8")) ITFAILS;
    if (!Is_Keyword_Text("This is a test")) ITFAILS;
    if (Is_Keyword_Text("This is 1 test")) ITFAILS;
    if (!Is_Real_Text("+1.56e-3")) ITFAILS;
    if (Is_Real_Text("1.39d-3")) ITFAILS;
    if (!Is_String_Text("\"This is a test.\"")) ITFAILS;
    if (Is_String_Text("\"This is a test")) ITFAILS;

    if (real_token == end_token)
    {
	FAILMSG("unlike token equality test did NOT return false");
    }
    else
    {
	PASSMSG("unlike token equality test returned false");
    }

    if (real_token == real_token)
    {
	PASSMSG("like token equality test did returned true");
    }
    else
    {
	FAILMSG("like token equality test did NOT return true");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_parser::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	if ( rtt_c4::nodes() == 1 ) token_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstToken, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_parser_test::passed) 
    {
        cout << "**** tstToken Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    rtt_c4::global_barrier();

    cout << "Done testing tstToken." << endl;

    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstToken.cc
//---------------------------------------------------------------------------//
