//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstClass_Parser.cc
 * \author Kent Budge
 * \date   Mon Aug 28 07:36:50 2006
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "../String_Token_Stream.hh"
#include "../utilities.hh"

#include "../Class_Parser.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
class DummyClass
{
  public:
    DummyClass(double const insouciance)
        :
        insouciance(insouciance)
    {
        Require(insouciance>=0.0);
    }

    double Get_Insouciance() const { return insouciance; }

  private:

    double insouciance;
};

namespace
{
//---------------------------------------------------------------------------//
double parsed_insouciance;

//---------------------------------------------------------------------------//
void Parse_Insouciance(Token_Stream &tokens, int)
{
    if (parsed_insouciance>=0.0)
    {
        tokens.report_semantic_error("duplicate specification of insouciance");
    }
    parsed_insouciance = parse_real(tokens);
    if (parsed_insouciance<0)
    {
        tokens.report_semantic_error("insouciance must be nonnegative");
        parsed_insouciance = 1;
    }
}

//---------------------------------------------------------------------------//
const Keyword dummy_keywords[] =
{
    {"insouciance",           Parse_Insouciance, 0, ""},
};
const unsigned number_of_dummy_keywords =
sizeof(dummy_keywords)/sizeof(Keyword);

} // anonymous namespace

namespace rtt_parser
{
//---------------------------------------------------------------------------//
template<>
Parse_Table
Class_Parser<DummyClass>::parse_table_(dummy_keywords,
                                      number_of_dummy_keywords);

//---------------------------------------------------------------------------//
template<>
void Class_Parser<DummyClass>::post_sentinels_()
{
    parsed_insouciance = -1.0;  // sentinel value
}

//---------------------------------------------------------------------------//
template<>
void Class_Parser<DummyClass>::check_completeness_(Token_Stream &tokens)
{
    if (parsed_insouciance<0)
    {
        tokens.report_semantic_error("insouciance was not specified");
    }
}

//---------------------------------------------------------------------------//
template<>
SP<DummyClass> Class_Parser<DummyClass>::create_object_()
{
    return SP<DummyClass>(new DummyClass(parsed_insouciance));
}

}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstClass_Parser( UnitTest & ut)
{
    string text = "insouciance = 3.3\nend\n";
    String_Token_Stream tokens(text);
    
    SP<DummyClass> dummy = Class_Parser<DummyClass>::parse(tokens);

    if (dummy != SP<DummyClass>())
    {
        ut.passes("parsed the class object");

        if (dummy->Get_Insouciance() == 3.3)
        {
            ut.passes("parsed the insouciance correctly");
        }
        else
        {
            ut.failure("did NOT parse the insouciance correctly");
        }
    }
    else
    {
        cout << tokens.messages() << endl;
        ut.failure("did NOT parse the class object");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut(argc, argv, release);
        tstClass_Parser(ut);
        ut.status();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstClass_Parser, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstClass_Parser, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstClass_Parser.cc
//---------------------------------------------------------------------------//
