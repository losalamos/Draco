//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstAbstract_Class_Parser.cc
 * \author Kent G. Budge
 * \date   Tue Nov  9 14:34:11 2010
 * \brief  Test the Abstract_Class_Parser template
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "../Class_Parser.hh"
#include "../Abstract_Class_Parser.hh"
#include "../utilities.hh"
#include "../File_Token_Stream.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/path.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*
 * The following would typically be declared in a file associated with the
 * Parent class.
 */

class Parent
{
  public:
    virtual ~Parent(){}
    virtual string name() = 0;
    static SP<Parent> parse(Token_Stream &tokens);
    static void register_model(string const &keyword,
                               SP<Parent> parse_function(Token_Stream&) );

};

SP<Parent> parsed_child;

//---------------------------------------------------------------------------//
namespace rtt_parser
{

template<>
rtt_parser::Parse_Table
Class_Parser<Parent, Parent, true>::parse_table_(NULL, 0);

template<>
void Class_Parser<Parent, Parent, true>::post_sentinels_()
{
    parsed_child.reset();
}

template<>
void Class_Parser<Parent, Parent, true>::
check_completeness_(Token_Stream &tokens)
{
    if (parsed_child == SP<Parent>())
    {
        tokens.report_semantic_error("no parent specified");
    }
}

template<>
SP<Parent>
Class_Parser<Parent, Parent, true>::create_object_()
{
    return parsed_child;
}

} // namespace rtt_parser
//---------------------------------------------------------------------------//

Parse_Table &get_parse_table()
{
    return Class_Parser<Parent, Parent, true>::parse_table_;
}

SP<Parent> &get_parsed_object()
{
    return parsed_child;
}

SP<Parent> Parent::parse(Token_Stream &tokens)
{
    return Class_Parser<Parent, Parent, true>::parse(tokens);
}

void Parent::
register_model(string const &keyword,
               SP<Parent> parse_function(Token_Stream&) )
{
    Abstract_Class_Parser<Parent, get_parse_table, get_parsed_object>::
        register_child(keyword, parse_function);
}

//---------------------------------------------------------------------------//
/*
 * The following is what you would expect to find in a file associated with
 * the Son class.
 */

class Son : public Parent
{
  public:

    virtual string name(){ return "son"; }

    Son(double /*snip_and_snails*/){}

    static SP<Son> parse(Token_Stream &tokens);

    static SP<Parent> parse_parent(Token_Stream &tokens)
    {
        return parse(tokens);
    }
};

SP<Parent> parse_son(Token_Stream &tokens)
{
    return Son::parse(tokens);
}

double parsed_snips_and_snails;
 
void parse_snips_and_snails(Token_Stream &tokens, int)
{
    if (parsed_snips_and_snails>=0.0)
    {
        tokens.report_semantic_error("snips and snails already specified");
    }

    parsed_snips_and_snails = parse_real(tokens);
    if (parsed_snips_and_snails<0.0)
    {
        tokens.report_semantic_error("snips and snails must not be "
                                     "negative");
        parsed_snips_and_snails = 2;
    }
}

const Keyword son_keywords[] =
{
    {"snips and snails", parse_snips_and_snails, 0, ""},
};

const unsigned number_of_son_keywords = sizeof(son_keywords)/sizeof(Keyword);


namespace rtt_parser
{
template<>
Parse_Table
Class_Parser<Son>::
parse_table_(son_keywords,
             number_of_son_keywords);

template<>
void Class_Parser<Son>::post_sentinels_()
{
    parsed_snips_and_snails = -1;
}

template<>
void
Class_Parser<Son>::
check_completeness_(Token_Stream &tokens)
{
    if (parsed_snips_and_snails == -1) 
    {
        tokens.report_semantic_error("no snips and snails specified");
    }
}

template<>
SP<Son>
Class_Parser<Son>::create_object_()
{
    SP<Son> Result(new Son(parsed_snips_and_snails));
    
    return Result;
}

} // end namespace rtt_parser

SP<Son> Son::parse(Token_Stream &tokens)
{
    return Class_Parser<Son>::parse(tokens);
}


//---------------------------------------------------------------------------//
/*
 * The following is what you would expect to find in a file associated with
 * the Daughter class.
 */

class Daughter : public Parent
{
  public:

    virtual string name(){ return "daughter"; }

    Daughter(double /*sugar_and_spice*/){}

    static SP<Daughter> parse(Token_Stream &tokens);

    static SP<Parent> parse_parent(Token_Stream &tokens)
    {
        return parse(tokens);
    }
};

SP<Parent> parse_daughter(Token_Stream &tokens)
{
    return Daughter::parse(tokens);
}

double parsed_sugar_and_spice;
 
void parse_sugar_and_spice(Token_Stream &tokens, int)
{
    if (parsed_sugar_and_spice>=0.0)
    {
        tokens.report_semantic_error("sugar and spice already specified");
    }

    parsed_sugar_and_spice = parse_real(tokens);
    if (parsed_sugar_and_spice<0.0)
    {
        tokens.report_semantic_error("sugar and spice must not be "
                                     "negative");
        parsed_sugar_and_spice = 2;
    }
}

const Keyword daughter_keywords[] =
{
    {"sugar and spice", parse_sugar_and_spice, 0, ""},
};

const unsigned number_of_daughter_keywords = sizeof(daughter_keywords)/sizeof(Keyword);


namespace rtt_parser
{
template<>
Parse_Table
Class_Parser<Daughter>::
parse_table_(daughter_keywords,
             number_of_daughter_keywords);

template<>
void Class_Parser<Daughter>::post_sentinels_()
{
    parsed_sugar_and_spice = -1;
}

template<>
void
Class_Parser<Daughter>::
check_completeness_(Token_Stream &tokens)
{
    if (parsed_sugar_and_spice == -1) 
    {
        tokens.report_semantic_error("no sugar and spice specified");
    }
}

template<>
SP<Daughter>
Class_Parser<Daughter>::create_object_()
{
    SP<Daughter> Result(new Daughter(parsed_sugar_and_spice));
    
    return Result;
}

} // end namespace rtt_parser

SP<Daughter> Daughter::parse(Token_Stream &tokens)
{
    return Class_Parser<Daughter>::parse(tokens);
}

/* the followingn would typically live in some client file for Parent. */

SP<Parent> parent;

static void parse_parent(Token_Stream &tokens, int)
{
    parent = Parent::parse(tokens);
}

const Keyword top_keywords[] =
{
    {"parent", parse_parent, 0, ""},
};

const unsigned number_of_top_keywords = sizeof(top_keywords)/sizeof(Keyword);

Parse_Table top_parse_table(top_keywords, number_of_top_keywords);


//---------------------------------------------------------------------------//

/* Test all the above */

void test(UnitTest &ut)
{
    Parent::register_model("son",
                           Son::parse_parent );

    Parent::register_model("daughter",
                           Daughter::parse_parent );

    // Build path for the input file
    string const sadInputFile(ut.getTestInputPath()
                           + std::string("sons_and_daughters.inp") );

    File_Token_Stream tokens( sadInputFile );

    top_parse_table.parse(tokens);

    cout << parent->name() << endl;

    if (tokens.error_count()==0 &&
        parent!=SP<Parent>() &&
        parent->name()=="son")
    {
        ut.passes("Parsed son correctly");
    }
    else
    {
        ut.failure("Did NOT parse son correctly");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try { test(ut); }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstAbstract_Class_Parser.cc
//---------------------------------------------------------------------------//
