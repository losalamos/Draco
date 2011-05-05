//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstExpression.cc
 * \author Kent Budge
 * \date   Wed Jul 26 08:15:18 2006
 * \brief  Test the Expression class and expression parsing.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <cmath>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"
#include "../String_Token_Stream.hh"
#include "ds++/square.hh"

#include "../Expression.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstExpression(UnitTest &ut)
{
    // Create an expression as a String_Token_Stream.

    string const expression_text =
        "(((+1 && 1.3)||!(y<-m))/5+(2>1)*(r/m)*(2.7-1.1*(z/m))^2)*(t/s)";
    
    String_Token_Stream tokens(expression_text);

    typedef pair<unsigned, Unit> vd;
    map<string, vd > variable_map;

    variable_map["r"] = vd(0, m);
    variable_map["y"] = vd(1, m);
    variable_map["z"] = vd(2, m);
    variable_map["t"] = vd(3, s);

    vector<string> vars(4);
    vars[0] = "r";
    vars[1] = "y";
    vars[2] = "z";
    vars[3] = "t";
    
    SP<Expression const> expression = Expression::parse(4,
                                                        variable_map,
                                                        tokens);

    if (tokens.error_count()==0 && tokens.lookahead().type()==EXIT)
    {
        ut.passes("expression successfully parsed");
    }
    else
    {
        ut.failure("expression NOT successfully parsed");
        cerr << tokens.messages() << endl;
    }

    ostringstream expression_text_copy;
    expression->write(vars, expression_text_copy);

    char const * expression_text_raw =
        "((1&&1.3||!(y<-m))/5+(2>1)*r/m*pow(2.7-1.1*z/m,2))*t/s";
    // changes slightly due to stripping of extraneous whitespace,
    // parentheses, and positive prefix
    if (expression_text_copy.str()==expression_text_raw)
    {
        ut.passes("expression successfully rendered as text");
    }
    else
    {
        ut.failure("expression NOT successfully rendered as text");
        cerr << expression_text_raw << endl;
        cerr << expression_text_copy.str() << endl;
    }

#if 0
    // Test calculus
    
    SP<Expression const> deriv = expression->pd(3);

    ostringstream deriv_text;
    deriv->write(vars, deriv_text);

    expression_text_raw =
        "((1&&1.3||!(y<-m))/5+(2>1)*r/m*pow(2.7-1.1*z/m,2))*t/s";

    if (deriv_text.str()==expression_text_raw)
    {
        ut.passes("expression successfully derived");
    }
    else
    {
        ut.failure("expression NOT successfully derived");
        cerr << expression_text_raw << endl;
        cerr << deriv_text.str() << endl;
    }
#endif

    vector<double> xs;
    
    double x = 1.2;
    double r = x;
    double y = 3.1;
    double z = 0.0;
    double t = 2.8;

    xs.resize(4);
    xs[0] = r;
    xs[1] = y;
    xs[2] = z;
    xs[3] = t;
    
    if (soft_equiv((*expression)(xs),
                   (((1 && 1.3) || !(y<-1))/5.
                    + (2>1) * r*square(2.7-1.1*z))*t))
    {
        ut.passes("expression successfully evaluated");
    }
    else
    {
        ut.failure("expression NOT successfully evaluated");
    }

    tokens = string("20*(r>=1.1*m && z<=1.5*m || r>=2.0*m)");

    expression = Expression::parse(4, variable_map, tokens);

    if (expression != SP<Expression>())
    {
        ut.passes("expression successfully parsed");
    }
    else
    {
        ut.failure("expression NOT successfully parsed");
    }        

    if (soft_equiv((*expression)(xs),
                   20.0*((r>=1.1) && ((z<=1.5) || (r>=2.0)))))
    {
        ut.passes("expression successfully evaluated");
    }
    else
    {
        ut.failure("expression NOT successfully evaluated");
    }

    tokens = String_Token_Stream("(1 && (4>=6 || 4>6 || 6<4 || 6<=4 || !0))"
                                 "* ( (r/m)^(t/s) + -3 - z/m)");
    
    expression = Expression::parse(4,
                                   variable_map,
                                   tokens);

    if (tokens.error_count()==0 && tokens.lookahead().type()==EXIT)
    {
        ut.passes("expression successfully parsed");
    }
    else
    {
        ut.failure("expression NOT successfully parsed");
        cerr << tokens.messages() << endl;
    }        

    if (soft_equiv((*expression)(xs),
                    pow(r, t) + -3 - z))
    {
        ut.passes("expression successfully evaluated");
    }
    else
    {
        ut.failure("expression NOT successfully evaluated");
    }

    if (!expression->is_constant() &&
        !expression->is_constant(0) &&
        expression->is_constant(1))
    {
        ut.passes("is_constant good");
    }
    else
    {
        ut.failure("is_constant NOT good");
    }

    tokens = String_Token_Stream("exp(-0.5*r/m)*(3*cos(2*y/m) + 5*sin(3*y/m))");

    expression = Expression::parse(4, variable_map, tokens);

    if (expression != SP<Expression>())
    {
        ut.passes("expression successfully parsed");
    }
    else
    {
        ut.failure("expression NOT successfully parsed");
    }        

    if (soft_equiv((*expression)(xs),
                   exp(-0.5*r)*(3*cos(2*y) + 5*sin(3*y))))
    {
        ut.passes("expression successfully evaluated");
    }
    else
    {
        ut.failure("expression NOT successfully evaluated");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstExpression(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing " << argv[0] << ", "
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing " << argv[0] << ", "
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails? EXIT_FAILURE : EXIT_SUCCESS;
}   

//---------------------------------------------------------------------------//
//                        end of tstExpression.cc
//---------------------------------------------------------------------------//
