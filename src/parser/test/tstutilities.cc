//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstutilities.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief
 * \note   Copyright (C) 2003-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../utilities.hh"
#include "../File_Token_Stream.hh"
#include "../String_Token_Stream.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/path.hh"
#include "units/PhysicalConstants.hh"
#include "units/PhysicalConstantsSI.hh"

#include <limits>

using namespace std;
using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstutilities(UnitTest &ut)
{
    std::cout << "Running test tstutilities()..." << std::endl;

    // Build path for the input file "utilities.inp"
    string const inputFile(ut.getTestInputPath()
                           + std::string("utilities.inp") );
    
    // Generate a File_Token_Stream
    File_Token_Stream tokens( inputFile );

    // Try to read a real number.

    double d = parse_real(tokens);
    if (tokens.error_count() != 0 || d != 5.)
    {
        ut.failure("real NOT successfully parsed");
    }
    else
    {
        ut.passes("real successfully parsed");
    }

    // Try to read an integer.

    int i = parse_integer(tokens);
    if (tokens.error_count() != 0 || i != 1)
    {
        ut.failure("integer NOT successfully parsed");
    }
    else
    {
        ut.passes("integer successfully parsed");
    }

    // Try to read a negative integer.

    i = parse_integer(tokens);
    if (tokens.error_count() != 0 || i != -3)
    {
        ut.failure("integer NOT successfully parsed");
    }
    else
    {
        ut.passes("integer successfully parsed");
    }

    // Try to read an unsigned integer.

    i = parse_unsigned_integer(tokens);
    if (tokens.error_count() != 0 || i != 4)
    {
        ut.failure("integer NOT successfully parsed");
    }
    else
    {
        ut.passes("integer successfully parsed");
    }

    // Try to read a positive integer.

    i = parse_positive_integer(tokens);
    if (tokens.error_count() != 0 || i != 1198)
    {
        ut.failure("positive integer NOT successfully parsed");
    }
    else
    {
        ut.passes("positive integer successfully parsed");
    }

    // Try to read an integer as a real.

    d = parse_real(tokens);
    if (tokens.error_count() != 0 || d != 2.)
    {
        ut.failure("integer NOT successfully parsed as real");
    }
    else
    {
        ut.passes("integer successfully parsed as real");
    }

    // Try to read some vectors.

    double v[3];
    parse_vector(tokens, v);
    Token token = tokens.shift();
    if (v[0] == 3. && v[1] == 0.0 && v[2] == 0.0 && 
        token.type() == KEYWORD && token.text() == "stop")
    {
        ut.passes("1-D vector successfully parsed");
    }
    else
    {
        ut.failure("1-D vector NOT successfully parsed");
    }

    parse_vector(tokens, v);
    token = tokens.shift();
    if (v[0] == 1. && v[1] == 2.0 && v[2] == 0.0 && 
        token.type() == KEYWORD && token.text() == "stop")
    {
        ut.passes("2-D vector successfully parsed");
    }
    else
    {
        ut.failure("2-D vector NOT successfully parsed");
    }

    parse_vector(tokens, v);
    if (v[0] == 4. && v[1] == 3.0 && v[2] == 2.0
        && tokens.shift().text()=="stop")
    {
        ut.passes("3-D vector successfully parsed");
    }
    else
    {
        ut.failure("3-D vector NOT successfully parsed");
    }
    unsigned w[3];
    parse_unsigned_vector(tokens, w, 3);
    token = tokens.shift();
    if (w[0] == 3 && w[1] == 2 && w[2] == 1 &&
        token.type() == KEYWORD && token.text() == "stop")
    {
        ut.passes("Vector of unsigned successfully parsed");
    }
    else
    {
        ut.failure("Vector of unsigned NOT successfully parsed");
    }

    // Try to read some unit expressions

    Unit one = {0,0,0,0,0,0,0,0,0, 1};

    Unit left = parse_unit(tokens);
    if (left!=J) ut.failure("test fails");

    left = parse_unit(tokens);
    Unit right = parse_unit(tokens);
    if (left!=right || left!=C) ut.failure("test fails");

    left = parse_unit(tokens);
    cout << left << endl;
    if (left!=1/s) ut.failure("test fails");

    left = parse_unit(tokens);
    right = parse_unit(tokens);
    if (left!=right || left!=N) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=J/K) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=J/mol) ut.failure("test fails");

    left = parse_unit(tokens);
    right = parse_unit(tokens);
    if (left!=right || left!=lm) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=rad/s) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=one)
    {
        ut.failure("dyne definition did NOT check out");
    }
    else
    {
        ut.passes("dyne definition checks out");
    }

    left = parse_unit(tokens);
    if (left!=one)
    {
        ut.failure("erg definition did NOT check out");
    }
    else
    {
        ut.passes("erg definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, cm) || !soft_equiv(left.conv, 0.0254))
    {
        ut.failure("inch definition did NOT check out");
    }
    else
    {
        ut.passes("inch definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, one) || !soft_equiv(left.conv, 12.0))
    {
        ut.failure("foot definition did NOT check out");
    }
    else
    {
        ut.passes("foot definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, one) || !soft_equiv(left.conv, 4.448221615))
    {
        ut.failure("pound definition did NOT check out");
    }
    else
    {
        ut.passes("pound definition checks out");
    }

    left = parse_unit(tokens);
    if (!is_compatible(left, one))
    {
        ut.failure("keV definition did NOT check out");
    }
    else
    {
        ut.passes("keV definition checks out");
    }

    left = parse_unit(tokens);
    if (left!=J) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=J) ut.failure("test fails");

    left = parse_unit(tokens);
    if (left!=K)
    {
        ut.failure("K definition did NOT check out");
    }
    else
    {
        ut.passes("K definition checks out");
    }

    left = parse_unit(tokens);
    if (left!=sr)
    {
        ut.failure("sr definition did NOT check out");
    }
    else
    {
        ut.passes("sr definition checks out");
    }

    // Now see if we catch a bogus unit expression.
    try
    {
        left = parse_unit(tokens);
        ut.failure("did NOT successfully catch bogus unit");
    }
    catch (const Syntax_Error &)
    {
        ut.passes("successfully caught bogus unit");
    }



    // Try to read some dimensioned quantities.
    double length = parse_quantity(tokens, rtt_parser::m, "length");
    if (fabs(length-3.0)<=std::numeric_limits<double>::epsilon())
    {
        ut.passes("length successfully parsed");
    }
    else
    {
        ut.failure("length NOT successfully parsed");
    }

    double energy = parse_quantity(tokens, rtt_parser::J, "energy");
    if (fabs(energy-2.3e-7)<=std::numeric_limits<double>::epsilon())
    {
        ut.passes("cgs energy successfully parsed");
    }
    else
    {
        ut.failure("cgs energy NOT successfully parsed");
    }

    unsigned old_error_count = tokens.error_count();
    length = parse_quantity(tokens, rtt_parser::m, "length");
    if (tokens.error_count()==old_error_count)
    {
        ut.failure("bad length NOT successfully detected");
    }
    else
    {
        ut.passes("bad length successfully detected");
    }

    old_error_count = tokens.error_count();
    double T = parse_temperature(tokens);
    if (tokens.error_count()!=old_error_count ||
        !soft_equiv(T, 273.16))
    {
        ut.failure("temperature NOT successfully parsed");
    }
    else
    {
        ut.passes("temperature successfully parsed");
    }
    T = parse_temperature(tokens);
    if (tokens.error_count()!=old_error_count ||
        !soft_equiv(T, rtt_units::EV2K))
    {
        ut.failure("temperature NOT successfully parsed");
    }
    else
    {
        ut.passes("temperature successfully parsed");
    }
    {
        String_Token_Stream tokens("3.0 m");
        parse_temperature(tokens);
        if (tokens.error_count()==0)
        {
            ut.failure("did NOT detect bad temperature units");
        }
        else
        {
            ut.passes("correctly detected bad temperature units");
        }
    }
    {
        String_Token_Stream tokens("-3.0 K");
        double const T = parse_temperature(tokens);
        if (tokens.error_count()==0 || T!=0)
        {
            ut.failure("did NOT detect negative temperature");
        }
        else
        {
            ut.passes("correctly detected negative temperature");
        }
    }

    // Try reading sequence of quantities with signs
    // Remember(double T=)
    parse_quantity(tokens, J, "energy");
    if (tokens.error_count()!=old_error_count)
    {
        ut.failure("second negative quantity NOT successfully parsed");
    }
    else
    {
        ut.passes("second negative quantity successfully parsed");
    }

    // Try reading a manifest string.
    string parsed_string = parse_manifest_string(tokens);
    if (parsed_string!="manifest string")
    {
        ut.failure("manifest string NOT successfully parsed");
    }
    else
    {
        ut.passes("manifest string successfully parsed");
    }

    {
        String_Token_Stream tokens("bad");
        try
        {
            parse_manifest_string(tokens);
            ut.failure("did NOT detect bad manifest string");
        }
        catch(...)
        {
            ut.passes("correctly detected bad manifest string");
        }
    }

    // Try reading a geometry.
    rtt_mesh_element::Geometry geometry = rtt_mesh_element::END_GEOMETRY;
    parse_geometry(tokens, geometry);
    if (geometry != rtt_mesh_element::AXISYMMETRIC)
    {
        ut.failure("axisymmetric geometry NOT successfully parsed");
    }
    else
    {
        ut.passes("geometry successfully parsed");
    }
    geometry = rtt_mesh_element::END_GEOMETRY;
    parse_geometry(tokens, geometry);
    if (geometry != rtt_mesh_element::CARTESIAN)
    {
        ut.failure("cartesian geometry NOT successfully parsed");
    }
    else
    {
        ut.passes("geometry successfully parsed");
    }
    geometry = rtt_mesh_element::END_GEOMETRY;
    parse_geometry(tokens, geometry);
    if (geometry != rtt_mesh_element::SPHERICAL)
    {
        ut.failure("spherical geometry NOT successfully parsed");
    }
    else
    {
        ut.passes("geometry successfully parsed");
    }

    {
        String_Token_Stream string("4.5");
        if (soft_equiv(parse_positive_real(string), 4.5))
        {
            ut.passes("read positive real correctly");
        }
        else
        {
            ut.failure("did NOT read positive real correctly");
        }
    }
    {
        String_Token_Stream tokens("cylindrical, cartesian, xy, bad");
        rtt_mesh_element::Geometry parsed_geometry =
            rtt_mesh_element::AXISYMMETRIC;
        parse_geometry(tokens, parsed_geometry);
        if (tokens.error_count()==0)
        {
            ut.failure("did NOT detect duplicate geometry definition");
        }
        if (parsed_geometry != rtt_mesh_element::AXISYMMETRIC)
        {
            ut.failure("did NOT parse cylindrical geometry correctly");
        }
        parse_geometry(tokens, parsed_geometry);
        if (parsed_geometry != rtt_mesh_element::CARTESIAN)
        {
            ut.failure("did NOT parse cartesian geometry correctly");
        }
        parse_geometry(tokens, parsed_geometry);
        if (parsed_geometry != rtt_mesh_element::CARTESIAN)
        {
            ut.failure("did NOT parse xy geometry correctly");
        }
        try
        {
            parse_geometry(tokens, parsed_geometry);
            ut.failure("did NOT catch bad geometry");
        }
        catch(...)
        {
        }
    }

    {
        String_Token_Stream string("4.5");
        try
        {
            parse_unsigned_integer(string);
            ut.failure("did NOT detect bad uinteger correctly");
        }
        catch (...)
        {        
            ut.passes("detected bad uinteger correctly");
        }
    }
    {
        String_Token_Stream string("333333333333333333333333");
        parse_unsigned_integer(string);
        if (string.error_count()==0)
        {
            ut.failure("did NOT detect bad uinteger correctly");
        }
        else
        {        
            ut.passes("detected bad uinteger correctly");
        }
    }
    {
        String_Token_Stream string("0");
        parse_positive_integer(string);
        if (string.error_count()==0)
        {
            ut.failure("did NOT detect bad positive integer correctly");
        }
        else
        {        
            ut.passes("detected bad positive integer correctly");
        }
    }
    {
        String_Token_Stream string("+3");
        if (parse_integer(string) == 3)
        {        
            ut.passes("correctly parsed +3");
        }
        else
        {
            ut.failure("did NOT correctly parse +3");
        }
    }
    {
        String_Token_Stream string("4.5");
        try
        {
            parse_integer(string);
            ut.failure("did NOT detect bad integer correctly");
        }
        catch (...)
        {        
            ut.passes("detected bad integer correctly");
        }
    }
    {
        String_Token_Stream string("333333333333333333333333");
        parse_integer(string);
        if (string.error_count()==0)
        {
            ut.failure("did NOT detect bad integer correctly");
        }
        else
        {        
            ut.passes("detected bad integer correctly");
        }
    }
    {
        String_Token_Stream string("-6.5");
        if (at_real(string))
        {
            ut.passes("detected at real correctly");
        }
        else
        {        
            ut.failure("did NOT detect at real for -6.5 correctly");
        }
    }
    {
        String_Token_Stream string("+3f");
        if (at_real(string))
        {
            ut.passes("detected at real correctly");
        }
        else
        {        
            ut.failure("did NOT detect at real for +3 correctly");
        }
    }
    {
        String_Token_Stream string("+3f");
        if (parse_real(string)==3.0)
        {
            ut.passes("parsed real correctly");
        }
        else
        {        
            ut.failure("did NOT parse real +3 correctly");
        }
    }
    {
        String_Token_Stream string("abba");
        try
        {
            parse_real(string);
            ut.failure("did NOT detect bad real correctly");
        }
        catch (...)
        {        
            ut.passes("detected bad real correctly");
        }
    }
    {
        String_Token_Stream string("1.8e10000");
        parse_real(string);
        if (string.error_count()==0)
        {
            ut.failure("did NOT detect real overflow correctly");
        }
        else
        {        
            ut.passes("detected real overflow correctly");
        }
    }
    {
        String_Token_Stream string("-8");
        parse_positive_real(string);
        if (string.error_count()==0)
        {
            ut.failure("did NOT detect nonpositive real correctly");
        }
        else
        {        
            ut.passes("detected nonpositive real correctly");
        }
    }
    {
        String_Token_Stream string("1 2 3");
        unsigned x[4];
        parse_unsigned_vector(string, x, 4);
        if (string.error_count()==0)
        {
            ut.failure("did NOT detect too short vector correctly");
        }
        else
        {        
            ut.passes("detected too short vector correctly");
        }
    }
    {
        String_Token_Stream string("fred");
        Unit const unit = parse_unit(string);
        if (!is_compatible(unit, dimensionless))
        {
            ut.failure("did NOT detect no unit");
        }
        else
        {        
            ut.passes("detected no unit correctly");
        }

    }
    {
        String_Token_Stream string("true");
        bool const b = parse_bool(string);
        if (!b)
        {
            ut.failure("did NOT parse bool");
        }
        else
        {        
            ut.passes("parsed bool correctly");
        }

    }
    {
        String_Token_Stream string("278 K");
        double const T = parse_temperature(string);
        if (soft_equiv(T, 278.))
        {
            ut.passes("parsed temperature correctly");
        }
        else
        {        
            ut.failure("did NOT parse temperature correctly");
        }

    }
    {
        String_Token_Stream string("0.0");
        double const T = parse_nonnegative_real(string);
        if (T==0.0)
        {
            ut.passes("parsed nonnegative real correctly");
        }
        else
        {        
            ut.failure("did NOT parse nonnegative real correctly");
        }

    }
    {
        String_Token_Stream string("5 cm");
        double const T = parse_quantity(string, m, "length");
        if (soft_equiv(T, 0.05))
        {
            ut.passes("parsed quantity correctly");
        }
        else
        {        
            ut.failure("did NOT parse quantity correctly");
        }
    }
    {
        String_Token_Stream string("278*K");
        map<std::string, pair<unsigned, Unit> > variable_map;
        Unit one = {0,0,0,0,0,0,0,0,0, 1};
        variable_map[std::string("x")] = std::make_pair(0,one);
        SP<Expression> T = parse_temperature(string, 1, variable_map);
        vector<double> x(1,0.0);
        if (soft_equiv((*T)(x), 278.))
        {
            ut.passes("parsed temperature correctly");
        }
        else
        {        
            ut.failure("did NOT parse temperature correctly");
        }

    }
    {
        String_Token_Stream string("J");
        map<std::string, pair<unsigned, Unit> > variable_map;
        Unit one = {0,0,0,0,0,0,0,0,0, 1};
        variable_map[std::string("x")] = std::make_pair(0,one);
        SP<Expression> T = parse_temperature(string, 1, variable_map);
        vector<double> x(1,0.0);
        if (soft_equiv((*T)(x), 1.0/rtt_units::boltzmannSI))
        {
            ut.passes("parsed temperature correctly");
        }
        else
        {        
            ut.failure("did NOT parse temperature correctly");
        }

    }
    {
        String_Token_Stream string("278*K");
        map<std::string, pair<unsigned, Unit> > variable_map;
        Unit one = {0,0,0,0,0,0,0,0,0, 1};
        variable_map[std::string("x")] = std::make_pair(0,one);
        SP<Expression> T = parse_quantity(string, K, "temperature", 1, variable_map);
        vector<double> x(1,0.0);
        if (soft_equiv((*T)(x), 278.))
            ut.passes("parsed temperature correctly");
        else
            ut.failure("did NOT parse temperature correctly");
    }

    // Screw around with the internal unit system and with optional unit
    // expressions.
    {
        using namespace rtt_units;
        
        String_Token_Stream quantity_with_units("3e10 cm/s");
        String_Token_Stream bare_quantity("3e10");

        // Check defaults
        double c = parse_quantity(quantity_with_units, m/s, "velocity");
        if (quantity_with_units.error_count()==0 && soft_equiv(c, 3e8))
        {
            ut.passes("parsed quantity with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to SI correctly");
        }
        quantity_with_units.rewind();

        c = parse_quantity(bare_quantity, m/s, "velocity");
        if (bare_quantity.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare quantity");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare quantity");
        }
        bare_quantity.rewind();

        // Turn off mandatory units
        set_unit_expressions_are_required(false);
       
        c = parse_quantity(quantity_with_units, m/s, "velocity");
        if (quantity_with_units.error_count()==0 && soft_equiv(c, 3e8))
        {
            ut.passes("parsed quantity with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to SI correctly");
        }
        quantity_with_units.rewind();

        c = parse_quantity(bare_quantity, m/s, "velocity");
        if (bare_quantity.error_count()==0 && soft_equiv(c, 3e10))
        {
            ut.passes("parsed bare quantity to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse bare quantity to SI correctly");
        }
        bare_quantity.rewind();

        // Turn mandatory units back on but switch internal units to cgs
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().CGS()));
       
        c = parse_quantity(quantity_with_units, m/s, "velocity");
        if (quantity_with_units.error_count()==0 && soft_equiv(c, 3e10))
        {
            ut.passes("parsed quantity with units to cgs correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to cgs correctly");
        }
        quantity_with_units.rewind();

        c = parse_quantity(bare_quantity, m/s, "velocity");
        if (bare_quantity.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare quantity");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare quantity");
        }
        bare_quantity.rewind();

        // Turn mandatory units off again
        set_unit_expressions_are_required(false);
       
        c = parse_quantity(quantity_with_units, m/s, "velocity");
        if (quantity_with_units.error_count()==0 && soft_equiv(c, 3e10))
        {
            ut.passes("parsed quantity with units to cgs correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to cgs correctly");
        }
        quantity_with_units.rewind();

        c = parse_quantity(bare_quantity, m/s, "velocity");
        if (bare_quantity.error_count()==0 && soft_equiv(c, 3e10))
        {
            ut.passes("parsed bare quantity to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse bare quantity to SI correctly");
        }
        bare_quantity.rewind();
    }
    // Screw around with unit expression settings, as before, but for
    // temperature expressions in K.
    {
        using namespace rtt_units;
        
        String_Token_Stream quantity_with_units("273 K");
        String_Token_Stream bare_quantity("273");

        // Check defaults
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().SI()));

        double T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, 273.))
        {
            ut.passes("parsed T with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse T with units to SI correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare quantity");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare quantity");
        }
        bare_quantity.rewind();

        // Turn off mandatory units
        set_unit_expressions_are_required(false);
       
        T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, 273.))
        {
            ut.passes("parsed quantity with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to SI correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()==0 && soft_equiv(T, 273.))
        {
            ut.passes("parsed bare quantity to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse bare quantity to SI correctly");
        }
        bare_quantity.rewind();

        // Turn mandatory units back on but switch internal units to X4
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().X4()));
       
        T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, 273.e-3/EV2K))
        {
            ut.passes("parsed quantity with units to X4 correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to X4 correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare quantity");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare quantity");
        }
        bare_quantity.rewind();

        // Turn mandatory units off again
        set_unit_expressions_are_required(false);
       
        T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, 273.e-3/EV2K))
        {
            ut.passes("parsed quantity with units to X4 correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to X4 correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()==0 && soft_equiv(T, 273.))
        {
            ut.passes("parsed bare quantity to X4 correctly");
        }
        else
        {
            ut.failure("did NOT parse bare quantity to X4 correctly");
        }
        bare_quantity.rewind();
    }
    // Screw around with unit expression settings, as before, but for
    // temperature expressions in keV.
    {
        using namespace rtt_units;
        
        String_Token_Stream quantity_with_units("0.001 keV");
        String_Token_Stream bare_quantity("0.001");

        // Check defaults
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().SI()));

        double T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, EV2K))
        {
            ut.passes("parsed T with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse T with units to SI correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare quantity");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare quantity");
        }
        bare_quantity.rewind();

        // Turn off mandatory units
        set_unit_expressions_are_required(false);
       
        T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, EV2K))
        {
            ut.passes("parsed quantity with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to SI correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()==0 && soft_equiv(T, 0.001))
        {
            ut.passes("parsed bare quantity to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse bare quantity to SI correctly");
        }
        bare_quantity.rewind();

        // Turn mandatory units back on but switch internal units to X4
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().X4()));
       
        T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, 0.001))
        {
            ut.passes("parsed quantity with units to X4 correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to X4 correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare quantity");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare quantity");
        }
        bare_quantity.rewind();

        // Turn mandatory units off again
        set_unit_expressions_are_required(false);
       
        T = parse_temperature(quantity_with_units);
        if (quantity_with_units.error_count()==0 && soft_equiv(T, 0.001))
        {
            ut.passes("parsed quantity with units to X4 correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to X4 correctly");
        }
        quantity_with_units.rewind();

        T = parse_temperature(bare_quantity);
        if (bare_quantity.error_count()==0 && soft_equiv(T, 0.001))
        {
            ut.passes("parsed bare quantity to X4 correctly");
        }
        else
        {
            ut.failure("did NOT parse bare quantity to X4 correctly");
        }
        bare_quantity.rewind();
    }
    // Screw around, etc., but this time with quantity expressions
    {
        using namespace rtt_units;

        std::map<string, pair<unsigned, Unit> > vmap;
        vmap["t"] = pair<unsigned, Unit>(0, s);
        vmap["x"] = pair<unsigned, Unit>(1, m);

        vector<double> var(2, 1U);
        
        String_Token_Stream expression_with_units("0.5*(t+2*x*s/cm)*erg");
        String_Token_Stream expression_appending_units("3.7 erg-s");
        String_Token_Stream bare_expression("0.5*(t+2*x)");

        // Check defaults
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().SI()));

        SP<Expression> c =
            parse_quantity(expression_with_units, erg*s, "action", 2U, vmap);
        
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 0.5*(1.+2*1./0.01)*1e-7))
        {
            ut.passes("parsed expression with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT expression with units to SI correctly");
        }
        expression_with_units.rewind();

        c = parse_quantity(expression_appending_units, erg*s, "action", 2U, vmap);
        
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 3.7*1e-7))
        {
            ut.passes("parsed expression appending units to SI correctly");
        }
        else
        {
            ut.failure("did NOT expression appending units to SI correctly");
        }
        expression_appending_units.rewind();

        c = parse_quantity(bare_expression, erg*s, "action", 2U, vmap);
        if (bare_expression.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare expression");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare expression");
        }
        bare_expression.rewind();

        // Turn off mandatory units
        set_unit_expressions_are_required(false);
       
        c = parse_quantity(expression_with_units, erg*s, "action", 2U, vmap);
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 0.5*(1.+2*1./0.01)*1e-7))
        {
            ut.passes("parsed quantity with units to SI correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to SI correctly");
        }
        expression_with_units.rewind();

        c = parse_quantity(expression_appending_units, erg*s, "action", 2U, vmap);
        
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 3.7*1e-7))
        {
            ut.passes("parsed expression appending units to SI correctly");
        }
        else
        {
            ut.failure("did NOT expression appending units to SI correctly");
        }
        expression_appending_units.rewind();

        c = parse_quantity(bare_expression, erg*s, "action", 2U, vmap);
        if (bare_expression.error_count()==0 &&
            soft_equiv((*c)(var), 0.5*(1.+2*1.)))
        {
            ut.passes("parsed bare expression to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse bare expression to SI correctly");
        }
        bare_expression.rewind();

        // Turn mandatory units back on but switch internal units to cgs
        set_unit_expressions_are_required(true);
        set_internal_unit_system(UnitSystem(UnitSystemType().CGS()));
       
        c = parse_quantity(expression_with_units, erg*s, "action", 2U, vmap);
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 0.5*(1.+2*1./1.)*1.))
        {
            ut.passes("parsed quantity with units to cgs correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to cgs correctly");
        }
        expression_with_units.rewind();

        c = parse_quantity(expression_appending_units, erg*s, "action", 2U, vmap);
        
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 3.7))
        {
            ut.passes("parsed expression appending units to SI correctly");
        }
        else
        {
            ut.failure("did NOT expression appending units to SI correctly");
        }
        expression_appending_units.rewind();

        c = parse_quantity(bare_expression, erg*s, "action", 2U, vmap);
        if (bare_expression.error_count()>0)
        {
            ut.passes("correctly flagged missing units in bare expression");
        }
        else
        {
            ut.failure("did NOT correctly flag missing units in bare expression");
        }
        bare_expression.rewind();

        // Turn mandatory units off again
        set_unit_expressions_are_required(false);
       
        c = parse_quantity(expression_with_units, erg*s, "action", 2U, vmap);
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 0.5*(1.+2*1./1.)*1.))
        {
            ut.passes("parsed quantity with units to cgs correctly");
        }
        else
        {
            ut.failure("did NOT quantity with units to cgs correctly");
        }
        expression_with_units.rewind();

        c = parse_quantity(expression_appending_units, erg*s, "action", 2U, vmap);
        
        if (expression_with_units.error_count()==0 &&
            soft_equiv((*c)(var), 3.7))
        {
            ut.passes("parsed expression appending units to SI correctly");
        }
        else
        {
            ut.failure("did NOT expression appending units to SI correctly");
        }
        expression_appending_units.rewind();

        c = parse_quantity(bare_expression, erg*s, "action", 2U, vmap);
        if (bare_expression.error_count()==0 &&
            soft_equiv((*c)(var), 0.5*(1.+2*1.)))
        {
            ut.passes("parsed bare expression to SI correctly");
        }
        else
        {
            ut.failure("did NOT parse bare expression to SI correctly");
        }
        bare_expression.rewind();
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try { tstutilities(ut); }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstutilities.cc
//---------------------------------------------------------------------------//
