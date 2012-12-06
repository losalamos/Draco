//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/quadrature_test.cc
 * \author Kent G. Budge
 * \brief  Define class quadrature_test
 * \note   Copyright (C) 2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <numeric>

#include "quadrature_test.hh"

#include "parser/String_Token_Stream.hh"

namespace rtt_quadrature
{
using namespace std;
using namespace rtt_parser;

//---------------------------------------------------------------------------------------//
void test_either(UnitTest &ut,
                 SP<Ordinate_Space> const &ordinate_space,
                 Quadrature &quadrature, 
                 unsigned const expansion_order)
{
    vector<Ordinate> const &ordinates = ordinate_space->ordinates();
    unsigned const number_of_ordinates = ordinates.size();

    rtt_mesh_element::Geometry const geometry = ordinate_space->geometry();
    unsigned const dimension = ordinate_space->dimension();
        
    if (ordinate_space->moments()[0] == Moment(0, 0))
    {
        ut.passes("first moment is correct");
    }
    else
    {
        ut.failure("first moment is NOT correct");
    }

    if (number_of_ordinates == ordinate_space->alpha().size())
    {
        ut.passes("alpha size is correct");
    }
    else
    {
        ut.failure("alpha size is NOT correct");
    }

    if (ordinate_space->ordinates().size() == ordinate_space->tau().size())
    {
        ut.passes("tau size is correct");
    }
    else
    {
        ut.failure("tau size is NOT correct");
    }

    if (ordinate_space->expansion_order()==expansion_order)
    {
        ut.passes("expansion order is correct");
    }
    else
    {
        ut.failure("expansion_order is NOT correct");
    }

    vector<unsigned> const &first_angles = ordinate_space->first_angles();
    unsigned const number_of_levels = quadrature.number_of_levels();

    if (geometry==rtt_mesh_element::SPHERICAL)
    {
        if (first_angles.size()==1)
        {
            ut.passes("first angles is correct");
        }
        else
        {
            ut.failure("first angles is NOT correct");
        }

        if (ordinate_space->bookkeeping_coefficient(number_of_ordinates-1)<=0.0)
        {
            ut.failure("bookkeeping coefficient is NOT plausible");
        }
        
        ordinate_space->psi_coefficient(number_of_ordinates-1);
        ordinate_space->source_coefficient(number_of_ordinates-1);
        // check that throws no exception
    }
    else if (geometry==rtt_mesh_element::AXISYMMETRIC)
    {
        if ((dimension>1 && first_angles.size()==number_of_levels) ||
            (dimension==1 && 2*first_angles.size()==number_of_levels))
        {
            ut.passes("first angles is correct");
        }
        else
        {
            ut.failure("first angles is NOT correct");
        }

        if (ordinate_space->bookkeeping_coefficient(number_of_ordinates-1)<=0.0)
        {
            ut.failure("bookkeeping coefficient is NOT plausible");
        }
        
        ordinate_space->psi_coefficient(number_of_ordinates-1);
        ordinate_space->source_coefficient(number_of_ordinates-1);
        // check that throws no exception

        vector<unsigned> const &levels = ordinate_space->levels();
        if (levels.size()==number_of_ordinates)
        {
            ut.passes("levels size is correct");
        }
        else
        {
            ut.failure("levels size is NOT correct");
        }
        for (unsigned i=0; i<number_of_ordinates; ++i)
        {
            if (levels[i]>=number_of_levels)
            {
                ut.failure("levels is NOT in bounds");
                return;
            }
        }

        vector<unsigned> const &moments_per_order =
            ordinate_space->moments_per_order();
            
        if (moments_per_order.size()==expansion_order+1)
        {
            ut.passes("moments_per_order size is correct");
        }
        else
        {
            ut.failure("moments_per_order size is NOT correct");
        }
        for (unsigned i=0; i<=expansion_order; ++i)
        {
            if ((dimension ==1 && moments_per_order[i]!=i/2+1) ||
                (dimension>1 && moments_per_order[i]!=i+1))
            {
                ut.failure("moments_per_order is NOT correct");
                return;
            }
        }
            
        if ((dimension == 1 && number_of_levels == 2*ordinate_space->number_of_levels()) ||
            (dimension>1 && number_of_levels == ordinate_space->number_of_levels()))
        {
            ut.passes("number of levels is consistent");
        }
        else
        {
            ut.failure("number of levels is NOT consistent");
        }
    }
    else
    {
        if (ordinate_space->first_angles().size()==0)
        {
            ut.passes("first angles is correct");
        }
        else
        {
            ut.failure("first angles is NOT correct");
        }
    }

    vector<unsigned> const &reflect_mu = ordinate_space->reflect_mu();
    if (reflect_mu.size() == number_of_ordinates)
    {
        ut.passes("reflect_mu is correct size");
    }
    else
    {
        ut.failure("reflect_mu is NOT correct size");
    }
    for (unsigned i=0; i<number_of_ordinates; ++i)
    {
        if (reflect_mu[i]>=number_of_ordinates)
        {
            ut.failure("reflect_mu is out of bounds");
            return;
        }
        if (ordinates[i].wt() != 0.0 && reflect_mu[reflect_mu[i]] != i)
        {
            ut.failure("reflect_mu is inconsistent");
            return;
        }
    }

    if (dimension>1)
    {
        vector<unsigned> const &reflect_eta = ordinate_space->reflect_eta();
        if (reflect_eta.size() == number_of_ordinates)
        {
            ut.passes("reflect_eta is correct size");
        }
        else
        {
            ut.failure("reflect_eta is NOT correct size");
        }
        for (unsigned i=0; i<number_of_ordinates; ++i)
        {
            if (reflect_eta[i]>=number_of_ordinates)
            {
                ut.failure("reflect_eta is out of bounds");
                return;
            }
            if (reflect_eta[reflect_eta[i]] != i)
            {
                ut.failure("reflect_eta is inconsistent");
                return;
            }
        }

        if (dimension>2)
        {
            vector<unsigned> const &reflect_xi = ordinate_space->reflect_xi();
            if (reflect_xi.size() == number_of_ordinates)
            {
                ut.passes("reflect_xi is correct size");
            }
            else
            {
                ut.failure("reflect_xi is NOT correct size");
            }
            for (unsigned i=0; i<number_of_ordinates; ++i)
            {
                if (reflect_xi[i]>=number_of_ordinates)
                {
                    ut.failure("reflect_xi is out of bounds");
                    return;
                }
                if (reflect_xi[reflect_xi[i]] != i)
                {
                    ut.failure("reflect_xi is inconsistent");
                    return;
                }
            }
        }
    }

    // See if count matches class

    unsigned L = quadrature.number_of_levels();
    unsigned N = ordinate_space->ordinates().size();
    switch (quadrature.quadrature_class())
    {
        case Quadrature::TRIANGLE_QUADRATURE:
            if (dimension==1)
            {
                if (geometry==rtt_mesh_element::CARTESIAN && L != N)
                {
                    ut.failure("ordinate count is wrong for triangular quadrature");
                }
            }
            else if (dimension==3)
            {
                if (L*(L+2)!= N)
                {
                    ut.failure("ordinate count is wrong for triangular quadrature");
                }
            }
            break;
            
        case Quadrature::SQUARE_QUADRATURE:
            if (dimension==3)
            {
                if (2*L*L != N)
                {
                    ut.failure("ordinate count is wrong for square quadrature");
                }
            }
            break;
            
        default:
            if (dimension==3)
            {
                if (4*L>N)
                {
                    ut.failure("ordinate count is too small for level count");
                }
            }
            break;
    }

    // Test that mean and flux are correct
        
    {
        vector<Ordinate> const &ordinates = ordinate_space->ordinates();
        unsigned const N = ordinates.size();
        double J = 0.0;
        double Fx = 0.0, Fy=0.0, Fz=0.0;
        double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
        double const MAGIC = 2.32; // avoid numerical coincidences
        cout << "Ordinates:" << endl;
        for (unsigned i=0; i<N; ++i)
        {
            cout << "  mu = "
                 << setprecision(10) << ordinates[i].mu()<< "  eta = "
                 << setprecision(10) << ordinates[i].eta()<< "  xi = "
                 << setprecision(10) << ordinates[i].xi() << " weight = "
                 << setprecision(10) << ordinates[i].wt() << endl;
            
            J += MAGIC * ordinates[i].wt();
            Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
            Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
            Fy += MAGIC * ordinates[i].eta()*ordinates[i].wt();
            Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
            Fz += MAGIC * ordinates[i].xi()*ordinates[i].wt();
            Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
        }
        if (soft_equiv(J, MAGIC))
        {
            ut.passes("J okay");
        }
        else
        {
            ut.failure("J NOT okay");
        }
        if (soft_equiv(Fx, 0.0))
        {
            ut.passes("xF okay");
        }
        else
        {
            ut.failure("Fx NOT okay");
        }
        if (soft_equiv(Fx2, MAGIC/3.0))
        {
            ut.passes("Fx2 okay");
        }
        else
        {
            cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
            ut.failure("Fx2 NOT okay");
        }
        if (dimension>1)
        {
            if (soft_equiv(Fy, 0.0))
            {
                ut.passes("Fy okay");
            }
            else
            {
                ut.failure("Fy NOT okay");
            }
            if (soft_equiv(Fy2, MAGIC/3.0))
            {
                ut.passes("Fy2 okay");
            }
            else
            {
                cout << "Fy2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                ut.failure("Fy2 NOT okay");
            }
        }
        if (dimension>2)
        {
            if (soft_equiv(Fz, 0.0))
            {
                ut.passes("Fz okay");
            }
            else
            {
                ut.failure("Fz NOT okay");
            }
            if (soft_equiv(Fz2, MAGIC/3.0))
            {
                ut.passes("Fz2 okay");
            }
            else
            {
                cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                ut.failure("Fz2 NOT okay");
            }
        }

        // Look at the moment to discrete and discrete to moment operator
        
        vector<double> M = ordinate_space->M();
        vector<double> D = ordinate_space->D();
        
        unsigned number_of_moments = ordinate_space->number_of_moments();

        if (M.size() == number_of_moments*number_of_ordinates)
        {
            ut.passes("M has right size");
        }
        else
        {
            ut.failure("M does NOT have right size");
        }
        if (D.size() == number_of_moments*number_of_ordinates)
        {
            ut.passes("D has right size");
        }
        else
        {
            ut.failure("D does NOT have right size");
        }

        if (ordinate_space->quadrature_interpolation_model()==GQ)
        {
            for (unsigned m=0; m<number_of_moments; ++m)
            {
                for (unsigned n=0; n<number_of_moments; ++n)
                {
                    double sum = 0.0;
                    for (unsigned a=0; a<number_of_ordinates; ++a)
                    {
                        sum +=  D[a + number_of_ordinates*m] * M[n + a*number_of_moments];
                    }
                    if (m==n)
                    {
                        if (!soft_equiv(sum, 1.0))
                        {
                            ut.failure("diagonal element of M*D NOT 1");
                            return;
                        }
                    }
                    else
                    {
                        if (!soft_equiv(sum, 0.0))
                        {
                            ut.failure("off-diagonal element of M*D NOT 0");
                            return;
                        }
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------------------//
void test_no_axis(UnitTest &ut,
                  Quadrature &quadrature,
                  unsigned const dimension,
                  rtt_mesh_element::Geometry const geometry,
                  unsigned const expansion_order,
                  string const &ordinate_interpolation_model,
                  bool const add_extra_directions,
                  Ordinate_Set::Ordering const ordering)
{
    // Parse the interpolation model
    
    QIM qim = END_QIM;
    String_Token_Stream stokens(ordinate_interpolation_model);
    parse_quadrature_interpolation_model(stokens, qim);

    // Build an angle operator
        
    SP<Ordinate_Space> ordinate_space =
        quadrature.create_ordinate_space(dimension,
                                         geometry,
                                         expansion_order,
                                         add_extra_directions,
                                         ordering,
                                         qim);

    test_either(ut,
                ordinate_space,
                quadrature,
                expansion_order);
}

//---------------------------------------------------------------------------------------//
void test_axis(UnitTest &ut,
               Quadrature &quadrature,
               unsigned const dimension,
               rtt_mesh_element::Geometry const geometry,
               unsigned const expansion_order,
               string const &ordinate_interpolation_model,
               bool const add_extra_directions,
               Ordinate_Set::Ordering const ordering,
               unsigned const mu_axis,
               unsigned const eta_axis)
{
    // Parse the interpolation model
    
    QIM qim = END_QIM;
    String_Token_Stream stokens(ordinate_interpolation_model);
    parse_quadrature_interpolation_model(stokens, qim);

    // Build an angle operator
        
    SP<Ordinate_Space> ordinate_space =
        quadrature.create_ordinate_space(dimension,
                                         geometry,
                                         expansion_order,
                                         mu_axis,
                                         eta_axis,
                                         add_extra_directions,
                                         ordering,
                                         qim);

    test_either(ut,
                ordinate_space,
                quadrature,
                expansion_order);
}

//---------------------------------------------------------------------------------------//
void quadrature_test(UnitTest &ut,
                     Quadrature &quadrature)
{
    cout << "Testing quadrature " << quadrature.name() << endl;
    cout << "  Parse name: " << quadrature.parse_name() << endl;
    switch (quadrature.quadrature_class())
    {
        case Quadrature::INTERVAL_QUADRATURE:
            cout << "  This is an interval quadrature." << endl;
            break;

        case Quadrature::TRIANGLE_QUADRATURE:
            cout << "  This is a triangle quadrature." << endl;
            {
                unsigned L = quadrature.number_of_levels();
                if (L)
                {
                    cout << "  Number of level sets = "
                         << L << endl;
                }
                else
                {
                    ut.failure("no level sets are defined.");
                }
            }
            break;

        case Quadrature::SQUARE_QUADRATURE:
            cout << "  This is a square quadrature." << endl;
            {
                unsigned L = quadrature.number_of_levels();
                if (L)
                {
                    cout << "  Number of level sets = "
                         << L << endl;
                }
                else
                {
                    ut.failure("no level sets are defined.");
                }
            }
            break;

        case Quadrature::OCTANT_QUADRATURE:
            cout << "  This is an octant quadrature." << endl;
            if (quadrature.number_of_levels())
            {
                cout << "  Number of level sets = "
                     << quadrature.number_of_levels() << endl;
            }
            else
            {
                cout << "  No level sets are defined." << endl;
            }
            break;

        default:
            ut.failure("Bad value for quadrature class");
            return;
    }

    // Test moment comparison.

    if (Moment(1,1)==Moment(0,0))
    {
        ut.failure("moment comparison NOT correct");
    }

    // Test ordinate comparison.

    if (Ordinate(0.4, 0.3, sqrt(1.0-0.4*0.4-0.3*0.3), 0.5)==
        Ordinate(0.3, 0.4, sqrt(1.0-0.4*0.4-0.3*0.3), 0.5))
    {
        ut.failure("moment comparison NOT correct");
    }
    if (!(Ordinate(0.4, 0.3, sqrt(1.0-0.4*0.4-0.3*0.3), 0.5)==
          Ordinate(0.4, 0.3, sqrt(1.0-0.4*0.4-0.3*0.3), 0.5)))
    {
        ut.failure("moment comparison NOT correct");
    }

    // Test ordinate access.

    if (Ordinate(1.0, 0.0, 0.0, 0.0).cosines()[0] != 1.0)
    {
        ut.failure("Ordinate::cosines NOT right");
    }

    // Test textifying and parsing.

    string text = quadrature.as_text("\n");
    String_Token_Stream tokens(text);
    SP<Quadrature> parsed_quadrature = Quadrature::parse(tokens);

    if (tokens.error_count())
    {
        ut.failure("Textification and parse did NOT succeed");
    }

    string text2 = parsed_quadrature->as_text("\n");
    if (text2 != text)
    {
        ut.failure("Textification and parse did NOT give identical results");
    }

    // Test various options

    if (!quadrature.has_axis_assignments())
    {
        unsigned dimension = 1;
        rtt_mesh_element::Geometry geometry = rtt_mesh_element::CARTESIAN;
        double norm = 1.0;
        bool add_starting_directions = false;
        bool add_extra_directions = false;
        
        // Build an ordinate set
        
        SP<Ordinate_Set> ordinate_set =
            quadrature.create_ordinate_set(dimension,
                                           geometry,
                                           norm,
                                           add_starting_directions,
                                           add_extra_directions,
                                           Ordinate_Set::LEVEL_ORDERED);
        
        if (ordinate_set->ordinates().size()>=2)
        {
            ut.passes("Ordinate count is plausible");
        }
        else
        {
            ut.failure("Ordinate count is NOT plausible");
        }
        
        if (soft_equiv(ordinate_set->norm(), norm))
        {
            ut.passes("Ordinate norm is correct");
        }
        else
        {
            ut.failure("Ordinate norm is NOT correct");
        }

        ordinate_set->display();
        
        test_no_axis(ut,
                     quadrature,
                     1U, // dimension,
                     rtt_mesh_element::CARTESIAN,
                     1U, // expansion_order,
                     "SN",
                     false, // add_extra_directions,
                     Ordinate_Set::LEVEL_ORDERED);

        if (0)
        {
            // broken
            test_no_axis(ut,
                         quadrature,
                         1U, // dimension,
                         rtt_mesh_element::CARTESIAN,
                         1U, // expansion_order,
                         "GALERKIN",
                         false, // add_extra_directions,
                         Ordinate_Set::LEVEL_ORDERED);
        }
        
        if (quadrature.is_open_interval())
        {
            test_no_axis(ut,
                         quadrature,
                         1U, // dimension,
                         rtt_mesh_element::SPHERICAL,
                         1U, // expansion_order,
                         "SN",
                         false, // add_extra_directions,
                         Ordinate_Set::LEVEL_ORDERED);
        }
    }
    
    if (quadrature.quadrature_class() !=  Quadrature::INTERVAL_QUADRATURE)
    {
        if (0)
        {
            // broken
            test_no_axis(ut,
                         quadrature,
                         2U, // dimension,
                         rtt_mesh_element::CARTESIAN,
                         8U, // expansion_order,
                         "GALERKIN",
                         false, // add_extra_directions,
                         Ordinate_Set::OCTANT_ORDERED);

            test_no_axis(ut,
                         quadrature,
                         2U, // dimension,
                         rtt_mesh_element::CARTESIAN,
                         8U, // expansion_order,
                         "GALERKIN",
                         false, // add_extra_directions,
                         Ordinate_Set::OCTANT_ORDERED);
        }
        
        test_no_axis(ut,
                     quadrature,
                     3U, // dimension,
                     rtt_mesh_element::CARTESIAN,
                     8U, // expansion_order,
                     "SN",
                     false, // add_extra_directions,
                     Ordinate_Set::OCTANT_ORDERED);

        if (quadrature.quadrature_class() ==  Quadrature::TRIANGLE_QUADRATURE)
        {
            test_no_axis(ut,
                         quadrature,
                         3U, // dimension,
                         rtt_mesh_element::CARTESIAN,
                         quadrature.number_of_levels()-1, // expansion_order,
                         "GALERKIN",
                         false, // add_extra_directions,
                         Ordinate_Set::LEVEL_ORDERED);
        }
        
        if (!quadrature.has_axis_assignments())
        {
            if (0)
            {
                // Broken
                test_no_axis(ut,
                             quadrature,
                             1U, // dimension,
                             rtt_mesh_element::AXISYMMETRIC,
                             8U, // expansion_order,
                             "GALERKIN",
                             false, // add_extra_directions,
                             Ordinate_Set::LEVEL_ORDERED);
            }
            
            test_no_axis(ut,
                         quadrature,
                         1U, // dimension,
                         rtt_mesh_element::AXISYMMETRIC,
                         8U, // expansion_order,
                         "SN",
                         false, // add_extra_directions,
                         Ordinate_Set::LEVEL_ORDERED);
            
            test_no_axis(ut,
                         quadrature,
                         2U, // dimension,
                         rtt_mesh_element::AXISYMMETRIC,
                         8U, // expansion_order,
                         "SN",
                         false, // add_extra_directions,
                         Ordinate_Set::LEVEL_ORDERED);
        }

        // Test overriding axis assignments

        test_axis(ut,
                  quadrature,
                  3U, // dimension,
                  rtt_mesh_element::CARTESIAN,
                  8U, // expansion_order,
                  "SN",
                  false, // add_extra_directions,
                  Ordinate_Set::LEVEL_ORDERED,
                  1, // mu_axis
                  2);
    }
    
#if 0
    // **** Build a set of ordinates in 1-D Cartesian geometry, unit norm.

    unsigned dimension = 1;
    rtt_mesh_element::Geometry geometry = rtt_mesh_element::CARTESIAN;
    double norm = 1.0;
    unsigned mu_axis = 0;
    unsigned eta_axis = 1;
    bool add_starting_directions = false;
    bool add_extra_directions = false;
    
    vector<Ordinate> ordinates =
        quadrature.create_ordinates(dimension,
                                    geometry,
                                    norm,
                                    mu_axis,
                                    eta_axis,
                                    add_starting_directions,
                                    add_extra_directions);

    // Build an ordinate set

    SP<Ordinate_Set> ordinate_set =
        quadrature.create_ordinate_set(dimension,
                                       geometry,
                                       norm,
                                       add_starting_directions,
                                       add_extra_directions,
                                       Ordinate_Set::LEVEL_ORDERED);

    if (ordinate_set->ordinates().size()>=2)
    {
        ut.passes("Ordinate count is plausible");
    }
    else
    {
        ut.failure("Ordinate count is NOT plausible");
    }

    if (soft_equiv(ordinate_set->norm(), norm))
    {
        ut.passes("Ordinate norm is correct");
    }
    else
    {
        ut.failure("Ordinate norm is NOT correct");
    }

    // Build an ordinate space

    String_Token_Stream stokens("SN");
    QIM qim = END_QIM;
    parse_quadrature_interpolation_model(stokens, qim);

    SP<Ordinate_Space> ordinate_space =
        quadrature.create_ordinate_space(dimension,
                                         geometry,
                                         1, // expansion order
                                         add_extra_directions,
                                         Ordinate_Set::LEVEL_ORDERED,
                                         qim);

    if (ordinate_set->ordinates().size() == ordinate_space->ordinates().size())
    {
        ut.passes("Ordinate counts agree");
    }
    else
    {
        ut.failure("Ordinate counts do NOT agree");
    }

    if (ordinate_space->ordinates().size() == ordinate_space->alpha().size())
    {
        ut.passes("alpha size is correct");
    }
    else
    {
        ut.failure("alpha size is NOT correct");
    }

    if (ordinate_space->expansion_order()==1)
    {
        ut.passes("expansion order is correct");
    }
    else
    {
        ut.failure("expansion_order is NOT correct");
    }

    // Test that mean and flux are correct

    {
        vector<Ordinate> const &ordinates = ordinate_space->ordinates();
        unsigned const N = ordinates.size();
        double J = 0.0, F = 0.0, F2 = 0.0;
        double const MAGIC = 2.32; // avoid numerical coincidences
        cout << "Ordinates:" << endl;
        for (unsigned i=0; i<N; ++i)
        {
            cout << "  mu = "
                 << setprecision(10) << ordinates[i].mu() << " weight = "
                 << setprecision(10) << ordinates[i].wt() << endl;
            
            J += MAGIC * ordinates[i].wt();
            F += MAGIC * ordinates[i].mu()*ordinates[i].wt();
            F2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
        }
        if (soft_equiv(J, MAGIC))
        {
            ut.passes("J okay");
        }
        else
        {
            ut.failure("J NOT okay");
        }
        if (soft_equiv(F, 0.0))
        {
            ut.passes("F okay");
        }
        else
        {
            ut.failure("F NOT okay");
        }
        if (soft_equiv(F2, MAGIC/3.0))
        {
            ut.passes("F2 okay");
        }
        else
        {
            cout << "F2 = " << F2 << ", expected " << (MAGIC/3.0) << endl;
            ut.failure("F2 NOT okay");
        }
    }

    // **** Build a set of ordinates in 1-D Cartesian geometry, unit norm, Galerkin

#if 0
    // Presently broken
    qim = END_QIM;
    {
        String_Token_Stream stokens("GALERKIN");
        parse_quadrature_interpolation_model(stokens, qim);
    }
    
    ordinate_space =
        quadrature.create_ordinate_space(dimension,
                                         geometry,
                                         1, // expansion order
                                         add_extra_directions,
                                         Ordinate_Set::LEVEL_ORDERED,
                                         qim);

    if (ordinate_set->ordinates().size() == ordinate_space->ordinates().size())
    {
        ut.passes("Ordinate counts agree");
    }
    else
    {
        ut.failure("Ordinate counts do NOT agree");
    }

    if (ordinate_space->ordinates().size() == ordinate_space->alpha().size())
    {
        ut.passes("alpha size is correct");
    }
    else
    {
        ut.failure("alpha size is NOT correct");
    }

    if (ordinate_space->expansion_order()==1)
    {
        ut.passes("expansion order is correct");
    }
    else
    {
        ut.failure("expansion_order is NOT correct");
    }

    // Test that mean and flux are correct

    {
        vector<Ordinate> const &ordinates = ordinate_space->ordinates();
        unsigned const N = ordinates.size();
        double J = 0.0, F = 0.0, F2 = 0.0;
        double const MAGIC = 2.32; // avoid numerical coincidences
        cout << "Ordinates:" << endl;
        for (unsigned i=0; i<N; ++i)
        {
            cout << "  mu = "
                 << setprecision(10) << ordinates[i].mu() << " weight = "
                 << setprecision(10) << ordinates[i].wt() << endl;
            
            J += MAGIC * ordinates[i].wt();
            F += MAGIC * ordinates[i].mu()*ordinates[i].wt();
            F2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
        }
        if (soft_equiv(J, MAGIC))
        {
            ut.passes("J okay");
        }
        else
        {
            ut.failure("J NOT okay");
        }
        if (soft_equiv(F, 0.0))
        {
            ut.passes("F okay");
        }
        else
        {
            ut.failure("F NOT okay");
        }
        if (soft_equiv(F2, MAGIC/3.0))
        {
            ut.passes("F2 okay");
        }
        else
        {
            cout << "F2 = " << F2 << ", expected " << (MAGIC/3.0) << endl;
            ut.failure("F2 NOT okay");
        }
    }
#endif
    
    // **** Try a 3-D if it's an octant quadrature.
        
    if (quadrature.quadrature_class() !=  Quadrature::INTERVAL_QUADRATURE)
    {
        unsigned dimension = 3;
        rtt_mesh_element::Geometry geometry = rtt_mesh_element::CARTESIAN;

        // Build an angle operator
        
        SP<Ordinate_Space> ordinate_space =
            quadrature.create_ordinate_space(dimension,
                                             geometry,
                                             3, // expansion order
                                             add_extra_directions,
                                             Ordinate_Set::LEVEL_ORDERED,
                                             SN);

        unsigned const number_of_ordinates = ordinate_space->ordinates().size();
        
        if (number_of_ordinates == ordinate_space->alpha().size())
        {
            ut.passes("alpha size is correct");
        }
        else
        {
            ut.failure("alpha size is NOT correct");
        }
        
        if (ordinate_space->expansion_order()==3)
        {
            ut.passes("expansion order is correct");
        }
        else
        {
            ut.failure("expansion_order is NOT correct");
        }

        vector<unsigned> const &reflect_mu = ordinate_space->reflect_mu();
        if (reflect_mu.size() == number_of_ordinates)
        {
            ut.passes("reflect_mu is correct size");
        }
        else
        {
            ut.failure("reflect_mu is NOT correct size");
        }
        for (unsigned i=0; i<number_of_ordinates; ++i)
        {
            if (reflect_mu[i]>=number_of_ordinates)
            {
                ut.failure("reflect_mu is out of bounds");
                return;
            }
            if (reflect_mu[reflect_mu[i]] != i)
            {
                ut.failure("reflect_mu is inconsistent");
                return;
            }
        }

        vector<unsigned> const &reflect_eta = ordinate_space->reflect_eta();
        if (reflect_eta.size() == number_of_ordinates)
        {
            ut.passes("reflect_eta is correct size");
        }
        else
        {
            ut.failure("reflect_eta is NOT correct size");
        }
        for (unsigned i=0; i<number_of_ordinates; ++i)
        {
            if (reflect_eta[i]>=number_of_ordinates)
            {
                ut.failure("reflect_eta is out of bounds");
                return;
            }
            if (reflect_eta[reflect_eta[i]] != i)
            {
                ut.failure("reflect_eta is inconsistent");
                return;
            }
        }

        vector<unsigned> const &reflect_xi = ordinate_space->reflect_xi();
        if (reflect_xi.size() == number_of_ordinates)
        {
            ut.passes("reflect_xi is correct size");
        }
        else
        {
            ut.failure("reflect_xi is NOT correct size");
        }
        for (unsigned i=0; i<number_of_ordinates; ++i)
        {
            if (reflect_xi[i]>=number_of_ordinates)
            {
                ut.failure("reflect_xi is out of bounds");
                return;
            }
            if (reflect_xi[reflect_xi[i]] != i)
            {
                ut.failure("reflect_xi is inconsistent");
                return;
            }
        }

        // See if count matches class

        unsigned L = quadrature.number_of_levels();
        unsigned N = ordinate_space->ordinates().size();
        switch (quadrature.quadrature_class())
        {
            case Quadrature::TRIANGLE_QUADRATURE:
                if (L*(L+2) != N)
                {
                    ut.failure("ordinate count is wrong for triangular quadrature");
                }
                break;
            
            case Quadrature::SQUARE_QUADRATURE:
                if (2*L*L != N)
                {
                    ut.failure("ordinate count is wrong for square quadrature");
                }
                break;
            
            default:
                if (4*L>N)
                {
                    ut.failure("ordinate count is too small for level count");
                }
                break;
        }

        // Test that mean and flux are correct
        
        {
            vector<Ordinate> const &ordinates = ordinate_space->ordinates();
            unsigned const N = ordinates.size();
            double J = 0.0;
            double Fx = 0.0, Fy=0.0, Fz=0.0;
            double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
            double const MAGIC = 2.32; // avoid numerical coincidences
            cout << "Ordinates:" << endl;
            for (unsigned i=0; i<N; ++i)
            {
                cout << "  mu = "
                     << setprecision(10) << ordinates[i].mu()<< "  eta = "
                     << setprecision(10) << ordinates[i].eta()<< "  xi = "
                     << setprecision(10) << ordinates[i].xi() << " weight = "
                     << setprecision(10) << ordinates[i].wt() << endl;
            
                J += MAGIC * ordinates[i].wt();
                Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
                Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
                Fy += MAGIC * ordinates[i].eta()*ordinates[i].wt();
                Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
                Fz += MAGIC * ordinates[i].xi()*ordinates[i].wt();
                Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
            }
            if (soft_equiv(J, MAGIC))
            {
                ut.passes("J okay");
            }
            else
            {
                ut.failure("J NOT okay");
            }
            if (soft_equiv(Fx, 0.0))
            {
                ut.passes("xF okay");
            }
            else
            {
                ut.failure("Fx NOT okay");
            }
            if (soft_equiv(Fx2, MAGIC/3.0))
            {
                ut.passes("Fx2 okay");
            }
            else
            {
                cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
                ut.failure("Fx2 NOT okay");
            }
            if (soft_equiv(Fy, 0.0))
            {
                ut.passes("Fy okay");
            }
            else
            {
                ut.failure("Fy NOT okay");
            }
            if (soft_equiv(Fy2, MAGIC/3.0))
            {
                ut.passes("Fy2 okay");
            }
            else
            {
                cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                ut.failure("Fz2 NOT okay");
            }
            if (soft_equiv(Fz, 0.0))
            {
                ut.passes("Fz okay");
            }
            else
            {
                ut.failure("Fz NOT okay");
            }
            if (soft_equiv(Fz2, MAGIC/3.0))
            {
                ut.passes("Fz2 okay");
            }
            else
            {
                cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                ut.failure("Fz2 NOT okay");
            }
        }

        // Build a Galerkin angle operator

        if (quadrature.quadrature_class() == Quadrature::TRIANGLE_QUADRATURE)
        {
            unsigned expansion_order = quadrature.number_of_levels()-1;

            String_Token_Stream stokens("GALERKIN");
            QIM qim = END_QIM;
            parse_quadrature_interpolation_model(stokens, qim);
            
            SP<Ordinate_Space> ordinate_space =
                quadrature.create_ordinate_space(dimension,
                                                 geometry,
                                                 expansion_order,
                                                 add_extra_directions,
                                                 Ordinate_Set::LEVEL_ORDERED,
                                                 qim);

            if (ordinate_space->ordinates().size() == ordinate_space->alpha().size())
            {
                ut.passes("alpha size is correct");
            }
            else
            {
                ut.failure("alpha size is NOT correct");
            }

            if (ordinate_space->expansion_order()==expansion_order)
            {
                ut.passes("expansion order is correct");
            }
            else
            {
                ut.failure("expansion_order is NOT correct");
            }

            if (ordinate_space->first_angles().size()==0)
            {
                ut.passes("first angles is correct");
            }
            else
            {
                ut.failure("first angles is NOT correct");
            }

            // Test that mean and flux are correct
        
            {
                vector<Ordinate> const &ordinates = ordinate_space->ordinates();
                unsigned const N = ordinates.size();
                double J = 0.0;
                double Fx = 0.0, Fy=0.0, Fz=0.0;
                double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
                double const MAGIC = 2.32; // avoid numerical coincidences
                cout << "Ordinates:" << endl;
                for (unsigned i=0; i<N; ++i)
                {
                    cout << "  mu = "
                         << setprecision(10) << ordinates[i].mu()<< "  eta = "
                         << setprecision(10) << ordinates[i].eta()<< "  xi = "
                         << setprecision(10) << ordinates[i].xi() << " weight = "
                         << setprecision(10) << ordinates[i].wt() << endl;
            
                    J += MAGIC * ordinates[i].wt();
                    Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
                    Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
                    Fy += MAGIC * ordinates[i].eta()*ordinates[i].wt();
                    Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
                    Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
                }
                if (soft_equiv(J, MAGIC))
                {
                    ut.passes("J okay");
                }
                else
                {
                    ut.failure("J NOT okay");
                }
                if (soft_equiv(Fx, 0.0))
                {
                    ut.passes("xF okay");
                }
                else
                {
                    ut.failure("Fx NOT okay");
                }
                if (soft_equiv(Fx2, MAGIC/3.0))
                {
                    ut.passes("Fx2 okay");
                }
                else
                {
                    cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
                    ut.failure("Fx2 NOT okay");
                }
                if (soft_equiv(Fy, 0.0))
                {
                    ut.passes("Fy okay");
                }
                else
                {
                    ut.failure("Fy NOT okay");
                }
                if (soft_equiv(Fy2, MAGIC/3.0))
                {
                    ut.passes("Fy2 okay");
                }
                else
                {
                    cout << "Fy2 = " << Fy2 << ", expected " << (MAGIC/3.0) << endl;
                    ut.failure("Fy2 NOT okay");
                }
                if (soft_equiv(Fz, 0.0))
                {
                    ut.passes("Fz okay");
                }
                else
                {
                    ut.failure("Fz NOT okay");
                }
                if (soft_equiv(Fz2, MAGIC/3.0))
                {
                    ut.passes("Fz2 okay");
                }
                else
                {
                    cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                    ut.failure("Fz2 NOT okay");
                }
            }

            // Look at the moment to discrete and discrete to moment operator

            vector<double> M = ordinate_space->M();
            vector<double> D = ordinate_space->D();

            unsigned number_of_moments = ordinate_space->number_of_moments();

            if (M.size() == number_of_moments*number_of_ordinates)
            {
                ut.passes("M has right size");
            }
            else
            {
                ut.failure("M does NOT have right size");
            }
            if (D.size() == number_of_moments*number_of_ordinates)
            {
                ut.passes("D has right size");
            }
            else
            {
                ut.failure("D does NOT have right size");
            }

            for (unsigned m=0; m<number_of_moments; ++m)
            {
                for (unsigned n=0; n<number_of_moments; ++n)
                {
                    double sum = 0.0;
                    for (unsigned a=0; a<number_of_ordinates; ++a)
                    {
                        sum +=  D[a + number_of_ordinates*m] * M[n + a*number_of_moments];
                    }
                    if (m==n)
                    {
                        if (!soft_equiv(sum, 1.0))
                        {
                            ut.failure("diagonal element of M*D NOT 1");
                            return;
                        }
                    }
                    else
                    {
                        if (!soft_equiv(sum, 0.0))
                        {
                            ut.failure("off-diagonal element of M*D NOT 0");
                            return;
                        }
                    }
                }
            }

            // Now do a 2-D axisymmetric space

            dimension = 2;
            geometry = rtt_mesh_element::AXISYMMETRIC;

            qim = END_QIM;
            {
                String_Token_Stream stokens("GALERKIN");
                parse_quadrature_interpolation_model(stokens, qim);
            }
            
            ordinate_space =
                quadrature.create_ordinate_space(dimension,
                                                 geometry,
                                                 expansion_order,
                                                 add_extra_directions,
                                                 Ordinate_Set::LEVEL_ORDERED,
                                                 qim);

            vector<double> const &alpha = ordinate_space->alpha();
            unsigned const number_of_ordinates = ordinate_space->ordinates().size();

            if (alpha.size()==number_of_ordinates)
            {
                ut.passes("alpha size is correct");
            }
            else
            {
                ut.failure("alpha size is NOT correct");
            }

            if (ordinate_space->ordinates().size() == ordinate_space->tau().size())
            {
                ut.passes("tau size is correct");
            }
            else
            {
                ut.failure("tau size is NOT correct");
            }

            if (ordinate_space->expansion_order()==expansion_order)
            {
                ut.passes("expansion order is correct");
            }
            else
            {
                ut.failure("expansion_order is NOT correct");
            }

            vector<unsigned> const &first_angles = ordinate_space->first_angles();
            unsigned const number_of_levels = quadrature.number_of_levels();
            
            if (first_angles.size()==number_of_levels)
            {
                ut.passes("first angles is correct");
            }
            else
            {
                ut.failure("first angles is NOT correct");
            }
            
            if (number_of_levels == ordinate_space->number_of_levels())
            {
                ut.passes("number of levels is consistent");
            }
            else
            {
                ut.failure("number of levels is NOT consistent");
            }
            
            ordinate_space->is_dependent(number_of_ordinates-1);
            // check for exception

            for (unsigned i=0; i<number_of_levels; ++i)
            {
                unsigned const a = first_angles[i];
                if (a>=alpha.size())
                {
                    ut.failure("first_angles is NOT in bounds");
                    return;
                }
                if (!soft_equiv(alpha[a], 0.0))
                {
                    ut.failure("final level alpha is NOT zero");
                    return;
                }
                
            }

            vector<unsigned> const &levels = ordinate_space->levels();
            if (levels.size()==number_of_ordinates)
            {
                ut.passes("levels size is correct");
            }
            else
            {
                ut.failure("levels size is NOT correct");
            }
            for (unsigned i=0; i<number_of_ordinates; ++i)
            {
                if (levels[i]>=number_of_levels)
                {
                    ut.failure("levels is NOT in bounds");
                    return;
                }
            }

            vector<unsigned> const &moments_per_order =
                ordinate_space->moments_per_order();
            
            if (moments_per_order.size()==expansion_order+1)
            {
                ut.passes("moments_per_order size is correct");
            }
            else
            {
                ut.failure("moments_per_order size is NOT correct");
            }
            for (unsigned i=0; i<=expansion_order; ++i)
            {
                if (moments_per_order[i]!=i+1)
                {
                    ut.failure("moments_per_order is NOT correct");
                    return;
                }
            }

            number_of_moments = ordinate_space->number_of_moments();
            
            if (accumulate(moments_per_order.begin(),
                           moments_per_order.end(),
                           0U)
                == number_of_moments)
            {
                ut.passes("number of moments is correct");
            }
            else
            {
                ut.failure("number of moments is NOT correct");
            }

            if (ordinate_space->bookkeeping_coefficient(number_of_ordinates-1)<=0.0)
            {
                ut.failure("bookkeeping coefficient is NOT plausible");
            }

            ordinate_space->psi_coefficient(number_of_ordinates-1);
            ordinate_space->source_coefficient(number_of_ordinates-1);
            // check that throws no exception

            // Test that mean and flux are correct
        
            {
                vector<Ordinate> const &ordinates = ordinate_space->ordinates();
                unsigned const N = ordinates.size();
                double J = 0.0;
                double Fx = 0.0, Fy=0.0;
                double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
                double const MAGIC = 2.32; // avoid numerical coincidences
                cout << "Ordinates:" << endl;
                for (unsigned i=0; i<N; ++i)
                {
                    cout << "  mu = "
                         << setprecision(10) << ordinates[i].mu()<< "  eta = "
                         << setprecision(10) << ordinates[i].eta()<< "  xi = "
                         << setprecision(10) << ordinates[i].xi() << " weight = "
                         << setprecision(10) << ordinates[i].wt() << endl;
            
                    J += MAGIC * ordinates[i].wt();
                    Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
                    Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
                    Fy += MAGIC * ordinates[i].eta()*ordinates[i].wt();
                    Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
                    Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
                }
                if (soft_equiv(J, MAGIC))
                {
                    ut.passes("J okay");
                }
                else
                {
                    ut.failure("J NOT okay");
                }
                if (soft_equiv(Fx, 0.0))
                {
                    ut.passes("xF okay");
                }
                else
                {
                    ut.failure("Fx NOT okay");
                }
                if (soft_equiv(Fx2, MAGIC/3.0))
                {
                    ut.passes("Fx2 okay");
                }
                else
                {
                    cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
                    ut.failure("Fx2 NOT okay");
                }
                if (soft_equiv(Fy, 0.0))
                {
                    ut.passes("Fy okay");
                }
                else
                {
                    ut.failure("Fy NOT okay");
                }
                if (soft_equiv(Fy2, MAGIC/3.0))
                {
                    ut.passes("Fy2 okay");
                }
                else
                {
                    cout << "Fy2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                    ut.failure("Fy2 NOT okay");
                }
                if (soft_equiv(Fz2, MAGIC/3.0))
                {
                    ut.passes("Fz2 okay");
                }
                else
                {
                    cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                    ut.failure("Fz2 NOT okay");
                }
            }

            // Now do a 1-D axisymmetric space

            {
                dimension = 1;
                geometry = rtt_mesh_element::AXISYMMETRIC;

                qim = END_QIM;
                {
                    String_Token_Stream stokens("SN");
                    parse_quadrature_interpolation_model(stokens, qim);
                }
            
                ordinate_space =
                    quadrature.create_ordinate_space(dimension,
                                                     geometry,
                                                     expansion_order,
                                                     add_extra_directions,
                                                     Ordinate_Set::LEVEL_ORDERED,
                                                     qim);

                vector<double> const &alpha = ordinate_space->alpha();
                unsigned const number_of_ordinates = ordinate_space->ordinates().size();

                if (alpha.size()==number_of_ordinates)
                {
                    ut.passes("alpha size is correct");
                }
                else
                {
                    ut.failure("alpha size is NOT correct");
                }

                if (ordinate_space->ordinates().size() == ordinate_space->tau().size())
                {
                    ut.passes("tau size is correct");
                }
                else
                {
                    ut.failure("tau size is NOT correct");
                }

                if (ordinate_space->expansion_order()==expansion_order)
                {
                    ut.passes("expansion order is correct");
                }
                else
                {
                    ut.failure("expansion_order is NOT correct");
                }

                vector<unsigned> const &first_angles = ordinate_space->first_angles();
                unsigned const number_of_levels = quadrature.number_of_levels();
            
                if (2*first_angles.size()==number_of_levels)
                {
                    ut.passes("first angles is correct");
                }
                else
                {
                    ut.failure("first angles is NOT correct");
                }
            
                if (number_of_levels == 2*ordinate_space->number_of_levels())
                {
                    ut.passes("number of levels is consistent");
                }
                else
                {
                    ut.failure("number of levels is NOT consistent");
                }
            
                ordinate_space->is_dependent(number_of_ordinates-1);
                // check for exception

                for (unsigned i=0; 2*i<number_of_levels; ++i)
                {
                    unsigned const a = first_angles[i];
                    if (a>=alpha.size())
                    {
                        ut.failure("first_angles is NOT in bounds");
                        return;
                    }
                    if (!soft_equiv(alpha[a], 0.0))
                    {
                        ut.failure("final level alpha is NOT zero");
                        return;
                    }
                
                }

                vector<unsigned> const &levels = ordinate_space->levels();
                if (levels.size()==number_of_ordinates)
                {
                    ut.passes("levels size is correct");
                }
                else
                {
                    ut.failure("levels size is NOT correct");
                }
                for (unsigned i=0; i<number_of_ordinates; ++i)
                {
                    if (levels[i]>=number_of_levels)
                    {
                        ut.failure("levels is NOT in bounds");
                        return;
                    }
                }

                vector<unsigned> const &moments_per_order =
                    ordinate_space->moments_per_order();
            
                if (moments_per_order.size()==expansion_order+1)
                {
                    ut.passes("moments_per_order size is correct");
                }
                else
                {
                    ut.failure("moments_per_order size is NOT correct");
                }
                for (unsigned i=0; i<=expansion_order; ++i)
                {
                    if (moments_per_order[i]!=(i/2)+1)
                    {
                        ut.failure("moments_per_order is NOT correct");
                        return;
                    }
                }

                unsigned const number_of_moments = ordinate_space->number_of_moments();
            
                if (accumulate(moments_per_order.begin(),
                               moments_per_order.end(),
                               0U)
                    == number_of_moments)
                {
                    ut.passes("number of moments is correct");
                }
                else
                {
                    ut.failure("number of moments is NOT correct");
                }

                // Test that mean and flux are correct
        
                {
                    vector<Ordinate> const &ordinates = ordinate_space->ordinates();
                    unsigned const N = ordinates.size();
                    double J = 0.0;
                    double Fx = 0.0;
                    double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
                    double const MAGIC = 2.32; // avoid numerical coincidences
                    cout << "Ordinates:" << endl;
                    for (unsigned i=0; i<N; ++i)
                    {
                        cout << "  mu = "
                             << setprecision(10) << ordinates[i].mu()<< "  eta = "
                             << setprecision(10) << ordinates[i].eta()<< "  xi = "
                             << setprecision(10) << ordinates[i].xi() << " weight = "
                             << setprecision(10) << ordinates[i].wt() << endl;
            
                        J += MAGIC * ordinates[i].wt();
                        Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
                        Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
                        Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
                        Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
                    }
                    if (soft_equiv(J, MAGIC))
                    {
                        ut.passes("J okay");
                    }
                    else
                    {
                        ut.failure("J NOT okay");
                    }
                    if (soft_equiv(Fx, 0.0))
                    {
                        ut.passes("xF okay");
                    }
                    else
                    {
                        ut.failure("Fx NOT okay");
                    }
                    if (soft_equiv(Fx2, MAGIC/3.0))
                    {
                        ut.passes("Fx2 okay");
                    }
                    else
                    {
                        cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fx2 NOT okay");
                    }
                    if (soft_equiv(Fy2, MAGIC/3.0))
                    {
                        ut.passes("Fy2 okay");
                    }
                    else
                    {
                        cout << "Fy2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fy2 NOT okay");
                    }
                    if (soft_equiv(Fz2, MAGIC/3.0))
                    {
                        ut.passes("Fz2 okay");
                    }
                    else
                    {
                        cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fz2 NOT okay");
                    }
                }
            }

            // Now do a 1-D axisymmetric space, Galerkin

            if (0)
                // Broken
            {
                dimension = 1;
                geometry = rtt_mesh_element::AXISYMMETRIC;

                qim = END_QIM;
                {
                    String_Token_Stream stokens("GALERKIN");
                    parse_quadrature_interpolation_model(stokens, qim);
                }
            
                ordinate_space =
                    quadrature.create_ordinate_space(dimension,
                                                     geometry,
                                                     expansion_order,
                                                     add_extra_directions,
                                                     Ordinate_Set::LEVEL_ORDERED,
                                                     qim);

                vector<double> const &alpha = ordinate_space->alpha();
                unsigned const number_of_ordinates = ordinate_space->ordinates().size();

                if (alpha.size()==number_of_ordinates)
                {
                    ut.passes("alpha size is correct");
                }
                else
                {
                    ut.failure("alpha size is NOT correct");
                }

                if (ordinate_space->ordinates().size() == ordinate_space->tau().size())
                {
                    ut.passes("tau size is correct");
                }
                else
                {
                    ut.failure("tau size is NOT correct");
                }

                if (ordinate_space->expansion_order()==expansion_order)
                {
                    ut.passes("expansion order is correct");
                }
                else
                {
                    ut.failure("expansion_order is NOT correct");
                }

                vector<unsigned> const &first_angles = ordinate_space->first_angles();
                unsigned const number_of_levels = quadrature.number_of_levels();
            
                if (2*first_angles.size()==number_of_levels)
                {
                    ut.passes("first angles is correct");
                }
                else
                {
                    ut.failure("first angles is NOT correct");
                }
            
                if (number_of_levels == 2*ordinate_space->number_of_levels())
                {
                    ut.passes("number of levels is consistent");
                }
                else
                {
                    ut.failure("number of levels is NOT consistent");
                }
            
                ordinate_space->is_dependent(number_of_ordinates-1);
                // check for exception

                for (unsigned i=0; 2*i<number_of_levels; ++i)
                {
                    unsigned const a = first_angles[i];
                    if (a>=alpha.size())
                    {
                        ut.failure("first_angles is NOT in bounds");
                        return;
                    }
                    if (!soft_equiv(alpha[a], 0.0))
                    {
                        ut.failure("final level alpha is NOT zero");
                        return;
                    }
                
                }

                vector<unsigned> const &levels = ordinate_space->levels();
                if (levels.size()==number_of_ordinates)
                {
                    ut.passes("levels size is correct");
                }
                else
                {
                    ut.failure("levels size is NOT correct");
                }
                for (unsigned i=0; i<number_of_ordinates; ++i)
                {
                    if (levels[i]>=number_of_levels)
                    {
                        ut.failure("levels is NOT in bounds");
                        return;
                    }
                }

                vector<unsigned> const &moments_per_order =
                    ordinate_space->moments_per_order();
            
                if (moments_per_order.size()==expansion_order+1)
                {
                    ut.passes("moments_per_order size is correct");
                }
                else
                {
                    ut.failure("moments_per_order size is NOT correct");
                }
                for (unsigned i=0; i<=expansion_order; ++i)
                {
                    if (moments_per_order[i]!=(i/2)+1)
                    {
                        ut.failure("moments_per_order is NOT correct");
                        return;
                    }
                }

                unsigned const number_of_moments = ordinate_space->number_of_moments();
            
                if (accumulate(moments_per_order.begin(),
                               moments_per_order.end(),
                               0U)
                    == number_of_moments)
                {
                    ut.passes("number of moments is correct");
                }
                else
                {
                    ut.failure("number of moments is NOT correct");
                }

                // Test that mean and flux are correct
        
                {
                    vector<Ordinate> const &ordinates = ordinate_space->ordinates();
                    unsigned const N = ordinates.size();
                    double J = 0.0;
                    double Fx = 0.0;
                    double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
                    double const MAGIC = 2.32; // avoid numerical coincidences
                    cout << "Ordinates:" << endl;
                    for (unsigned i=0; i<N; ++i)
                    {
                        cout << "  mu = "
                             << setprecision(10) << ordinates[i].mu()<< "  eta = "
                             << setprecision(10) << ordinates[i].eta()<< "  xi = "
                             << setprecision(10) << ordinates[i].xi() << " weight = "
                             << setprecision(10) << ordinates[i].wt() << endl;
            
                        J += MAGIC * ordinates[i].wt();
                        Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
                        Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
                        Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
                        Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
                    }
                    if (soft_equiv(J, MAGIC))
                    {
                        ut.passes("J okay");
                    }
                    else
                    {
                        ut.failure("J NOT okay");
                    }
                    if (soft_equiv(Fx, 0.0))
                    {
                        ut.passes("xF okay");
                    }
                    else
                    {
                        ut.failure("Fx NOT okay");
                    }
                    if (soft_equiv(Fx2, MAGIC/3.0))
                    {
                        ut.passes("Fx2 okay");
                    }
                    else
                    {
                        cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fx2 NOT okay");
                    }
                    if (soft_equiv(Fy2, MAGIC/3.0))
                    {
                        ut.passes("Fy2 okay");
                    }
                    else
                    {
                        cout << "Fy2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fy2 NOT okay");
                    }
                    if (soft_equiv(Fz2, MAGIC/3.0))
                    {
                        ut.passes("Fz2 okay");
                    }
                    else
                    {
                        cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fz2 NOT okay");
                    }
                }
            }

            {
                // Now do a 2-D Cartesian space

                dimension = 2;
                geometry = rtt_mesh_element::CARTESIAN;

                qim = END_QIM;
                {
                    String_Token_Stream stokens("SN");
                    parse_quadrature_interpolation_model(stokens, qim);
                }
            
                ordinate_space =
                    quadrature.create_ordinate_space(dimension,
                                                     geometry,
                                                     expansion_order,
                                                     add_extra_directions,
                                                     Ordinate_Set::LEVEL_ORDERED,
                                                     qim);

                vector<double> const &alpha = ordinate_space->alpha();
                unsigned const number_of_ordinates = ordinate_space->ordinates().size();

                if (alpha.size()==number_of_ordinates)
                {
                    ut.passes("alpha size is correct");
                }
                else
                {
                    ut.failure("alpha size is NOT correct");
                }

                if (ordinate_space->ordinates().size() == ordinate_space->tau().size())
                {
                    ut.passes("tau size is correct");
                }
                else
                {
                    ut.failure("tau size is NOT correct");
                }

                if (ordinate_space->expansion_order()==expansion_order)
                {
                    ut.passes("expansion order is correct");
                }
                else
                {
                    ut.failure("expansion_order is NOT correct");
                }

                unsigned const number_of_levels = quadrature.number_of_levels();
            
                ordinate_space->is_dependent(number_of_ordinates-1);
                // check for exception

                vector<unsigned> const &levels = ordinate_space->levels();
                if (levels.size()==number_of_ordinates)
                {
                    ut.passes("levels size is correct");
                }
                else
                {
                    ut.failure("levels size is NOT correct");
                }
                for (unsigned i=0; i<number_of_ordinates; ++i)
                {
                    if (levels[i]>=number_of_levels)
                    {
                        ut.failure("levels is NOT in bounds");
                        return;
                    }
                }

                vector<unsigned> const &moments_per_order =
                    ordinate_space->moments_per_order();
            
                if (moments_per_order.size()==expansion_order+1)
                {
                    ut.passes("moments_per_order size is correct");
                }
                else
                {
                    ut.failure("moments_per_order size is NOT correct");
                }
                for (unsigned i=0; i<=expansion_order; ++i)
                {
                    if (moments_per_order[i]!=i+1)
                    {
                        ut.failure("moments_per_order is NOT correct");
                        return;
                    }
                }

                unsigned const number_of_moments = ordinate_space->number_of_moments();
            
                if (accumulate(moments_per_order.begin(),
                               moments_per_order.end(),
                               0U)
                    == number_of_moments)
                {
                    ut.passes("number of moments is correct");
                }
                else
                {
                    ut.failure("number of moments is NOT correct");
                }

                // Test that mean and flux are correct
        
                {
                    vector<Ordinate> const &ordinates = ordinate_space->ordinates();
                    unsigned const N = ordinates.size();
                    double J = 0.0;
                    double Fx = 0.0, Fy=0.0;
                    double Fx2 = 0.0, Fy2=0.0, Fz2=0.0;
                    double const MAGIC = 2.32; // avoid numerical coincidences
                    cout << "Ordinates:" << endl;
                    for (unsigned i=0; i<N; ++i)
                    {
                        cout << "  mu = "
                             << setprecision(10) << ordinates[i].mu()<< "  eta = "
                             << setprecision(10) << ordinates[i].eta()<< "  xi = "
                             << setprecision(10) << ordinates[i].xi() << " weight = "
                             << setprecision(10) << ordinates[i].wt() << endl;
            
                        J += MAGIC * ordinates[i].wt();
                        Fx += MAGIC * ordinates[i].mu()*ordinates[i].wt();
                        Fx2 += MAGIC * ordinates[i].mu()*ordinates[i].mu()*ordinates[i].wt();
                        Fy += MAGIC * ordinates[i].eta()*ordinates[i].wt();
                        Fy2 += MAGIC * ordinates[i].eta()*ordinates[i].eta()*ordinates[i].wt();
                        Fz2 += MAGIC * ordinates[i].xi()*ordinates[i].xi()*ordinates[i].wt();
                    }
                    if (soft_equiv(J, MAGIC))
                    {
                        ut.passes("J okay");
                    }
                    else
                    {
                        ut.failure("J NOT okay");
                    }
                    if (soft_equiv(Fx, 0.0))
                    {
                        ut.passes("xF okay");
                    }
                    else
                    {
                        ut.failure("Fx NOT okay");
                    }
                    if (soft_equiv(Fx2, MAGIC/3.0))
                    {
                        ut.passes("Fx2 okay");
                    }
                    else
                    {
                        cout << "Fx2 = " << Fx2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fx2 NOT okay");
                    }
                    if (soft_equiv(Fy, 0.0))
                    {
                        ut.passes("Fy okay");
                    }
                    else
                    {
                        ut.failure("Fy NOT okay");
                    }
                    if (soft_equiv(Fy2, MAGIC/3.0))
                    {
                        ut.passes("Fy2 okay");
                    }
                    else
                    {
                        cout << "Fy2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fy2 NOT okay");
                    }
                    if (soft_equiv(Fz2, MAGIC/3.0))
                    {
                        ut.passes("Fz2 okay");
                    }
                    else
                    {
                        cout << "Fz2 = " << Fz2 << ", expected " << (MAGIC/3.0) << endl;
                        ut.failure("Fz2 NOT okay");
                    }
                }
            }
        }
    }
#endif
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//              end of quadrature/quadrature_test.cc
//---------------------------------------------------------------------------//
