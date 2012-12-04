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

#include "quadrature_test.hh"

#include "parser/String_Token_Stream.hh"

namespace rtt_quadrature
{
using namespace std;
using namespace rtt_parser;

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

    // Build an angle operator

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
    
    // **** Try a 3-D if it's an octant quadrature.
        
    if (quadrature.quadrature_class() ==  Quadrature::OCTANT_QUADRATURE)
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
        }
    }
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//              end of quadrature/quadrature_test.cc
//---------------------------------------------------------------------------//
