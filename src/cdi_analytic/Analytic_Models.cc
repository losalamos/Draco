//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Models.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 21 14:36:15 2001
 * \brief  Analytic_Models implementation file.
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Analytic_Models.hh"
#include "ds++/Packing_Utils.hh"

namespace rtt_cdi_analytic
{

//===========================================================================//
// CONSTANT_ANALYTIC_MODEL MEMBER DEFINITIONS
//===========================================================================//
// Unpacking constructor.

Constant_Analytic_Opacity_Model::Constant_Analytic_Opacity_Model(
    const sf_char &packed)
    : sigma(0)
{
    // size of stream
    int size( sizeof(int) + sizeof(double) );

    Require (packed.size() == static_cast<size_t>(size));

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    
    // set the unpacker
    unpacker.set_buffer(size, &packed[0]);

    // unpack the indicator
    int indicator;
    unpacker >> indicator;
    Insist (indicator == CONSTANT_ANALYTIC_OPACITY_MODEL,
	    "Tried to unpack the wrong type in Constant_Analytic_Opacity_Model");
	
    // unpack the data
    unpacker >> sigma;
    Check (sigma >= 0.0);

    Ensure (unpacker.get_ptr() == unpacker.end());
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_Opacity_Model::sf_char
Constant_Analytic_Opacity_Model::pack() const 
{
    // get the registered indicator 
    int indicator = CONSTANT_ANALYTIC_OPACITY_MODEL;

    // caculate the size in bytes: indicator + 1 double
    int size = sizeof(int) +  sizeof(double);

    // make a vector of the appropriate size
    sf_char pdata(size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[0]);

    // pack the indicator
    packer << indicator;
	
    // pack the data
    packer << sigma;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
// Return the model parameters

Analytic_Opacity_Model::sf_double
Constant_Analytic_Opacity_Model::get_parameters() const
{
    return sf_double(1, sigma);
}

//===========================================================================//
// POLYNOMIAL_ANALYTIC_OPACITY_MODEL DEFINITIONS
//===========================================================================//
// Unpacking constructor.

Polynomial_Analytic_Opacity_Model::Polynomial_Analytic_Opacity_Model(
    const sf_char &packed)
    : a(0.0), b(0.0), c(0.0), d(0.0), e(0.0)
{
    // size of stream
    size_t size = sizeof(int) + 4 * sizeof(double);

    Require (packed.size() == size);

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    
    // set the unpacker
    unpacker.set_buffer(size, &packed[0]);

    // unpack the indicator
    int indicator;
    unpacker >> indicator;
    Insist (indicator == POLYNOMIAL_ANALYTIC_OPACITY_MODEL,
	    "Tried to unpack the wrong type in Polynomial_Analytic_Opacity_Model");
	
    // unpack the data
    unpacker >> a >> b >> c >> d;

    Ensure (unpacker.get_ptr() == unpacker.end());
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_Opacity_Model::sf_char
Polynomial_Analytic_Opacity_Model::pack() const 
{
    // get the registered indicator 
    int indicator = POLYNOMIAL_ANALYTIC_OPACITY_MODEL;

    // caculate the size in bytes: indicator + 4 * double
    int size = sizeof(int) + 4 * sizeof(double);

    // make a vector of the appropriate size
    sf_char pdata(size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[0]);

    // pack the indicator
    packer << indicator;
	
    // pack the data
    packer << a;
    packer << b;
    packer << c;
    packer << d;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
// Return the model parameters

Analytic_Opacity_Model::sf_double
Polynomial_Analytic_Opacity_Model::get_parameters() const
{
    sf_double p(4);
    p[0] = a;
    p[1] = b;
    p[2] = c;
    p[3] = d;

    return p;
}

//===========================================================================//
// POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL DEFINITIONS
//===========================================================================//
// Unpacking constructor.

Polynomial_Specific_Heat_Analytic_EoS_Model::
Polynomial_Specific_Heat_Analytic_EoS_Model(const sf_char &packed)
    : a(0.0), b(0.0), c(0.0), d(0.0), e(0.0), f(0.0)
{
    // size of stream
    size_t size = sizeof(int) + 6 * sizeof(double);

    Require (packed.size() == size);

    // make an unpacker
    rtt_dsxx::Unpacker unpacker;
    
    // set the unpacker
    unpacker.set_buffer(size, &packed[0]);

    // unpack the indicator
    int indicator;
    unpacker >> indicator;
    Insist (indicator == POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL,
	    "Tried to unpack the wrong type in Polynomial_Specific_Heat_Analytic_EoS_Model");
	
    // unpack the data
    unpacker >> a >> b >> c >> d >> e >> f;

    Ensure (unpacker.get_ptr() == unpacker.end());
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_Opacity_Model::sf_char
Polynomial_Specific_Heat_Analytic_EoS_Model::pack() const 
{
    // get the registered indicator 
    int indicator = POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL;

    // caculate the size in bytes: indicator + 6 * double
    int size = sizeof(int) + 6 * sizeof(double);

    // make a vector of the appropriate size
    sf_char pdata(size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[0]);

    // pack the indicator
    packer << indicator;
	
    // pack the data
    packer << a;
    packer << b;
    packer << c;
    packer << d;
    packer << e;
    packer << f;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
// Return the model parameters

Analytic_EoS_Model::sf_double
Polynomial_Specific_Heat_Analytic_EoS_Model::get_parameters() const
{
    sf_double p(6);
    p[0] = a;
    p[1] = b;
    p[2] = c;
    p[3] = d;
    p[4] = e;
    p[5] = f;

    return p;
}


//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_Opacity_Model::
Pseudo_Line_Analytic_Opacity_Model(double continuum,
                                   unsigned number_of_lines,
                                   double peak,
                                   double width,
                                   double emin,
                                   double emax,
                                   unsigned seed)
    :
    continuum_(continuum),
    seed_(seed),
    number_of_lines_(number_of_lines),
    peak_(peak),
    width_(width),
    center_(number_of_lines)
{
    Require(continuum>=0.0);
    Require(peak>=0.0);
    Require(width>=0.0);
    Require(emin>=0.0);
    Require(emax>emin);
    
    srand(seed);

    for (unsigned i=0; i<number_of_lines; ++i)
    {
        center_[i] = (emax-emin)*static_cast<double>(rand())/RAND_MAX + emin;
    }
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_Opacity_Model::sf_char
Pseudo_Line_Analytic_Opacity_Model::pack() const 
{
    // get the registered indicator 
    int indicator = PSEUDO_LINE_ANALYTIC_OPACITY_MODEL;

    // caculate the size in bytes: indicator + 4 * double
    int size =
        sizeof(int) + 3 * sizeof(double) + 2 * sizeof(int);

    // make a vector of the appropriate size
    sf_char pdata(size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[0]);

    // pack the indicator
    packer << indicator;
	
    // pack the data
    packer << continuum_;
    packer << seed_;
    packer << number_of_lines_;
    packer << peak_;
    packer << width_;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
double Pseudo_Line_Analytic_Opacity_Model::calculate_opacity(double T,
                                                             double rho,
                                                             double nu) const
{
    double Result = continuum_;
    for (unsigned i=0; i<number_of_lines_; ++i)
    {
        double const nu0 = center_[i];
        double const d = nu - nu0;
        Result += peak_*exp(-d*d/(width_*width_*nu0*nu0));
    }

    return Result;
}

//---------------------------------------------------------------------------//
double Pseudo_Line_Analytic_Opacity_Model::calculate_opacity(double T,
                                                             double rho) const
{
    return continuum_;
    // Really shouldn't be used as a gray model
}

//---------------------------------------------------------------------------//
// Return the model parameters

Analytic_Opacity_Model::sf_double
Pseudo_Line_Analytic_Opacity_Model::get_parameters() const
{
    sf_double p(5);
    p[0] = continuum_;
    p[1] = seed_;
    p[2] = number_of_lines_;
    p[3] = peak_;
    p[4] = width_;

    return p;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Analytic_Models.cc
//---------------------------------------------------------------------------//
