//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Models.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 21 14:36:15 2001
 * \brief  Analytic_Models implementation file.
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
    int size = sizeof(int) + sizeof(double);

    Require (packed.size() == size);

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
{
    // size of stream
    int size = sizeof(int) + 5 * sizeof(double);

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
    unpacker >> a >> b >> c >> d >> e;

    Ensure (unpacker.get_ptr() == unpacker.end());
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_Opacity_Model::sf_char
Polynomial_Analytic_Opacity_Model::pack() const 
{
    // get the registered indicator 
    int indicator = POLYNOMIAL_ANALYTIC_OPACITY_MODEL;

    // caculate the size in bytes: indicator + 5 * double
    int size = sizeof(int) + 5 * sizeof(double);

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

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
// Return the model parameters

Analytic_Opacity_Model::sf_double
Polynomial_Analytic_Opacity_Model::get_parameters() const
{
    sf_double p(5);
    p[0] = a;
    p[1] = b;
    p[2] = c;
    p[3] = d;
    p[4] = e;

    return p;
}

//===========================================================================//
// POLYNOMIAL_SPECIFIC_HEAT_ANALYTIC_EOS_MODEL DEFINITIONS
//===========================================================================//
// Unpacking constructor.

Polynomial_Specific_Heat_Analytic_EoS_Model::
Polynomial_Specific_Heat_Analytic_EoS_Model(const sf_char &packed)
{
    // size of stream
    int size = sizeof(int) + 6 * sizeof(double);

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

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Analytic_Models.cc
//---------------------------------------------------------------------------//
