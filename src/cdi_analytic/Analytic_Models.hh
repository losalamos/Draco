//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Models.hh
 * \author Thomas M. Evans
 * \date   Wed Aug 29 16:46:52 2001
 * \brief  Analytic_Model definitions
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_Models_hh__
#define __cdi_analytic_Analytic_Models_hh__

#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_cdi_analytic
{

//===========================================================================//
/*!
 * \class Analytic_Opacity_Model
 * \brief Analytic_Opacity_Model base class.
 *
 * This is a base class that defines the interface given to
 * Analytic_Gray_Opacity constructors.  The user can define any derived class
 * to give to Analytic_Gray_Opacity as long as it contains the following
 * function: (declared virtual in this class).
 *
 * This class is a pure virtual base class. 
 *
 * \arg double calculate_opacity(double T, double rho)
 *
 * The returned opacity should have units of cm^2/g.
 *
 */
//===========================================================================//

class Analytic_Opacity_Model
{
  public:
    //! Virtual destructor for proper inheritance destruction.
    ~Analytic_Opacity_Model() {/*...*/}

    //! Interface for derived analytic opacity models.
    virtual double calculate_opacity(double T, double rho) const = 0;
};

//---------------------------------------------------------------------------//
/*!
 * \class Constant_Analytic_Opacity_Model
 * \brief Derived Analytic_Opacity_Model class that defines a constant
 *  opacity.
 *
 * The opacity is defined:
 *
 * \arg opacity = a
 *
 * where the coefficient has the following units:
 *
 * \arg a = [cm^2/g]
 *
 */
class Constant_Analytic_Opacity_Model : public Analytic_Opacity_Model
{
  private:
    // Constant opacity.
    double sigma;

  public:
    //! Constructor, sig has units of cm^2/g.
    Constant_Analytic_Opacity_Model(double sig) 
	: sigma(sig) 
    { 
	Require (sigma > 0.0); 
    }

    //! Calculate the opacity in units of c^2/g.
    double calculate_opacity(double T, double rho) const
    {
	return sigma;
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class Polynomial_Analytic_Opacity_Model
 * \brief Derived Analytic_Opacity_Model class that defines a polymomial
 * function for the opacity.
 *
 * The opacity is defined:
 *
 * \arg opacity = (a + b * T^c) * d * rho^e
 *
 * where the coefficients have the following units:
 *
 * \arg a = [cm^2/g]
 * \arg b = [keV^(-c) cm^2/g]
 * \arg d = [cm^3/g]
 *
 */
class Polynomial_Analytic_Opacity_Model : public Analytic_Opacity_Model
{
  private: 
    // Coefficients
    double a;  // constant [cm^2/g]
    double b;  // temperature multiplier [keV^(-c) cm^2/g]
    double c;  // temperature power
    double d;  // density multiplier [cm^3/g]
    double e;  // density power

  public:
    /*!
     * \brief Constructor.
     * \param a_ constant [cm^2/g]
     * \param b_ temperature multiplier [keV^(-c) cm^2/g]
     * \param c_ temperature power
     * \param d_ density multiplier [cm^3/g]
     * \param e_ density power
     */
    Polynomial_Analytic_Opacity_Model(double a_, double b_, double c_,
				      double d_, double e_)
	: a(a_), b(b_), c(c_), d(d_), e(e_)
    {
	/*...*/
    }

    //! Calculate the opacity in units of c^2/g
    double calculate_opacity(double T, double rho) const
    {
	Require (T >= 0.0);
	Require (rho >= 0.0);

	double T_power   = std::pow(T, c);
	double rho_power = std::pow(rho, e);

	double opacity   = (a + b * T_power) * d * rho_power;

	Ensure (opacity >= 0.0);
	return opacity;
    }
};

} // end namespace rtt_cdi_analytic

#endif                          // __cdi_analytic_Analytic_Models_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_analytic/Analytic_Models.hh
//---------------------------------------------------------------------------//
