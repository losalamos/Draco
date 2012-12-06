//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Double_Gauss.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval double Gauss-Legendre quadrature set.
 * \note   Copyright 2000-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef __quadrature_Double_Gauss_hh__
#define __quadrature_Double_Gauss_hh__

#include "Interval_Quadrature.hh"

namespace rtt_quadrature
{

//=======================================================================================//
/*!
 * \class Double_Gauss
 *
 * \brief A class representing an interval double Gauss-Legendre quadrature set.
 *
 * This is an interval (e.g. 1D) angle quadrature set in which the
 * Gauss-Legendre quadrature of order n/2 is mapped separately onto the
 * intervals [-1,0] and [0,1] assuming N>2. The case N=2 does not integrate
 * the flux and we substitute ordinary Gauss-Legendre quadrature (N=2) on the
 * full interval [-1,1].
 */
//=======================================================================================//

class Double_Gauss : public Interval_Quadrature 
{
  public:
    
    // CREATORS

    explicit Double_Gauss(unsigned sn_order)
        : sn_order_(sn_order) 
    {
        Require(sn_order>0 && sn_order%2==0);

        Ensure(check_class_invariants());
        Ensure(this->sn_order()==sn_order);
    }

    // ACCESSORS

    unsigned sn_order() const
    {
        return sn_order_;
    }

    // SERVICES

    virtual string name() const;

    virtual string parse_name() const;

    virtual unsigned number_of_levels() const;

    virtual string as_text(string const &indent) const;

    bool check_class_invariants() const;

    // STATICS

    static SP<Quadrature> parse(Token_Stream &tokens);

  protected:
    
    virtual vector<Ordinate> create_level_ordinates_(double norm) const;

    // DATA

    unsigned sn_order_;
};

} // end namespace rtt_quadrature

#endif // __quadrature_Quadrature_hh__

//---------------------------------------------------------------------------------------//
//                       end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------------------//
