//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Tri_Chebyshev_Legendre.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------------------//
// $Id: Tri_Chebyshev_Legendre.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_Tri_Chebyshev_Legendre_hh
#define quadrature_Tri_Chebyshev_Legendre_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature
{

//=======================================================================================//
/*!
 * \class Tri_Chebyshev_Legendre
 * \brief A class to encapsulate a triangular Chebyshev-Legendre quadrature set.
 */
//=======================================================================================//

class Tri_Chebyshev_Legendre : public Octant_Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    explicit Tri_Chebyshev_Legendre( unsigned sn_order)
        :
        sn_order_( sn_order )
        
    {
        Require(sn_order>0 && sn_order%2==0);
    }

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    explicit Tri_Chebyshev_Legendre( unsigned sn_order,
                                     unsigned const mu_axis,
                                     unsigned const eta_axis)
        :
        Octant_Quadrature(mu_axis, eta_axis),
        sn_order_( sn_order )
        
    {
        Require(sn_order>0 && sn_order%2==0);
    }

    Tri_Chebyshev_Legendre();    // disable default construction

    // ACCESSORS
    
    unsigned sn_order()     const { return sn_order_; }

    // SERVICES
    
    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.

    string name()        const;
    
    string parse_name()  const;
    
    Quadrature_Class quadrature_class() const;

    unsigned number_of_levels() const;
    
    string as_text(string const &indent) const;
    
    // STATICS

    static SP<Quadrature> parse(Token_Stream &tokens);

  private:

    // IMPLEMENTATION
    
    //! Virtual hook for create_ordinate_set
    virtual void create_octant_ordinates_(vector<double> &mu,
                                          vector<double> &eta,
                                          vector<double> &wt) const;


    // DATA
    unsigned sn_order_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Tri_Chebyshev_Legendre_hh

//---------------------------------------------------------------------------------------//
//              end of quadrature/Tri_Chebyshev_Legendre.hh
//---------------------------------------------------------------------------------------//
