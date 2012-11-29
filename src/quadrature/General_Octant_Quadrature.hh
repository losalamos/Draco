//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/General_Octant_Quadrature.hh
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  A class to encapsulate a 3D Level Symmetric quadrature set.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------------------//
// $Id: General_Octant_Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_General_Octant_Quadrature_hh
#define quadrature_General_Octant_Quadrature_hh

#include "Octant_Quadrature.hh"

namespace rtt_quadrature
{

//=======================================================================================//
/*!
 * \class General_Octant_Quadrature
 * \brief A class to encapsulate a client-defined ordinate set.
 */
//=======================================================================================//

class General_Octant_Quadrature : public Octant_Quadrature
{
  public:

    // CREATORS

    // The default values for snOrder_ and norm_ were set in QuadCreator.
    explicit General_Octant_Quadrature(vector<double> const &mu,
                                       vector<double> const &eta,
                                       vector<double> const &xi,
                                       vector<double> const &wt,
                                       unsigned number_of_levels,
                                       QIM const qim)
        :
        Octant_Quadrature(qim),
        mu_(mu), eta_(eta), xi_(xi), wt_(wt), number_of_levels_(number_of_levels)
        
    {
        Require(mu.size()>0 && eta.size()==mu.size() && xi.size()==mu.size() &&
                wt.size()==mu.size());
    }

    General_Octant_Quadrature();    // disable default construction

    // ACCESSORS
    
    vector<double> const & mu()     const { return mu_; }
    vector<double> const & eta()    const { return eta_; }
    vector<double> const & xi()     const { return xi_; }
    vector<double> const & wt()     const { return wt_; }

    // SERVICES
    
    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.

    string name()        const;
    
    string parse_name()  const;
    
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
    vector<double> mu_, eta_, xi_, wt_;
    unsigned number_of_levels_;
};

} // end namespace rtt_quadrature

#endif // quadrature_General_Octant_Quadrature_hh

//---------------------------------------------------------------------------------------//
//              end of quadrature/General_Octant_Quadrature.hh
//---------------------------------------------------------------------------------------//
