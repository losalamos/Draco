//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices_SN.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide Moment-to-Discrete and Discrete-to-Moment operations.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id: QuadServices_SN.hh 6499 2012-03-15 20:19:33Z kgbudge $
//---------------------------------------------------------------------------//

#ifndef quadrature_QuadServices_SN_hh
#define quadrature_QuadServices_SN_hh

#include "ds++/SP.hh"
#include "Ordinate.hh"
#include "Quadrature.hh"
#include "QuadServices.hh"

namespace rtt_quadrature
{

class QuadServices_SN : public QuadServices
{
  public:
    
    // NESTED CLASSES AND TYPEDEFS
    
    // CREATORS
    
    //! Default constructor assumes that only isotropic scattering is used. 
    QuadServices_SN( std::vector<Ordinate> const &ordinates,
                     double                const norm,
                     unsigned              const dimension,
                     unsigned              const expansionOrder = 0,
                     rtt_mesh_element::Geometry const geometry = rtt_mesh_element::CARTESIAN);
    
    //! Copy constructor (the long doxygen description is in the .cc file).
    QuadServices_SN( QuadServices_SN const & rhs );

    //! Destructor.
    virtual ~QuadServices_SN() { /* empty */ }

    // MANIPULATORS
    
    //! Assignment operator for QuadServices_SN.
    QuadServices_SN& operator=( QuadServices_SN const & rhs );

    // ACCESSORS

    //! \brief Provide the maximum available expansion order.
    virtual unsigned getMaxExpansionOrder() const { return maxExpansionOrder; }

    //! \brief Provide the number of moments used by QuadServices_SN.
    virtual unsigned getNumMoments() const { return n2lk.size(); }

    //! \brief Provide the number of moments in each expansion order
    virtual std::vector<unsigned> getMoments() const { return moments; }

    //! \brief Provide the number of moments used by QuadServices_SN.
    virtual std::vector< lk_index > const &get_n2lk() const { return n2lk; }

    //! \brief Return the (l,k) index pair associated with moment index n.
    virtual lk_index lkPair( unsigned n ) const
    {
        Check(n<n2lk.size());
        return n2lk[n];
    }

    //! \brief Return the moment-to-discrete operator.
    virtual std::vector< double > getM_() const { return M; }

    //! \brief Return the discrete-to-moment operator.
    virtual std::vector< double > getD_() const { return D; }

  private:

    std::vector< lk_index > compute_n2lk( unsigned const expansionOrder,
                                          unsigned const dim);
    
    std::vector< double > computeM();    
    std::vector< double > computeM(std::vector<Ordinate> const &ordinates,
                                   std::vector< lk_index > const &n2lk,
                                   unsigned const dim,
                                   double const sumwt);
    
    std::vector< double > computeD();
    std::vector< double > computeD(std::vector<Ordinate> const &ordinates,
                                   std::vector< lk_index > const &n2lk,
                                   unsigned const dim,
                                   double const sumwt);

    std::vector< double > computeD(std::vector<Ordinate> const &ordinates,
                                   std::vector< lk_index > const &n2lk,
                                   std::vector<double> const &Mm);
    
    static std::vector< lk_index > compute_n2lk_2D( unsigned const L);
    static std::vector< lk_index > compute_n2lk_3D( unsigned const L);

    // DATA
    rtt_mesh_element::Geometry const geometry;
    std::vector< lk_index >   const n2lk;
    unsigned                  const maxExpansionOrder;
    std::vector< unsigned >   const moments;
    std::vector< double >     const M;
    std::vector< double >     const D;
};

} // end namespace rtt_quadrature

#endif // quadrature_QuadServices_SN_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices_SN.hh
//---------------------------------------------------------------------------//
