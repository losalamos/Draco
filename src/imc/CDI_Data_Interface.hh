//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/CDI_Data_Interface.hh
 * \author Thomas M. Evans
 * \date   Fri Nov 16 11:56:40 2001
 * \brief  Interface definition for CDI_Mat_State_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_CDI_Data_Interface_hh__
#define __imc_CDI_Data_Interface_hh__

#include "cdi/CDI.hh"
#include "ds++/SP.hh"
#include <vector>
#include <utility>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class CDI_Data_Interface
 *
 * \brief Additional interface functions required by the
 * CDI_Mat_State_Builder.
 *
 * This interface class contains additional definitions (beyond those
 * required by rtt_imc::Interface) required by the CDI_Mat_State_Builder.
 * To make a class an interface for the CDI_Mat_State_Builder it must
 * inherit from both rtt_imc::Interface and rtt_imc::CDI_Data_Interface.
 *
 * There are no name collisions between this interface definition,
 * Flat_Data_Interface, and rtt_imc::Interface.
 *
 * The interface data required by CDI_Mat_State_Builder is:
 * - a list of CDI objects for each material in the problem
 * - a cell-length map that maps the CDI objects to each cell
 * - a list of pairs that gives the (absorption,scattering) model used for
 * each CDI object
 * .
 * For example, assume a 4 cell mesh with 2 materials (hydrogen and
 * aluminum).  The CDI list has the following form:
 * \code
 *     cdi_list[0] = cdi_Al;
 *     cdi_list[1] = cdi_H;
 * \endcode
 * If cells 1 and 3 contain aluminum and cells 2 and 4 contain hydrogen, then
 * the cell map should have the following form:
 * \code
 *     cdi_map[0] = 1;
 *     cdi_map[1] = 2;
 *     cdi_map[2] = 1;
 *     cdi_map[3] = 2;
 * \endcode
 * The cdi indices returned by cdi_map should be numbered [1,num_materials].
 * Thus, the CDI object for variable cell, where cell > 0 and cell <=
 * num_cells, is cdi_list[cdi_map[cell-1]-1].
 *
 * Finally, a list of pairs is given that tells what model to use for
 * absorption and scattering.  In our example above, if we wished to
 * rtt_cdi::ROSSELAND for absorption and rtt_cdi::ANALYTIC for scattering for
 * the Al and rtt_cdi::PLANCK for absorption for rtt_cdi::ROSSELAND for
 * scattering we would have:
 * \code
 *     cdi_models[0].first  = rtt_cdi::ROSSELAND;
 *     cdi_models[0].second = rtt_cdi::ANALYTIC;
 *     cdi_models[1].first  = rtt_cdi::PLANCK;
 *     cdi_models[1].second = rtt_cdi::ROSSELAND;  
 * \endcode
 * Thus, the first member of the pair is the absorption model and the second
 * member is the scattering model.
 */
// revision history:
// -----------------
// 0) original
// 1) 21-FEB-2003 : added cdi models description
// 
//===========================================================================//

class CDI_Data_Interface 
{
  private:
    // Useful typedefs.
    typedef std::pair<int, int>        model_pair;
    typedef rtt_dsxx::SP<rtt_cdi::CDI> SP_CDI;

  public:
    //! Constructor.
    CDI_Data_Interface() { /* no data to construct */ }
  
    //! Virtual constructor to make life happy down the inhertiance chain.
    virtual ~CDI_Data_Interface() { }

    // >>> PUBLIC INTERFACE REQUIRED BY CDI_MAT_STATE_BUILDER
    
    //! Get a list of CDIs for each material in the problem.
    virtual std::vector<SP_CDI>     get_CDIs()       const = 0;
 
    //! Get a cell-length field mapping CDI objects to cells.
    virtual std::vector<int>        get_CDI_map()    const = 0;

    //! Get a list of pairs giving (absorption,scattering) models.
    virtual std::vector<model_pair> get_CDI_models() const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_CDI_Data_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/CDI_Data_Interface.hh
//---------------------------------------------------------------------------//
