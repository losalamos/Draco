//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "CDI.hh"
#include "GrayOpacity.hh"
#include "MultigroupOpacity.hh"
#include "EoS.hh"
#include "ds++/Assert.hh"

namespace rtt_cdi
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS AND DESTRUCTORS
//---------------------------------------------------------------------------//

CDI::CDI(const std_string &id)
    : matID(id),
      grayOpacities(constants::num_Models, 
		    SF_GrayOpacity(constants::num_Reactions)),
      multigroupOpacities(constants::num_Reactions)
{
    Ensure (grayOpacities.size() == constants::num_Models);
    Ensure (multigroupOpacities.size() == constants::num_Reactions);
}

//---------------------------------------------------------------------------//
    
CDI::~CDI() 
{
    // empty
}
    
//---------------------------------------------------------------------------//
// SET FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Register a gray opacity (rtt_cdi::GrayOpacity) with CDI.
 *
 * This function sets a gray opacity object of type rtt_cdi::GrayOpacity with
 * the CDI object.  It stores the gray opacity object based upon its
 * rtt_cdi::Model and rtt_cdi::Reaction types.  If a GrayOpacity with these
 * types has already been registered an exception is thrown.  To register a
 * new set of GrayOpacity objects call CDI::reset() first.  You cannot
 * overwrite registered objects with the setGrayOpacity() function!
 *
 * \param spGOp smart pointer (rtt_dsxx::SP) to a GrayOpacity object
 */
void CDI::setGrayOpacity(const SP_GrayOpacity &spGOp)
{
    Require (spGOp);

    // determine the model and reaction type (these MUST be in the correct
    // range because the Model and Reaction are constrained by the
    // rtt_cdi::Model and rtt_cdi::Reaction enumerations, assuming nobody
    // hosed these)
    int model    = spGOp->getModelType();
    int reaction = spGOp->getReactionType();

    Insist (!grayOpacities[model][reaction], 
	    "Tried to overwrite a set GrayOpacity object!");

    // assign the smart pointer
    grayOpacities[model][reaction] = spGOp;

    Ensure (grayOpacities[model][reaction]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Register a multigroup opacity (rtt_cdi::MultigroupOpacity) with
 * CDI.
 *
 * This function sets a multigroup opacity object of type
 * rtt_cdi::MultigroupOpacity with the CDI object.  It stores the multigroup
 * opacity object based upon its rtt_cdi::Reaction type.  If a
 * MultigroupOpacity with this type has already been registered an exception
 * is thrown.  To register a new set of MultigroupOpacity objects call
 * CDI::reset() first.  You cannot overwrite registered objects with the
 * setMultigroupOpacity() function!
 *
 * \param spGOp smart pointer (rtt_dsxx::SP) to a MultigroupOpacity object
 */
void CDI::setMultigroupOpacity(const SP_MultigroupOpacity &spMGOp) 
{
    Require (spMGOp);

    // determine the reaction type
    int reaction = spMGOp->getReactionType();

    Insist (!multigroupOpacities[reaction],
	    "Tried to overwrite a set MultigroupOpacity object!");\

    // assign the smart pointer
    multigroupOpacities[reaction] = spMGOp;

    Ensure (multigroupOpacities[reaction]);
}

//---------------------------------------------------------------------------//

void CDI::setEoS( const SP_EoS &in_spEoS )
{
    Require (in_spEoS);

    Insist (!spEoS, "Tried to overwrite a set EoS object.!");

    // set the smart pointer
    spEoS = in_spEoS;

    Ensure (spEoS);
}

//---------------------------------------------------------------------------//
// GET FUNCTIONS
//---------------------------------------------------------------------------//

// Provide CDI with access to the full interfaces defined by
// GrayOpacity.hh and MultigroupOpacity.hh
    
CDI::SP_GrayOpacity CDI::gray(rtt_cdi::Model    m, 
			      rtt_cdi::Reaction r) const 
{ 
    Insist (grayOpacities[m][r], "Undefined GrayOpacity!");
    return grayOpacities[m][r]; 
}
    
CDI::SP_MultigroupOpacity CDI::mg(rtt_cdi::Reaction r) const 
{
    Insist (multigroupOpacities[r], "Undefined MultigroupOpacity!");
    return multigroupOpacities[r]; 
}

//---------------------------------------------------------------------------//

// Provide CDI with access to the full interfaces defined by
// EoS.hh
    
CDI::SP_EoS CDI::eos() const 
{ 
    Insist (spEoS, "Undefined EoS!");
    return spEoS; 
}

//---------------------------------------------------------------------------//
// RESET THE CDI OBJECT
//---------------------------------------------------------------------------//

void CDI::reset()
{
    Check (grayOpacities.size() == constants::num_Models);
    Check (multigroupOpacities.size() == constants::num_Reactions);

    // reset the gray opacities
    for (int i = 0; i < constants::num_Models; i++)
    {
	Check (grayOpacities.size() == constants::num_Reactions);

	for (int j = 0; j < constants::num_Reactions; j++)
	{
	    // reassign the GrayOpacity SP
	    grayOpacities[i][j] = SP_GrayOpacity();

	    // check it
	    Check (!grayOpacities[i][j]);
	}
    }

    // reset the multigroup opacities
    for (int i = 0; i < constants::num_Reactions; i++)
    {
	// reassign the MultigroupOpacity SP
	multigroupOpacities[i] = SP_MultigroupOpacity();

	// check it
	Check (!multigroupOpacities[i]);
    }

    // reset the EoS SP
    spEoS = SP_EoS();
    Check (!spEoS);
}


} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
