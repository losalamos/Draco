/*-----------------------------------*-C-*-----------------------------------*/
/* UserVecTraits.hh */
/* Randy M. Roberts */
/* Tue Apr 20 15:03:30 1999 */
/*---------------------------------------------------------------------------*/
/* @> Specialized Expression Engine Traits for the UserVec class */
/*---------------------------------------------------------------------------*/

#ifndef __ExpEngineTraits_UserVecTraits_hh__
#define __ExpEngineTraits_UserVecTraits_hh__

#include "UserVec.hh"
#include "expTraits.hh"

namespace rtt_expTraits
{

template<class T>
class ExpEngineTraits< UserVec<T> >
{
  public:
    typedef UserVec<T>  ExpEnabledContainer;

    static UserVec<T> &Glom(UserVec<T> &rct)
    {
	return (rct);
    }
};

} // end namespace rtt_expTraits

#endif    /* __ExpEngineTraits_UserVecTraits_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of UserVecTraits.hh */
/*---------------------------------------------------------------------------*/
