/*-----------------------------------*-C-*-----------------------------------*/
/* expTraits.hh */
/* Randy M. Roberts */
/* Mon Apr 19 13:14:26 1999 */
/*---------------------------------------------------------------------------*/
/* @> Default Expression Engine Traits */
/*---------------------------------------------------------------------------*/

#ifndef __ExpEngineTraits_expTraits_hh__
#define __ExpEngineTraits_expTraits_hh__

#include "xm/xm.hh"

namespace rtt_expTraits
{
    
template<class RCT>
class ExpEngineTraits
{
    // This nested class gloms on the XM expression templates
    // to the class, RCT (Random Access Container Type)
    // Notice that this class publicly inherits from the original class.
    // This enables the expression enabled class to act like the original
    // class.
    
    class ERCT
	: public xm::Indexable<typename RCT::value_type, ERCT>,
	  public RCT
    {
      public:

	// Typedefs required by random access containers.
	
	typedef typename RCT::pointer          iterator;
	typedef typename RCT::const_pointer    const_iterator;
	typedef typename RCT::reference        reference;
	typedef typename RCT::const_reference  const_reference;
	typedef typename RCT::size_type        size_type;
	typedef typename RCT::difference_type  difference_type;
	typedef typename RCT::value_type       value_type;
	typedef typename RCT::allocator_type   allocator_type;
	typedef typename RCT::pointer          pointer;
	typedef typename RCT::const_pointer    const_pointer;
	typedef typename RCT::const_reverse_iterator     const_reverse_iterator;
	typedef typename RCT::reverse_iterator           reverse_iterator;

      private:
	
	// No Data
	
      public:

	// This method required to disambiguate the bracket operator.
	
        const value_type &operator[]( int i ) const
	{
	    return RCT::operator[](i);
	}
	
	// This method required to disambiguate the bracket operator.
	
        value_type &operator[]( int i )
	{
	    return RCT::operator[](i);
	}

	// This method required for XM to work.
	// Or else you wont be able to say...
	//    ec = val;
	// where ec is an expression enabled container
	
	ERCT &operator=( const value_type &val)
	{
	    std::fill(begin(), end(), val);
	    return *this;
	}

	// This method required for XM to work.
	
        template<class X>
        ERCT &operator=( const xm::Xpr< value_type, X, ERCT >& x )
        {
            return assign_from( x );
        }
    };
    
  public:

    // The typedef that is used in applications.
    
    typedef ERCT ExpEnabledContainer;

    // This static method converts an instance of a container
    // to a reference of an expression enabled container.
    
    static ERCT &Glom(RCT &rct)
    {
	return static_cast<ERCT &>(rct);
    }
};

} // end namespace rtt_expTraits

#endif    /* __ExpEngineTraits_expTraits_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of expTraits.hh */
/*---------------------------------------------------------------------------*/
