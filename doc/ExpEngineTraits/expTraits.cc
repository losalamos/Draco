//----------------------------------*-C++-*----------------------------------//
// expTraits.cc
// Randy M. Roberts
// Tue Apr 20 17:01:56 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "expTraits.hh"
#include "UserVec.hh"
#include "UserVecTraits.hh"
#include <vector>
#include <iostream>

// The function doit is templated on the container.
// The container may or may not have expression templates associated
// with it.
// If it does not have expression templates (the default assumption)
// then the ECT::Glom static method will glom Draco's XM expression
//  machinery onto it.
// If it does have expression templates already, and has defined its
// own specialization of ExpEngineTraits (as UserVec has done in
// UserVecTraits.hh) then ECT::Glom will do just about nothing.

template <class Container>
void doit()
{
    using rtt_expTraits::ExpEngineTraits;

    // Create some typedef's to save my typing fingers.
    // ECT is the Expression Engine Traits class for the container
    // EC is the nested type that (optionally) gloms expressions
    // onto the original container class.
    
    typedef ExpEngineTraits< Container > ECT;
    typedef ECT::ExpEnabledContainer EC;

    // Create some containers.

    const int sz = 5;
    
    Container c1(sz), c2(sz), c3(sz);

    // Glom (optionally) expressions onto the containers.
    
    EC &ec1 = ECT::Glom(c1);
    EC &ec2 = ECT::Glom(c2);
    EC &ec3 = ECT::Glom(c3);

    // Prove that we have expressions.
    
    ec1 = 1.0;
    ec2 = 20.0;
    ec1[1] = std::sqrt(2.0)/2.0;

    ec3 = asin(ec1) + 3.0*ec2;

    // *** NOTE *** NOTE *** NOTE ***
    //
    // Note that we are now using the original container and not
    // the (optionally) expression-glommed version in the std::copy.
    //
    // *** NOTE *** NOTE *** NOTE ***
    
    std::copy(c3.begin(), c3.end(),
	      std::ostream_iterator<Container::value_type>(std::cout, "\n"));
}

int main()
{

    // Call doit with std::vector<double>.
    // This container does not have expression templates, and by
    // default will have Draco's XM expression template engine
    // glommed onto it inside doit.
    
    std::cout << "doit< std::vector<double> > ..." << std::endl;
    doit< std::vector<double> >();

    std::cout << std::endl;
    
    // Call doit with UserVec<double>.
    // This container already has expression templates, and
    // since it defines its own specialization of
    // ExpEngineTraits< UserVec<T> > that has been designed
    // to do no glomming, i.e. EC is just Container,
    // and ECT::Glom is a no-op, doit will use the UserVec's own
    // expression template machinery.
    
    std::cout << "doit< UserVec<double> > ..." << std::endl;
    doit< UserVec<double> >();
    
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of expTraits.cc
//---------------------------------------------------------------------------//
