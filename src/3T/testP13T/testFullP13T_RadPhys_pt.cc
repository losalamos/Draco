#include "3T/testP13T/testFullP13T.hh"

#include "traits/ContainerTraits.hh"

#include "radphys/RadiationPhysics.t.cc"

using XTM::ContainerTraits;
using XTM::RadiationPhysics;

typedef Mesh_XYZ MT;

template<>
class ContainerTraits<MT::ccsf >
{
  public:
    typedef MT::ccsf::iterator iterator;
    typedef MT::ccsf::const_iterator const_iterator;
    static inline iterator begin(MT::ccsf &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const MT::ccsf &a)
    {
	return a.begin();
    }
    static inline iterator end(MT::ccsf &a)
    {
	return a.end();
    }
    static inline const_iterator end(const MT::ccsf &a)
    {
	return a.end();
    }
    static inline bool conformal(const MT::ccsf &a,
				 const MT::ccsf &b)
    {
	const MT &msha = a.get_Mesh();
	const MT &mshb = b.get_Mesh();

	return msha == mshb;
    }
};

template<>
class ContainerTraits<MT::fcdsf >
{
  public:
    typedef MT::fcdsf::iterator iterator;
    typedef MT::fcdsf::const_iterator const_iterator;
    static inline iterator begin(MT::fcdsf &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const MT::fcdsf &a)
    {
	return a.begin();
    }
    static inline iterator end(MT::fcdsf &a)
    {
	return a.end();
    }
    static inline const_iterator end(const MT::fcdsf &a)
    {
	return a.end();
    }
    static inline bool conformal(const MT::fcdsf &a,
				 const MT::fcdsf &b)
    {
	const MT &msha = a.get_Mesh();
	const MT &mshb = b.get_Mesh();

	return msha == mshb;
    }
};

template void RadiationPhysics::getPlanck(const MT::ccsf &TElectron,
					  MT::ccsf &planckian) const;
template void RadiationPhysics::
    getPlanckTemperatureDerivative(const MT::ccsf &TElectron,
				   MT::ccsf &dplanckdT) const;

