#include "3T/testP13T/testFullP13T.hh"

#include "traits/ContainerTraits.hh"

#include "3T/Diffusion_P1.cc"
#include "3T/P13T.cc"
#include "radphys/RadiationPhysics.t.cc"
#include "matprops/InterpedMaterialProps.cc"

using namespace XTM;

typedef Mesh_XYZ MT;
typedef InterpedMaterialProps MP;
typedef Diffusion_P1<MT> DS;

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
	return a.size() == b.size();
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
	return a.size() == b.size();
    }
};

// template class MT;
// template class DS;
template class P13T<MT,MP,DS>;
template void RadiationPhysics::getPlanck(const MT::ccsf &TElectron,
					  MT::ccsf &planckian) const;
template void RadiationPhysics::
    getPlanckTemperatureDerivative(const MT::ccsf &TElectron,
				   MT::ccsf &dplanckdT) const;

typedef MT::ccsf ccsf;
typedef MT::fcdsf fcdsf;

typedef ccsf T1;
typedef ccsf T2;

template
MP::MaterialStateField<T1>
MP::getMaterialState<T1, T2>(const T1 &, const T1 &, const T1 &,
			     const T2 &) const;
typedef fcdsf TT1;
typedef fcdsf TT2;

template
MP::MaterialStateField<TT1>
MP::getMaterialState<TT1, TT2>(const TT1 &, const TT1 &, const TT1 &,
			       const TT2 &) const;
