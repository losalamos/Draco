//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   P1Diffusion/P1Momentum.t.hh
 * \author Randy M. Roberts
 * \date   Wed Feb  2 09:26:01 2000
 * \brief  Handles the velocity dependent parts of P1Diffusion
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "P1Momentum.hh"
#include "ds++/SP.hh"
#include "traits/MT_traits.hh"

namespace rtt_P1Diffusion
{
 
template<class MT,bool HV>
void P1Momentum<MT,HV>::discFluxToDiscMomentum(DiscMomentumField &result,
					       const DiscFluxField &flux) const
{
    // This method moves the flux-like field from the DiscFluxField location
    // to the DiscMomentumField location.

    typedef typename MT::fcdvsf NormalsField;

    NormalsField faceNormals(fCtor);
    faceNormals.get_Mesh().get_face_normals(faceNormals);

    // ConnFacesAroundVertices is a class that defines objects
    // that will iterate through a face-centered field around each vertex,
    // before going onto the next vertex's faces.

    typedef typename MT::ConnFacesAroundVertices<NormalsField> ConnNormals;
    typedef typename MT::ConnFacesAroundVertices<const DiscFluxField> ConnFlux;

    const ConnNormals connNormals(faceNormals);
    const ConnFlux connFlux(flux);

    // Make sure that everyone is the same size (in vertices).

    Assert(std::distance(connNormals.begin(), connNormals.end()) ==
	   std::distance(connFlux.begin(), connFlux.end()));
    Assert(std::distance(connNormals.begin(), connNormals.end()) ==
	   std::distance(result.begin(), result.end()));
     
    // Loop over vertices
     
    ConnNormals::const_iterator itNormVertex = connNormals.begin();
    ConnFlux::const_iterator itFluxVertex = connFlux.begin();
    DiscMomentumField::iterator itResult = result.begin();

    while (itResult != result.end())
    {
	// Make sure that everyone is the same size (in faces per vertex).
	 
	Assert(std::distance((*itNormVertex).begin(),
			     (*itNormVertex).end()) ==
	       std::distance((*itFluxVertex).begin(),
			     (*itFluxVertex).end()));
	 
	// Loop over faces per vertex
	 
	ConnNormals::value_type::const_iterator itNormFace =
	    (*itNormVertex).begin();
	ConnFlux::value_type::const_iterator itFluxFace =
	    (*itFluxVertex).begin();

	while (itNormFace != (*itNormVertex).end())
	{
	    typedef NormalsField::value_type NormVector;
	    typedef DiscMomentumField::value_type ResultVector;
	    typedef DiscFluxField::value_type FluxType;

	    ResultVector &res = *itResult;
	    const NormVector &norm = *itNormFace;
	    const FluxType &flux = *itFluxFace;

	    res[0] = norm[0]*flux;
	    res[1] = norm[1]*flux;
	    res[2] = norm[2]*flux;
	     
	    itNormFace++;
	    itFluxFace++;
	}

	itResult++;
	itNormVertex++;
	itFluxVertex++;
    }
}

template<class MT,bool HV>
void P1Momentum<MT,HV>::dotProduct(DiscKineticEnergyField &result,
				   const DiscMomentumField &vec1,
				   const DiscMomentumField &vec2) const
{
    Assert(result.size() == vec1.size());
     
    DiscMomentumField::const_iterator iv1 = vec1.begin();
    DiscMomentumField::const_iterator iv2 = vec2.begin();
    DiscKineticEnergyField::iterator ir = result.begin();

    while (ir != result.end())
    {
	typedef DiscMomentumField::value_type vector;	 
	typedef rtt_traits::vector_traits<vector> vtraits;

	*ir = vtraits::dot(*iv1, *iv2);

	ir++;
	iv1++;
	iv2++;
    }
}

template<class MT,bool HV>
void P1Momentum<MT,HV>::dotProduct(DiscKineticEnergyField &KEnergy,
				   const DiscFluxField &sigmaF,
				   const DiscMomentumField &velocity) const
{
    Assert(KEnergy.size() == velocity.size());

    // Move the flux-like field from the DiscFluxField location
    // to the DiscMomentumField location.
     
    DiscMomentumField sigmaFAtMomentum(fCtor);
    discFluxToDiscMomentum(sigmaFAtMomentum, sigmaF);

    // KEnergy is the
    // dot_product((sigma * Flux), velocity)
    // at each vertex.

    dotProduct(KEnergy, sigmaFAtMomentum, velocity);
     
}

} // end namespace rtt_P1Diffusion

//---------------------------------------------------------------------------//
//                        end of P1Diffusion/P1Momentum.t.hh
//---------------------------------------------------------------------------//
