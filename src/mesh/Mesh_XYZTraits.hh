//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/MeshXYZTraits.hh
 * \author Randy M. Roberts
 * \date   Mon Feb  7 17:02:22 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mesh_Mesh_XYZTraits_hh__
#define __mesh_Mesh_XYZTraits_hh__

#include "Mesh_XYZ.hh"
#include "traits/ContainerTraits.hh"

namespace rtt_traits
{
 
//===========================================================================//
/*!
 * \class ContainerTraits
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class T1>
class ContainerTraits<Mesh_XYZ::cctf<T1> >
{
  public:
    typedef typename Mesh_XYZ::cctf<T1>::iterator iterator;
    typedef typename Mesh_XYZ::cctf<T1>::const_iterator const_iterator;
    static inline iterator begin(Mesh_XYZ::cctf<T1> &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const Mesh_XYZ::cctf<T1> &a)
    {
	return a.begin();
    }
    static inline iterator end(Mesh_XYZ::cctf<T1> &a)
    {
	return a.end();
    }
    static inline const_iterator end(const Mesh_XYZ::cctf<T1> &a)
    {
	return a.end();
    }
    static inline bool conformal(const Mesh_XYZ::cctf<T1> &a,
				 const Mesh_XYZ::cctf<T1> &b)
    {
	return a.get_Mesh() == b.get_Mesh();
    }
};

template<class T1, class T2>
class ContainerTraits<Mesh_XYZ::cctf<T1>, Mesh_XYZ::cctf<T2> >
{
    static inline bool conformal(const Mesh_XYZ::cctf<T1> &a,
				 const Mesh_XYZ::cctf<T2> &b)
    {
	return a.get_Mesh() == b.get_Mesh();
    }
};

template<class T1>
class ContainerTraits<Mesh_XYZ::fcdtf<T1> >
{
  public:
    typedef typename Mesh_XYZ::fcdtf<T1>::iterator iterator;
    typedef typename Mesh_XYZ::fcdtf<T1>::const_iterator const_iterator;
    static inline iterator begin(Mesh_XYZ::fcdtf<T1> &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const Mesh_XYZ::fcdtf<T1> &a)
    {
	return a.begin();
    }
    static inline iterator end(Mesh_XYZ::fcdtf<T1> &a)
    {
	return a.end();
    }
    static inline const_iterator end(const Mesh_XYZ::fcdtf<T1> &a)
    {
	return a.end();
    }
    static inline bool conformal(const Mesh_XYZ::fcdtf<T1> &a,
				 const Mesh_XYZ::fcdtf<T1> &b)
    {
	return a.get_Mesh() == b.get_Mesh();
    }
};

template<class T1, class T2>
class ContainerTraits<Mesh_XYZ::fcdtf<T1>, Mesh_XYZ::fcdtf<T2> >
{
    static inline bool conformal(const Mesh_XYZ::fcdtf<T1> &a,
				 const Mesh_XYZ::fcdtf<T2> &b)
    {
	return a.get_Mesh() == b.get_Mesh();
    }
};

template<class T1>
class ContainerTraits<Mesh_XYZ::bstf<T1> >
{
  public:
    typedef typename Mesh_XYZ::bstf<T1>::iterator iterator;
    typedef typename Mesh_XYZ::bstf<T1>::const_iterator const_iterator;
    static inline iterator begin(Mesh_XYZ::bstf<T1> &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const Mesh_XYZ::bstf<T1> &a)
    {
	return a.begin();
    }
    static inline iterator end(Mesh_XYZ::bstf<T1> &a)
    {
	return a.end();
    }
    static inline const_iterator end(const Mesh_XYZ::bstf<T1> &a)
    {
	return a.end();
    }
    static inline bool conformal(const Mesh_XYZ::bstf<T1> &a,
				 const Mesh_XYZ::bstf<T1> &b)
    {
	return a.get_Mesh() == b.get_Mesh();
    }
};

template<class T1, class T2>
class ContainerTraits<Mesh_XYZ::bstf<T1>, Mesh_XYZ::bstf<T2> >
{
    static inline bool conformal(const Mesh_XYZ::bstf<T1> &a,
				 const Mesh_XYZ::bstf<T2> &b)
    {
	return a.get_Mesh() == b.get_Mesh();
    }
};

} // end namespace rtt_traits

#endif                          // __mesh_Mesh_XYZTraits_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Mesh_XYZTraits.hh
//---------------------------------------------------------------------------//
