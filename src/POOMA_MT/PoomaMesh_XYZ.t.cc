//----------------------------------*-C++-*----------------------------------//
// PoomaMesh_XYZ.t.cc
// Julian C. Cummings
// Fri Sep 25 11:02:08 1998
//---------------------------------------------------------------------------//
// @> Definitions for PoomaMesh_XYZ
//---------------------------------------------------------------------------//


template <class Mesh>
template <class T1, class T2, class Op> 
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
  const typename PoomaMesh_XYZ<Mesh>::cctf<T2>& from, const Op& op)
{
  // scatter to each local face
  PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
  PoomaMesh_XYZ<Mesh>::cctf<T2>::const_iterator fromit, fromend = from.end();
  for (fromit = from.begin(); fromit != fromend; ++fromit)
    for (int f=0; f<6; ++f, ++toit)
      PETE_apply(op, *toit, *fromit);

  // scatter to each adjacent face
  int loc[3];
  int ncx = from.get_Mesh().get_ncx();
  int ncy = from.get_Mesh().get_ncy();
  int ncz = from.get_Mesh().get_ncz();
  toit = to.begin();
  for (fromit = from.begin(); fromit != fromend; ++fromit) {
    fromit.GetCurrentLocation(loc);
    if (loc[0] != 0)
      PETE_apply(op, *toit, fromit.offset(-1,0,0));
    ++toit;
    if (loc[0] != ncx-1)
      PETE_apply(op, *toit, fromit.offset(+1,0,0));
    ++toit;
    if (loc[1] != 0)
      PETE_apply(op, *toit, fromit.offset(0,-1,0));
    ++toit;
    if (loc[1] != ncy-1)
      PETE_apply(op, *toit, fromit.offset(0,+1,0));
    ++toit;
    if (loc[2] != 0)
      PETE_apply(op, *toit, fromit.offset(0,0,-1));
    ++toit;
    if (loc[2] != ncz-1)
      PETE_apply(op, *toit, fromit.offset(0,0,+1));
    ++toit;
  }
}

template <class Mesh>
template <class T1, class T2, class Op> 
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::cctf<T1>& to,
  const typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>& from, const Op& op)
{
  // iterate over fields and scatter from each local face
  PoomaMesh_XYZ<Mesh>::cctf<T1>::iterator toit, toend = to.end();
  PoomaMesh_XYZ<Mesh>::fcdtf<T2>::const_iterator fromit = from.begin();
  for (toit = to.begin(); toit != toend; ++toit)
    for (int f=0; f<6; ++f, ++fromit)
      PETE_apply(op, *toit, *fromit);
}

template <class Mesh>
template <class T1, class T2, class Op> 
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
  const typename PoomaMesh_XYZ<Mesh>::cctf<T2>& from, const Op& op)
{
  // iterate over fields and gather to each local face
  PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
  PoomaMesh_XYZ<Mesh>::cctf<T2>::const_iterator fromit, fromend = from.end();
  for (fromit = from.begin(); fromit != fromend; ++fromit)
    for (int f=0; f<6; ++f, ++toit)
      PETE_apply(op, *toit, *fromit);
}

template <class Mesh>
template <class T1, class T2, class Op> 
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::bstf<T1>& to,
  const typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>& from, const Op& op)
{
  // iterate over boundary faces and gather to bstf
  PoomaMesh_XYZ<Mesh>::bstf<T1>::iterator toit, toend = to.end();
  PoomaMesh_XYZ<Mesh>::fcdtf<T2>::const_iterator fromit = from.begin();
  for (toit = to.begin(); toit != toend; ++toit) {
    fromit = toit;  // set fcdtf iterator to same location
    PETE_apply(op, *toit, *fromit);
  }
}

template <class Mesh>
template <class T1, class T2, class Op> 
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
  const typename PoomaMesh_XYZ<Mesh>::bstf<T2>& from, const Op& op)
{
  // iterate over boundary faces and gather to fcdtf
  PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
  PoomaMesh_XYZ<Mesh>::bstf<T2>::const_iterator fromit, fromend = from.end();
  for (fromit = from.begin(); fromit != fromend; ++fromit) {
    toit = fromit;  // set fcdtf iterator to same location
    PETE_apply(op, *toit, *fromit);
  }
}

template <class Mesh>
template <class T>
void 
PoomaMesh_XYZ<Mesh>::swap(typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& to,
                          const typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& from)
{
  // get BaseField references
  typedef typename PoomaMesh_XYZ<Mesh>::fcdtf<T>::BaseField_t BF_t;
  BF_t& bfto(dynamic_cast<BF_t&>(to));
  const BF_t& bffrom(dynamic_cast<const BF_t&>(from));

  // get BaseField iterators
  BF_t::iterator bftoit, bftoend = bfto.end(), bffromit = bffrom.begin(), 
    bffromend = bffrom.end();
  // swap adjacent faces from "from" into "to"
  int loc[3];
  int ncx = from.get_Mesh().get_ncx();
  int ncy = from.get_Mesh().get_ncy();
  int ncz = from.get_Mesh().get_ncz();
  for (bftoit = bfto.begin(); bftoit != bftoend; ++bffromit, ++bftoit) {
    bffromit.GetCurrentLocation(loc);
    if (loc[0] != 0)
      (*bftoit)(0) = bffromit.offset(-1,0,0)(1);
    else
      (*bftoit)(0) = 0;
    if (loc[0] != ncx-1)
      (*bftoit)(1) = bffromit.offset(+1,0,0)(0);
    else
      (*bftoit)(1) = 0;
    if (loc[1] != 0)
      (*bftoit)(2) = bffromit.offset(0,-1,0)(3);
    else
      (*bftoit)(2) = 0;
    if (loc[1] != ncy-1)
      (*bftoit)(3) = bffromit.offset(0,+1,0)(2);
    else
      (*bftoit)(3) = 0;
    if (loc[2] != 0)
      (*bftoit)(4) = bffromit.offset(0,0,-1)(5);
    else
      (*bftoit)(4) = 0;
    if (loc[2] != ncz-1)
      (*bftoit)(5) = bffromit.offset(0,0,+1)(4);
    else
      (*bftoit)(5) = 0;
  }
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::sum(const typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& f)
{
  Vektor<T,6> vsum = ::sum(f);
  T ssum = 0;
  for (int face=0; face<6; ++face)
    ssum += vsum(face);
  return ssum;
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::min(const typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& f)
{
  Vektor<T,6> vmin = ::min(f);
  T smin = vmin(0);
  for (int face=1; face<6; ++face)
    if (vmin(face) < smin) smin = vmin(face);
  return smin;
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::max(const typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& f)
{
  Vektor<T,6> vmax = ::max(f);
  T smax = vmax(0);
  for (int face=1; face<6; ++face)
    if (vmax(face) > smax) smax = vmax(face);
  return smax;
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::sum(const typename PoomaMesh_XYZ<Mesh>::cctf<T>& f)
{
  return ::sum(f);
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::min(const typename PoomaMesh_XYZ<Mesh>::cctf<T>& f)
{
  return ::min(f);
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::max(const typename PoomaMesh_XYZ<Mesh>::cctf<T>& f)
{
  return ::max(f);
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::sum(const typename PoomaMesh_XYZ<Mesh>::bstf<T>& f)
{
  Vektor<T,6> vsum = ::sum(f);
  T ssum = 0;
  for (int face=0; face<6; ++face)
    ssum += vsum(face);
  return ssum;
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::min(const typename PoomaMesh_XYZ<Mesh>::bstf<T>& f)
{
  Vektor<T,6> vmin = ::min(f);
  T smin = vmin(0);
  for (int face=1; face<6; ++face)
    if (vmin(face) < smin) smin = vmin(face);
  return smin;
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::max(const typename PoomaMesh_XYZ<Mesh>::bstf<T>& f)
{
  Vektor<T,6> vmax = ::max(f);
  T smax = vmax(0);
  for (int face=1; face<6; ++face)
    if (vmax(face) > smax) smax = vmax(face);
  return smax;
}

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.t.cc
//---------------------------------------------------------------------------//
