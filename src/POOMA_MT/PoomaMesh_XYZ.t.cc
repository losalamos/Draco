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
  to.Uncompress();
  PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
  PoomaMesh_XYZ<Mesh>::cctf<T2>::const_iterator fromit, fromend = from.end();
  for (fromit = from.begin(); fromit != fromend; ++fromit)
    for (int f=0; f<6; ++f, ++toit)
      PETE_apply(op, *toit, *fromit);

  // scatter to each adjacent face
  toit = to.begin();
  for (fromit = from.begin(); fromit != fromend; ++fromit) {
    PETE_apply(op, *toit, fromit.offset(-1,0,0));  ++toit;
    PETE_apply(op, *toit, fromit.offset(+1,0,0));  ++toit;
    PETE_apply(op, *toit, fromit.offset(0,-1,0));  ++toit;
    PETE_apply(op, *toit, fromit.offset(0,+1,0));  ++toit;
    PETE_apply(op, *toit, fromit.offset(0,0,-1));  ++toit;
    PETE_apply(op, *toit, fromit.offset(0,0,+1));  ++toit;
  }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::cctf<T1>& to,
  const typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>& from, const Op& op)
{
  // iterate over fields and scatter from each local face
  to.Uncompress();
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
  to.Uncompress();
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
  to.Uncompress();
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
  to.Uncompress();
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
PoomaMesh_XYZ<Mesh>::swap(typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& f1,
                          typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& f2)
{
  // make sure we're uncompressed for later iteration
  f1.Uncompress();
  f2.Uncompress();

  // first save a copy of one field to make this easier
  fcdtf<T> ft(f1.get_Mesh());
  ft = f1;

  // get BaseField references
  typedef typename PoomaMesh_XYZ<Mesh>::fcdtf<T>::BaseField_t BF_t;
  BF_t& bf1(dynamic_cast<BF_t&>(f1));
  BF_t& bf2(dynamic_cast<BF_t&>(f2));
  BF_t& bft(dynamic_cast<BF_t&>(ft));
  // get BaseField iterators
  BF_t::iterator bf1it, bf1end = bf1.end(), bf2it = bf2.begin(),
    bf2end = bf2.end(), bftit = bft.begin();
  // swap adjacent faces from f2 into f1
  for (bf1it = bf1.begin(); bf1it != bf1end; ++bf1it, ++bf2it) {
    (*bf1it)(0) = bf2it.offset(-1,0,0)(1);
    (*bf1it)(1) = bf2it.offset(+1,0,0)(0);
    (*bf1it)(2) = bf2it.offset(0,-1,0)(3);
    (*bf1it)(3) = bf2it.offset(0,+1,0)(2);
    (*bf1it)(4) = bf2it.offset(0,0,-1)(5);
    (*bf1it)(5) = bf2it.offset(0,0,+1)(4);
  }
  // swap adjacent faces from f1 copy into f2
  for (bf2it = bf2.begin(); bf2it != bf2end; ++bf2it, ++bftit) {
    (*bf2it)(0) = bftit.offset(-1,0,0)(1);
    (*bf2it)(1) = bftit.offset(+1,0,0)(0);
    (*bf2it)(2) = bftit.offset(0,-1,0)(3);
    (*bf2it)(3) = bftit.offset(0,+1,0)(2);
    (*bf2it)(4) = bftit.offset(0,0,-1)(5);
    (*bf2it)(5) = bftit.offset(0,0,+1)(4);
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
//                              end of PoomaMesh_XYZ.t.cc
//---------------------------------------------------------------------------//
