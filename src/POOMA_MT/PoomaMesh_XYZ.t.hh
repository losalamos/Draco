//----------------------------------*-C++-*----------------------------------//
// PoomaMesh_XYZ.t.hh
// Julian C. Cummings
// Wed Jan 27 1999
//---------------------------------------------------------------------------//
// @> Definitions for PoomaMesh_XYZ
//---------------------------------------------------------------------------//


template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
			     const typename PoomaMesh_XYZ<Mesh>::cctf<T2>& from,
			     const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // update guard cells on "from" field
    typename PoomaMesh_XYZ<Mesh>::cctf<T2>& ncfrom =
	const_cast<typename PoomaMesh_XYZ<Mesh>::cctf<T2>&>(from);
    ncfrom.fillGuardCells();

    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
    typename PoomaMesh_XYZ<Mesh>::cctf<T2>::const_iterator fromit,
	fromend = from.end();
    // scatter from cells to each local face
    for (fromit = from.begin(); fromit != fromend; ++fromit)
	for (int f=0; f<6; ++f, ++toit)
	    PETE_apply(op, *toit, *fromit);

    // get BaseField reference
    typedef typename PoomaMesh_XYZ<Mesh>::cctf<T2>::BaseField_t BF_t;
    const BF_t& bffrom(dynamic_cast<const BF_t&>(from));

    // get BaseField iterator
    typename BF_t::iterator bffromit, bffromend = bffrom.end();

    // scatter from cells to each adjacent face
    int loc[3];
    int ncx = from.get_Mesh().get_ncx();
    int ncy = from.get_Mesh().get_ncy();
    int ncz = from.get_Mesh().get_ncz();
    toit = to.begin();
    for (bffromit = bffrom.begin(); bffromit != bffromend; ++bffromit) {
	bffromit.GetCurrentLocation(loc);
	if (loc[0] != 0)
	    PETE_apply( op, *toit, bffromit.offset(-1, 0, 0) );
	++toit;
	if (loc[0] != ncx-1)
	    PETE_apply( op, *toit, bffromit.offset(+1, 0, 0) );
	++toit;
	if (loc[1] != 0)
	    PETE_apply( op, *toit, bffromit.offset( 0,-1, 0) );
	++toit;
	if (loc[1] != ncy-1)
	    PETE_apply( op, *toit, bffromit.offset( 0,+1, 0) );
	++toit;
	if (loc[2] != 0)
	    PETE_apply( op, *toit, bffromit.offset( 0, 0,-1) );
	++toit;
	if (loc[2] != ncz-1)
	    PETE_apply( op, *toit, bffromit.offset( 0, 0,+1) );
	++toit;
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::cctf<T1>& to,
			     const typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>& from,
			     const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::cctf<T1>::iterator toit, toend = to.end();
    typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>::const_iterator fromit =
	from.begin();
    // scatter from each local face to cell center
    for (toit = to.begin(); toit != toend; ++toit)
	for (int f=0; f<6; ++f, ++fromit)
	    PETE_apply(op, *toit, *fromit);
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
			     const typename PoomaMesh_XYZ<Mesh>::vctf<T2>& from,
			     const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get BaseField references
    typedef typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>::BaseField_t BF1_t;
    typedef typename PoomaMesh_XYZ<Mesh>::vctf<T2>::BaseField_t BF2_t;
    BF1_t& bfto = dynamic_cast<BF1_t&>(to);
    const BF2_t& bffrom = dynamic_cast<const BF2_t&>(from);

    // get BaseField iterators, which iterate over cells
    typename BF1_t::iterator toit, toend = bfto.end();
    typename BF2_t::iterator fromit = bffrom.begin();

    // loop over cells and scatter from vertices to adjoining faces
    for (toit = bfto.begin(); toit != toend; ++toit, ++fromit) {
	// face 0
	PETE_apply( op, (*toit)(0), (*fromit)(0) );
	PETE_apply( op, (*toit)(0), (*fromit)(2) );
	PETE_apply( op, (*toit)(0), (*fromit)(4) );
	PETE_apply( op, (*toit)(0), (*fromit)(6) );
	// face 1
	PETE_apply( op, (*toit)(1), (*fromit)(1) );
	PETE_apply( op, (*toit)(1), (*fromit)(3) );
	PETE_apply( op, (*toit)(1), (*fromit)(5) );
	PETE_apply( op, (*toit)(1), (*fromit)(7) );
	// face 2
	PETE_apply( op, (*toit)(2), (*fromit)(0) );
	PETE_apply( op, (*toit)(2), (*fromit)(1) );
	PETE_apply( op, (*toit)(2), (*fromit)(4) );
	PETE_apply( op, (*toit)(2), (*fromit)(5) );
	// face 3
	PETE_apply( op, (*toit)(3), (*fromit)(2) );
	PETE_apply( op, (*toit)(3), (*fromit)(3) );
	PETE_apply( op, (*toit)(3), (*fromit)(6) );
	PETE_apply( op, (*toit)(3), (*fromit)(7) );
	// face 4
	PETE_apply( op, (*toit)(4), (*fromit)(0) );
	PETE_apply( op, (*toit)(4), (*fromit)(1) );
	PETE_apply( op, (*toit)(4), (*fromit)(2) );
	PETE_apply( op, (*toit)(4), (*fromit)(3) );
	// face 5
	PETE_apply( op, (*toit)(5), (*fromit)(4) );
	PETE_apply( op, (*toit)(5), (*fromit)(5) );
	PETE_apply( op, (*toit)(5), (*fromit)(6) );
	PETE_apply( op, (*toit)(5), (*fromit)(7) );
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::vctf<T1>& to,
			     const typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>& from,
			     const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get BaseField references
    typedef typename PoomaMesh_XYZ<Mesh>::vctf<T1>::BaseField_t BF1_t;
    typedef typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>::BaseField_t BF2_t;
    BF1_t& bfto = dynamic_cast<BF1_t&>(to);
    const BF2_t& bffrom = dynamic_cast<const BF2_t&>(from);

    // get BaseField iterators, which iterate over cells
    typename BF1_t::iterator toit, toend = bfto.end();
    typename BF2_t::iterator fromit = bffrom.begin();

    // loop over cells and scatter from faces to their vertices
    for (toit = bfto.begin(); toit != toend; ++toit, ++fromit) {
	// vertex 0
	PETE_apply( op, (*toit)(0), (*fromit)(0) );
	PETE_apply( op, (*toit)(0), (*fromit)(2) );
	PETE_apply( op, (*toit)(0), (*fromit)(4) );
	// vertex 1
	PETE_apply( op, (*toit)(1), (*fromit)(1) );
	PETE_apply( op, (*toit)(1), (*fromit)(2) );
	PETE_apply( op, (*toit)(1), (*fromit)(4) );
	// vertex 2
	PETE_apply( op, (*toit)(2), (*fromit)(0) );
	PETE_apply( op, (*toit)(2), (*fromit)(3) );
	PETE_apply( op, (*toit)(2), (*fromit)(4) );
	// vertex 3
	PETE_apply( op, (*toit)(3), (*fromit)(1) );
	PETE_apply( op, (*toit)(3), (*fromit)(3) );
	PETE_apply( op, (*toit)(3), (*fromit)(4) );
	// vertex 4
	PETE_apply( op, (*toit)(4), (*fromit)(0) );
	PETE_apply( op, (*toit)(4), (*fromit)(2) );
	PETE_apply( op, (*toit)(4), (*fromit)(5) );
	// vertex 5
	PETE_apply( op, (*toit)(5), (*fromit)(1) );
	PETE_apply( op, (*toit)(5), (*fromit)(2) );
	PETE_apply( op, (*toit)(5), (*fromit)(5) );
	// vertex 6
	PETE_apply( op, (*toit)(6), (*fromit)(0) );
	PETE_apply( op, (*toit)(6), (*fromit)(3) );
	PETE_apply( op, (*toit)(6), (*fromit)(5) );
	// vertex 7
	PETE_apply( op, (*toit)(7), (*fromit)(1) );
	PETE_apply( op, (*toit)(7), (*fromit)(3) );
	PETE_apply( op, (*toit)(7), (*fromit)(5) );
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::nctf<T1>& to,
			     const typename PoomaMesh_XYZ<Mesh>::vctf<T2>& from,
			     const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // update guard cells on "from" field
    typename PoomaMesh_XYZ<Mesh>::vctf<T2>& ncfrom =
	const_cast<typename PoomaMesh_XYZ<Mesh>::vctf<T2>&>(from);
    ncfrom.fillGuardCells();

    // get BaseField references
    typedef typename PoomaMesh_XYZ<Mesh>::nctf<T1>::BaseField_t BF1_t;
    typedef typename PoomaMesh_XYZ<Mesh>::vctf<T2>::BaseField_t BF2_t;
    BF1_t& bfto = dynamic_cast<BF1_t&>(to);
    const BF2_t& bffrom = dynamic_cast<const BF2_t&>(from);

    // get BaseField iterators
    typename BF1_t::iterator toit, toend = bfto.end();
    typename BF2_t::iterator fromit = bffrom.begin();
 
    FieldLoc<3> loc;
    int iloc[3];
    int ncx = from.get_Mesh().get_ncx();
    int ncy = from.get_Mesh().get_ncy();
    int ncz = from.get_Mesh().get_ncz();
    // loop over cells and find the vertices contributing to each node
    for (toit = bfto.begin(); toit != toend; ++toit) {
	// find cell location for this node and set iterator
	toit.GetCurrentLocation(loc);
	fromit.SetCurrentLocation(loc);
	// get location so we can avoid boundary contributions
	fromit.GetCurrentLocation(iloc);
	// now gather into each vertex
	if (iloc[2] != ncz) {
	    if (iloc[1] != ncy) {
		if (iloc[0] != ncx)
		    PETE_apply( op, *toit, (*fromit)(0) );
		if (iloc[0] != 0)
		    PETE_apply( op, *toit, fromit.offset(-1, 0, 0)(1) );
	    }
	    if (iloc[1] != 0) {
		if (iloc[0] != ncx)
		    PETE_apply( op, *toit, fromit.offset( 0,-1, 0)(2) );
		if (iloc[0] != 0)
		    PETE_apply( op, *toit, fromit.offset(-1,-1, 0)(3) );
	    }
	}
	if (iloc[2] != 0) {
	    if (iloc[1] != ncy) {
		if (iloc[0] != ncx)
		    PETE_apply( op, *toit, fromit.offset( 0, 0,-1)(4) );
		if (iloc[0] != 0)
		    PETE_apply( op, *toit, fromit.offset(-1, 0,-1)(5) );
	    }
	    if (iloc[1] != 0) {
		if (iloc[0] != ncx)
		    PETE_apply( op, *toit, fromit.offset( 0,-1,-1)(6) );
		if (iloc[0] != 0)
		    PETE_apply( op, *toit, fromit.offset(-1,-1,-1)(7) );
	    }
	}
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::scatter(typename PoomaMesh_XYZ<Mesh>::cctf<T1>& to,
			     const typename PoomaMesh_XYZ<Mesh>::vctf<T2>& from,
			     const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::cctf<T1>::iterator toit, toend = to.end();
    typename PoomaMesh_XYZ<Mesh>::vctf<T2>::const_iterator fromit = from.begin();
	
    // loop over cells and scatter from each vertex to cell center
    for (toit = to.begin(); toit != toend; ++toit) {
	for (int v=0; v<8; ++v) {
	    PETE_apply(op, *toit, *fromit);
	    ++fromit;
	}
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
			    const typename PoomaMesh_XYZ<Mesh>::cctf<T2>& from,
			    const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
    typename PoomaMesh_XYZ<Mesh>::cctf<T2>::const_iterator fromit,
	fromend = from.end();
    // gather from cell center to each local face
    for (fromit = from.begin(); fromit != fromend; ++fromit)
	for (int f=0; f<6; ++f, ++toit)
	    PETE_apply(op, *toit, *fromit);
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::bstf<T1>& to,
			    const typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>& from,
			    const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::bstf<T1>::iterator toit, toend = to.end();
    typename PoomaMesh_XYZ<Mesh>::fcdtf<T2>::const_iterator fromit =
	from.begin();
    // loop over boundary faces and gather from face centers to boundary faces
    for (toit = to.begin(); toit != toend; ++toit) {
	fromit = toit; // set fcdtf iterator to same location
	PETE_apply(op, *toit, *fromit);
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>& to,
			    const typename PoomaMesh_XYZ<Mesh>::bstf<T2>& from,
			    const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::fcdtf<T1>::iterator toit = to.begin();
    typename PoomaMesh_XYZ<Mesh>::bstf<T2>::const_iterator fromit,
	fromend = from.end();
    // loop over boundary faces and gather from boundary faces to face centers
    for (fromit = from.begin(); fromit != fromend; ++fromit) {
	toit = fromit; // set fcdtf iterator to same location
	PETE_apply(op, *toit, *fromit);
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::vctf<T1>& to,
			    const typename PoomaMesh_XYZ<Mesh>::nctf<T2>& from,
			    const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // update guard cells on "from" field
    typename PoomaMesh_XYZ<Mesh>::nctf<T2>& ncfrom =
	const_cast<typename PoomaMesh_XYZ<Mesh>::nctf<T2>&>(from);
    ncfrom.fillGuardCells();

    // get BaseField references
    typedef typename PoomaMesh_XYZ<Mesh>::vctf<T1>::BaseField_t BF1_t;
    typedef typename PoomaMesh_XYZ<Mesh>::nctf<T2>::BaseField_t BF2_t;
    BF1_t& bfto = dynamic_cast<BF1_t&>(to);
    const BF2_t& bffrom = dynamic_cast<const BF2_t&>(from);

    // get BaseField iterators
    typename BF1_t::iterator toit, toend = bfto.end();
    typename BF2_t::iterator fromit = bffrom.begin();

    // loop over cells and find the nodes contributing to each vertex
    FieldLoc<3> loc;
    for (toit = bfto.begin(); toit != toend; ++toit) {
	// find node location for vertex 0
	toit.GetCurrentLocation(loc);
	fromit.SetCurrentLocation(loc);
	// now gather into each vertex
	PETE_apply( op, (*toit)(0), *fromit );
	PETE_apply( op, (*toit)(1), fromit.offset(1,0,0) );
	PETE_apply( op, (*toit)(2), fromit.offset(0,1,0) );
	PETE_apply( op, (*toit)(3), fromit.offset(1,1,0) );
	PETE_apply( op, (*toit)(4), fromit.offset(0,0,1) );
	PETE_apply( op, (*toit)(5), fromit.offset(1,0,1) );
	PETE_apply( op, (*toit)(6), fromit.offset(0,1,1) );
	PETE_apply( op, (*toit)(7), fromit.offset(1,1,1) );
    }
}

template <class Mesh>
template <class T1, class T2, class Op>
void
PoomaMesh_XYZ<Mesh>::gather(typename PoomaMesh_XYZ<Mesh>::vctf<T1>& to,
			    const typename PoomaMesh_XYZ<Mesh>::cctf<T2>& from,
			    const Op& op)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // get field iterators
    typename PoomaMesh_XYZ<Mesh>::vctf<T1>::iterator toit = to.begin();
    typename PoomaMesh_XYZ<Mesh>::cctf<T2>::const_iterator fromit,
	fromend = from.end();
    // loop over cells and gather to each vertex in cell
    for (fromit = from.begin(); fromit != fromend; ++fromit) {
	for (int v=0; v<8; ++v) {
	    PETE_apply(op, *toit, *fromit);
	    ++toit;
	}
    }
}

template <class Mesh>
template <class T>
void
PoomaMesh_XYZ<Mesh>::swap_faces(typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& to,
				const typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& from)
{
    PAssert( to.get_Mesh() == from.get_Mesh() );
    // update guard cells on "from" field
    typename PoomaMesh_XYZ<Mesh>::fcdtf<T>& ncfrom =
	const_cast<typename PoomaMesh_XYZ<Mesh>::fcdtf<T>&>(from);
    ncfrom.fillGuardCells();

    // get BaseField references
    typedef typename PoomaMesh_XYZ<Mesh>::fcdtf<T>::BaseField_t BF_t;
    BF_t& bfto(dynamic_cast<BF_t&>(to));
    const BF_t& bffrom(dynamic_cast<const BF_t&>(from));

    // get BaseField iterators
    typename BF_t::iterator bftoit, bftoend = bfto.end(),
	bffromit = bffrom.begin(), bffromend = bffrom.end();
    // swap adjacent faces from "from" into "to"
    int loc[3];
    int ncx = from.get_Mesh().get_ncx();
    int ncy = from.get_Mesh().get_ncy();
    int ncz = from.get_Mesh().get_ncz();
    for (bftoit = bfto.begin(); bftoit != bftoend; ++bffromit, ++bftoit) {
	bffromit.GetCurrentLocation(loc);
	if (loc[0] != 0)
	    (*bftoit)(0) = bffromit.offset(-1, 0, 0)(1);
	else
	    (*bftoit)(0) = 0;
	if (loc[0] != ncx-1)
	    (*bftoit)(1) = bffromit.offset(+1, 0, 0)(0);
	else
	    (*bftoit)(1) = 0;
	if (loc[1] != 0)
	    (*bftoit)(2) = bffromit.offset( 0,-1, 0)(3);
	else
	    (*bftoit)(2) = 0;
	if (loc[1] != ncy-1)
	    (*bftoit)(3) = bffromit.offset( 0,+1, 0)(2);
	else
	    (*bftoit)(3) = 0;
	if (loc[2] != 0)
	    (*bftoit)(4) = bffromit.offset( 0, 0,-1)(5);
	else
	    (*bftoit)(4) = 0;
	if (loc[2] != ncz-1)
	    (*bftoit)(5) = bffromit.offset( 0, 0,+1)(4);
	else
	    (*bftoit)(5) = 0;
    }
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
PoomaMesh_XYZ<Mesh>::sum(const typename PoomaMesh_XYZ<Mesh>::nctf<T>& f)
{
    return ::sum(f);
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::sum(const typename PoomaMesh_XYZ<Mesh>::vctf<T>& f)
{
    Vektor<T,8> vsum = ::sum(f);
    T ssum = 0;
    for (int vert=0; vert<8; ++vert)
	ssum += vsum(vert);
    return ssum;
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
PoomaMesh_XYZ<Mesh>::min(const typename PoomaMesh_XYZ<Mesh>::nctf<T>& f)
{
    return ::min(f);
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::min(const typename PoomaMesh_XYZ<Mesh>::vctf<T>& f)
{
    Vektor<T,8> vmin = ::min(f);
    T smin = vmin(0);
    for (int vert=1; vert<8; ++vert)
	if (vmin(vert) < smin) smin = vmin(vert);
    return smin;
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
PoomaMesh_XYZ<Mesh>::max(const typename PoomaMesh_XYZ<Mesh>::bstf<T>& f)
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
PoomaMesh_XYZ<Mesh>::max(const typename PoomaMesh_XYZ<Mesh>::nctf<T>& f)
{
    return ::max(f);
}

template <class Mesh>
template <class T>
inline T
PoomaMesh_XYZ<Mesh>::max(const typename PoomaMesh_XYZ<Mesh>::vctf<T>& f)
{
    Vektor<T,8> vmax = ::max(f);
    T smax = vmax(0);
    for (int vert=1; vert<8; ++vert)
	if (vmax(vert) > smax) smax = vmax(vert);
    return smax;
}

//---------------------------------------------------------------------------//
//                              end of PoomaMesh_XYZ.t.hh
//---------------------------------------------------------------------------//
