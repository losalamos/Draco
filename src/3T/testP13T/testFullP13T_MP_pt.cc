#include "3T/testP13T/testFullP13T.hh"
#include "matprops/InterpedMaterialProps.t.cc"

typedef Mesh_XYZ MT;
typedef XTM::InterpedMaterialProps MP;

typedef MT::ccsf ccsf;
typedef MT::fcdsf fcdsf;
typedef MT::ccif ccif;
typedef MT::fcdif fcdif;

#define T1 ccsf
#define T2 ccif

template
MP::MaterialStateField<T1>
MP::getMaterialState<T1, T2>(const T1 &, const T1 &, const T1 &,
			     const T2 &) const;

#undef T1
#undef T2
#define T1 fcdsf
#define T2 fcdif

template
MP::MaterialStateField<T1>
MP::getMaterialState<T1, T2>(const T1 &, const T1 &, const T1 &,
			     const T2 &) const;

typedef XTM::BilinearInterpTable BIT;

#undef T1
#undef T2
#define T1 ccsf
#define T2 MP::MaterialStateField<ccsf>::MultByDensity

template
void
MP::interpolate<T1, T2>(const MP::MaterialStateField<T1> &, int,
			const MP::GroupedTable &(MP::MaterialTables::*)() const,
			T2, T1 &) const;

template
void
MP::interpolate<T1, T2>(const MP::MaterialStateField<T1> &,
			const BIT &(MP::MaterialTables::*)() const,
			T2, T1 &) const;

#undef T1
#undef T2
#define T1 fcdsf
#define T2 MP::MaterialStateField<fcdsf>::MultByDensity

template
void
MP::interpolate<T1, T2>(const MP::MaterialStateField<T1> &, int,
			const MP::GroupedTable &(MP::MaterialTables::*)() const,
			T2, T1 &) const;
template
void
MP::interpolate<T1, T2>(const MP::MaterialStateField<T1> &,
			const BIT &(MP::MaterialTables::*)() const,
			T2, T1 &) const;

#undef T1
#undef T2
