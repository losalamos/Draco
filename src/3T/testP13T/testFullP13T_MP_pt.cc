#include "3T/testP13T/testFullP13T.hh"

#include "matprops/InterpedMaterialProps.t.cc"

typedef Mesh_XYZ MT;
typedef rtt_matprops::InterpedMaterialProps MP;

typedef MT::ccsf TS1;
typedef MT::ccif TI1;

template
class MP::MaterialStateField<TS1,TI1>;

typedef MT::fcdsf TS2;
typedef MT::fcdif TI2;

template
class MP::MaterialStateField<TS2,TI2>;

#include <vector>
typedef std::vector<double>  TS3;
typedef std::vector<int   >  TI3;

template
class MP::MaterialStateField<TS3,TI3>;

template
MP::MaterialStateField<TS3, TI3> MP::getMaterialState(const TS3 &, 
						      const TS3 &, 
						      const TS3 &, 
						      const TI3 &) const;

