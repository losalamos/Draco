#include "3T/testP13T/testFullP13T.hh"
#include "timestep/field_ts_advisor.t.cc"

typedef rtt_timestep::field_ts_advisor ftsa;

typedef Mesh_XYZ MT;
typedef MT::ccsf T1;

template
void ftsa::update_tstep<T1>(const rtt_timestep::ts_manager &, const T1 &,
			    const T1 &);
