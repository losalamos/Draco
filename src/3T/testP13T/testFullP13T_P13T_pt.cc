#include "mesh/Mesh_XYZ.hh"
#include "3T/testP13T/GmvDump.t.cc"
#include "3T/testP13T/utils.t.cc"

#include <functional>

typedef Mesh_XYZ MT;

namespace rtt_3T_testP13T
{
 template class GmvDump<MT>;
 template
 bool isContinuous<MT>(const MT::fcdsf &, const dsxx::SP<MT> &);

 typedef MT::fcdsf fcdsf;

 template
 fcdsf::value_type min(const fcdsf &,
		       std::pointer_to_unary_function<double,double> &);
 template
 fcdsf::value_type max(const fcdsf &,
		       std::pointer_to_unary_function<double,double> &);
 
 typedef MT::ccsf ccsf;

 template
 ccsf::value_type min(const ccsf &,
		      std::pointer_to_unary_function<double,double> &);
 template
 ccsf::value_type max(const ccsf &,
		      std::pointer_to_unary_function<double,double> &);

 template
 double sum(const ccsf &);
 
}
