//----------------------------------*-C++-*----------------------------------//
// dummy_package.cc
// John McGhee
// Thu Aug 27 07:48:41 1998
//---------------------------------------------------------------------------//
// @> A dummy package to exercize the time-step controller field advisors.
//---------------------------------------------------------------------------//

#include "timestep/test/dummy_package.hh"

#include "timestep/ts_manager.hh"

#include "timestep/field_ts_advisor.hh"

#include <vector>

#include "timestep/test/test_utils.hh"

using std::vector;

using namespace rtt_timestep;

vector<double> operator*(double lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs*rhs[i];

    return results;
}

dummy_package::dummy_package(ts_manager &tsm_)
    : tsm(tsm_)
{

// Set up a Electron Temperature advisor

    sp_te = new field_ts_advisor( "Electron Temperature",
				   ts_advisor::max,
				   field_ts_advisor::a_mean);
    tsm.add_advisor(sp_te);

// Set up a Ion-Temperature advisor

    sp_ti = new field_ts_advisor("Ion Temperature",
				  ts_advisor::max,
				  field_ts_advisor::rc_mean);
    tsm.add_advisor(sp_ti);

// Set up a Radiation-Intensity advisor

    sp_ri = new field_ts_advisor("Radiation Intensity",
				  ts_advisor::max,
				  field_ts_advisor::rcq_mean);
    tsm.add_advisor(sp_ri);
}

dummy_package::~dummy_package()
{
    tsm.remove_advisor(sp_te);
    tsm.remove_advisor(sp_ti);
    tsm.remove_advisor(sp_ri);
}

void dummy_package::advance_state()
{

// Create a set of dummy arrays to serve as control fields for
// use in exercizing the various advisors.

    const double a1[] = {1., 10., 11., 3., 2., 5., 5., 6.7};
    int sizea = sizeof(a1)/sizeof(a1[0]);
    vector<double> te_old(a1,a1+sizea);
    vector<double> te_new = 1.09*te_old;
    vector<double> ti_old=0.97*te_old;
    vector<double> ti_new=1.05*te_old;
    vector<double> ri_old=1.10*te_old;
    vector<double> ri_new=1.15*te_old;

// Set a floor for the electron temperature controller, to
// execcize this method. Just accelpt the default floor on the
// other controllers.

    sp_te -> set_floor(te_new,0.001);

// Get a new time-step from each of the advisors that
// belong to this package.

    sp_te -> update_tstep(tsm, te_old, te_new);
    sp_ti -> update_tstep(tsm, ti_old, ti_new);
    sp_ri -> update_tstep(tsm, ri_old, ri_new);
}

bool dummy_package::tests_passed() const
{
    const int nd = 5;
    return compare_reals(1.371742e+00, sp_te->get_dt_rec(tsm),nd)
        && compare_reals(1.496914e+00, sp_ti->get_dt_rec(tsm),nd)
        && compare_reals(2.716049e+00, sp_ri->get_dt_rec(tsm),nd);
}

//---------------------------------------------------------------------------//
//                              end of dummy_package.cc
//---------------------------------------------------------------------------//
