#include "vectorops.hh"

#include <cmath>

vector<double> pow(const vector<double> &lhs, double n)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = std::pow(lhs[i], n);

    return results;
}

vector<double> sqrt(const vector<double> &lhs)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = std::sqrt(lhs[i]);

    return results;
}

vector<double> log(const vector<double> &lhs)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = std::log(lhs[i]);

    return results;
}

vector<double> max(const vector<double> &lhs, const vector<double> &rhs)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = std::max(lhs[i], rhs[i]);

    return results;
}

vector<double> max(const vector<double> &lhs, double rhs)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = std::max(lhs[i], rhs);

    return results;
}

vector<double> max(double lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = std::max(lhs, rhs[i]);

    return results;
}

vector<double> operator/(const vector<double> &lhs,
				const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs[i]/rhs[i];

    return results;
}

vector<double> operator/(double lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs/rhs[i];

    return results;
}

vector<double> operator/(const vector<double> &lhs, double rhs)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = lhs[i]/rhs;

    return results;
}

vector<double> operator*(const vector<double> &lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs[i]*rhs[i];

    return results;
}

vector<double> operator*(double lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs*rhs[i];

    return results;
}

vector<double> operator*(const vector<double> &lhs, double rhs)
{
    return rhs * lhs;
}

vector<double> operator-(const vector<double> &lhs,
				const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs[i]-rhs[i];

    return results;
}

vector<double> operator-(double lhs, const vector<double> &rhs)
{
    vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs-rhs[i];

    return results;
}

vector<double> operator-(const vector<double> &lhs, double rhs)
{
    vector<double> results(lhs.size());

    for (int i=0; i<lhs.size(); i++)
	results[i] = lhs[i]-rhs;

    return results;
}

