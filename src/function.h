//
// Created by Hanxiao Wu on 7/15/23.
//

#ifndef MCMC_CPLUS_FUNCTION_H
#define MCMC_CPLUS_FUNCTION_H
#include "GlobalParameter.h"
#include <cmath>
#include <random>
#include <vector>
#include <string>
#include "MC.h"
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using namespace std;

double get_vp_brocher(const double &vs);
double get_rho_brocher(const double &v, const char &flag);
double get_rho_hacker(const double &vs);

double rand_uniform(const double &left,const double &right, const double &value);
AdditionalPert rand_uniform_add(const AdditionalPert &left, const AdditionalPert &right, const AdditionalPert &v);
double rand_gaussian(const double &mean, const double &sigma, const double &left, const double &right);
AdditionalPert rand_gaussian_add(const AdditionalPert &mean, const AdditionalPert &sigma, const AdditionalPert &left, const AdditionalPert &right);

double add_pertb(const vector<AdditionalPert> &add_pert, const double &Thk, const double &dep, double value);

int find_closest_value(const vector<double> &v, const double &x); // find the location of the value which is the closest to x

vector<double> joint_lyr(vector<double> v1, vector<double> v2);

double Gm_d(const vector<double> &pre, const vector<double> &obs, const vector<double> &err);
double reduceLbyhand(double mis);

void get_bound(const vector<double> &v, const vector<double> &radius, vector<double> &left, vector<double> &right);
void get_bound(const vector<AdditionalPert> &v, const vector<AdditionalPert> &radius, vector<AdditionalPert> &left, vector<AdditionalPert> &right);

bool check_monot(const vector<double> &v);

void Split(const string &s, vector<string> &v, const string& del);

Model average_m(const Model &m, const vector<vector<Group> > &groups);
#endif //MCMC_CPLUS_FUNCTION_H
