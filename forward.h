//
// Created by Hanxiao Wu on 8/6/23.
//

#ifndef MCMC_CPLUS_FORWARD_H
#define MCMC_CPLUS_FORWARD_H
#include "MC.h"
#include<complex>
using namespace std::complex_literals;

vector<double> surface(Lyr_model lyrm, vector<double> periods, int nmodes );

void mchdep(); // todo to be added;

void gcmdln(bool &verby, double &ccmin, double &ccmax); // todo: to be added



void disprs(int ilvry, double dt, int npts, int iret, bool verby, int nfval, double fval, double ccmin, double ccmax);

// examine the velocity model. order the vel model in terms of increasing vel.
// this is done twice, once for S only and then for S&P combined.
// the objective is to determine regions where a denser search in phase vel. should be made to avoid missing modes
void mdsrch(double cmin,double cmax, Lyr_model lyrm, bool allfluid);

void bsort(vector<double> &vth,int flag);

void gtsolh(double vp, double vs, double cmin);

double dltar(int mvno, double omega, int kk, Lyr_model lyrm);
double dltar1(int mvno, double omega, Lyr_model lyrm);
double dltar4(int mvno, double omega, Lyr_model lyrm);

void gtegsh(int m, int wvno, double omega, double rsh, bool lshimag);
void gtesh(vector<complex<double> > esh, vector<complex<double> > einvsh, complex<double> rsh, int wvno, double mu, bool lshimag);
#endif //MCMC_CPLUS_FORWARD_H
