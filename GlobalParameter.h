//
// Created by Hanxiao Wu on 7/14/23.
//

#ifndef MCMC_CPLUS_GLOBALPARAMETER_H
#define MCMC_CPLUS_GLOBALPARAMETER_H
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;
const double threshold_E = 0.9; // todo: to be changed
const double threshold_mis = 0.5;
const double threshold_sw_one = 1.;
const double earth_radius = 6371.;
const double pi = atan(1)*4;
const int N_L = 2; //#likelihood
//const int NUMPERT = 3; // current max #sub_perturb parameters: vp(vp/vs),vs,rho
const int Nbad_threshold =2000; // if have searched 2000+ models but still could find a model satisfying the prior constraints, exit

const double VPVS = 1.75; //default vpvs
const double ice_vs = 1.94; // Group_tp = -1
const double ice_vp = 3.87; // vp_tp = -1
const double ice_rho = 0.92; // rho_tp = -1

const double water_vs = 0.; // Group_tp = 0
const double water_vp = 1.45; // vp_tp = 0
const double water_vpvs = -1.;  //vp_tp = 0
const double water_rho = 1.02; // rho_tp = 0

const double ice_Qs = 80.; //Q_tp = -1
const double ice_Qp = 160.; //Q_tp = -1
const double water_Qs = 0.; //Q_tp = 0
const double water_Qp = 57822.; //Q_tp = 0
const double sedi_Qs = 80.; // Q_tp = 1
const double sedi_Qp = 160.; // Q_tp =1
const double cru_Qs = 599.99; // Q_tp = 2
const double cru_Qp = 1368.02; //Q_tp = 2
const double man_Qs = 417.59; // Q_tp = 3
const double man_Qp = 1008.71; // Q_tp = 3

const double man_vpvs = 1.789;

const int factor_E = 20; // used to calculate 'likelihood' of Energy
const double factor_M = 0.5; // used to calculate 'likelihood' of Misfit
const double MAX_VS = 4.9; //max Vs

const int SizeOfInt = sizeof(int);
const int SizeOfDouble = sizeof(double);
const double DEFAULTd = -1234.;
const int DEFAULTi = -12345;
const int N_property = 7;

const int Nrf = 2000; //max number of rf points

template <typename T>
vector<T> operator+(const vector<T> &a, const vector<T> &b){
    if(a.size() != b.size()) cout<<"length is not equal"<<endl;
    vector<T> re;re.clear();
    for (int i=0;i<a.size();i++) re.push_back(a[i]+b[i]);
    return re;
}
template <typename T>
vector<T> operator-(const vector<T> &a, const vector<T> &b){
    if(a.size() != b.size()) cout<<"length is not equal"<<endl;
    vector<T> re;re.clear();
    for (int i=0;i<a.size();i++) re.push_back(a[i]-b[i]);
    return re;
}
template <typename T>
vector<T> operator/(const vector<T> &a,const int &n){
    vector<T> res;res.clear();
    for(int i=0;i<a.size();i++) res.push_back(a[i]/n);
    //cout<<"vector divide:"<<res.size()<<endl;
    return res;
}
#endif //MCMC_CPLUS_GLOBALPARAMETER_H
