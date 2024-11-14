//
// Created by Hanxiao Wu on 7/11/23.
//

#ifndef MCMC_CPLUS_MC_H
#define MCMC_CPLUS_MC_H
#include <iostream>
#include <fstream>
//#include <cmath>
#include "GlobalParameter.h"
//#include "function.h"
//#include <random>
//#include <vector>
#include "SAC_IO/sac_stream.hpp"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using namespace std;
class swData;
class RFData;
class HkData;
class Lyr_model;
class AdditionalPert;
class Property;
class Group;
class subSample;
class Sample;
class Model;

class swData {
    vector<double>  T, data, err, ref_data; // period, data, data uncertainty
    char sw_tp; // type of sw, Rayleigh or Love
    char data_tp; // type of data
    int n;
public:
    swData(const char &sw_tp, const char &tp);
    ~swData();

    [[nodiscard]] char get_tp()  const;
    [[nodiscard]] vector<double> get_data() const;
    [[nodiscard]] vector<double> get_err() const;

    int read_from_file(const string &fname);
    int print() const;
    void Clear();
    friend class SWData;
    void write_to_txt(const string &nm) const;
};
class SWData{
    vector<swData> sw;
    int n; //#type of sw data
public:
    SWData();
    ~SWData();
    void add_sw(const swData &sw);
    void set_n(const int &n);
    void Clear();

    int Disp(const Lyr_model &lyrm);

    [[nodiscard]] int get_n() const;
    [[nodiscard]] vector<double> get_data(const int &idx);
    [[nodiscard]] vector<double> get_err(const int &idx);

    void print_swdata() const;
    void write_to_txt(const string &nm) const;
    vector<char> get_disp_tp() const;
    vector<char> get_sw_tp() const;
};

class RFData {
    vector<double> t, rf, err; // RF waveform
    double a, rayp, dist, baz,snr,fitness; //gaussian parameter, ray parameter
    int n;
    vector<double> E, Eps, Eppps, Epsps; //
    vector<double> TPs, TPpPs, TPsPs; //each stores the arrival time from different disco
    vector<int> disco; // indicate the location of discontinuity
public:
    RFData();
    RFData(const double &a, const double &rayp);
    RFData(const vector<double> &time, const vector<double> &amplitude, const vector<double> &err, const double &a, const double &rayp, const vector<int> &disco );
    ~RFData();
    int set_disco(const vector<int> &disco);
    int read_from_file(const string &fname);
    int read_from_sac(const string &fname);
    int RF(const Lyr_model &lyrm);

    [[nodiscard]] double get_a() const;
    [[nodiscard]] double get_rayp() const;
    [[nodiscard]] vector<double> get_t() const;
    [[nodiscard]] vector<double> get_rf() const ;
    [[nodiscard]] vector<double> get_err() const ;
    [[nodiscard]] vector<double> get_Eps() const ;
    [[nodiscard]] vector<double> get_Eppps() const ;
    [[nodiscard]] vector<double> get_Epsps() const ;
    [[nodiscard]] vector<double> get_Tps() const;
    [[nodiscard]] vector<double> get_Tppps() const;
    [[nodiscard]] vector<double> get_Tpsps() const;
    [[nodiscard]] vector<double> get_E() const ;
    int cal_arrival(const Lyr_model &lyrm);
    int cal_Energy(const double &wPs, const double &wPsPs, const double &wPpPs); // the arrival should store the arrival of different phases converted from the same discontinuity
    void print() const;
    void write_to_txt(const string &nm) const;
    void write_to_sac(const string &nm) const;
};

class HkData {
    double a, rayp, dist,baz; //gaussian parameter, ray parameter
    vector<double> E, Eps, Eppps, Epsps; //
    vector<double> TPs, TPpPs, TPsPs; //each stores the arrival time from different disco
    vector<int> disco; // indicate the location of discontinuity
public:
    HkData();
    ~HkData();

    void set_idsco(const vector<int> &idsco);
    void set(const double &a, const double &rayp, const double &dist, const double &baz);
    void set_E(vector<double> E, vector<double> Eps, vector<double> Eppps, vector<double> Epsps);
    void set_T(vector<double> Tps, vector<double> Tppps, vector<double> Tpsps);
    [[nodiscard]] vector<double> get_Eps() const ;
    [[nodiscard]] vector<double> get_Eppps() const ;
    [[nodiscard]] vector<double> get_Epsps() const ;
    [[nodiscard]] vector<double> get_E() const ;
    [[nodiscard]] double get_a() const;
    void write_to_txt(ofstream &f) const;
};


class Lyr_model {
    int nlyr; // #layers
    double Thk,dep_top; // need to be updated
    vector<double> vs;
    vector<double> vp;
    vector<double> rho;
    vector<double> vpvs;
    vector<double> Qs;
    vector<double> Qp;
    vector<double> thk;
    vector<double> T;
    vector<double> P;
    vector<double> depth;
//    vector<double> porosity;

public:
    Lyr_model();
    Lyr_model(double Thickness, double DepthOfTop,int num_lyr);
    ~Lyr_model();

    int set_Thk(double Thk);
    int set_topDep(double topdep);
    void set_nlyr(int n);

    vector<double> gradient(const vector<double> &v, const vector<AdditionalPert> &add);
    vector<double> layered(const vector<double> &v, const vector<AdditionalPert> &add);
    vector<double> bulk(const vector<double> &v, const vector<AdditionalPert> &add);
    vector<double> bspline(const vector<double> &v, const vector<AdditionalPert> &add);
    int update_vs(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_vs);
    int update_vp(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_vpvs);
    int update_rho(const int &flag, const vector<double> &v,const vector<AdditionalPert> &add_rho);
    int update_Qs(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_Qs); // set q according to Q_tp;
    int update_Qp(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_Qp); // set q according to Q_tp;
    int update_T(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_T);
    int update_P(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_P);

    [[nodiscard]] int get_nlyr() const ;
    [[nodiscard]] double get_Thk() const ;
    [[nodiscard]] double get_dep_top() const;
    [[nodiscard]] vector<double> get_vs() const;
    [[nodiscard]] vector<double> get_vpvs() const;
    [[nodiscard]] vector<double> get_vp() const;
    [[nodiscard]] vector<double> get_rho() const;
    [[nodiscard]] vector<double> get_thk() const;
    [[nodiscard]] vector<double> get_depth() const;
    [[nodiscard]] vector<double> get_Qs() const;
    [[nodiscard]] vector<double> get_Qp() const;
    [[nodiscard]] vector<double> get_T() const;
    [[nodiscard]] vector<double> get_P() const ;
    //vector<double> get_porosity() const;

    int set_vs(const vector<double> &vs);
    int set_thkdep(const vector<double> &thk,const vector<double> &dep);
    int set_vp(const vector<double> &vp);
    int set_vpvs(const vector<double> &vpvs);
    int set_rho(const vector<double> &rho);
    int set_Q(const vector<double> &Qs, const vector<double> &Qp);

    void print_lyr() const;
    void write_to_txt(ofstream &f) const;
//    int update_info();

};

class AdditionalPert{
    double top,bot,value;
    //int tp;
public:
    AdditionalPert();
    AdditionalPert(double top, double bot, double value);
    AdditionalPert operator+(const AdditionalPert& a) const;
    AdditionalPert operator- (const AdditionalPert& a) const;
    AdditionalPert operator/ (int n);
    int set(double top, double bot, double value);
    void set_top(double top);
    void set_bot(double bot);
    void set_value(double value);
    ~AdditionalPert();
    [[nodiscard]] double get_top() const ;
    [[nodiscard]] double get_bot() const ;
    [[nodiscard]] double get_value() const ;
    friend class Group;
};

class Property{
    vector<double> para; // parameters
    int n_para; // #parameters
    int tp; //interp type
    double thk, topdep;
    vector<AdditionalPert> add_pert;
    int n_pert;
public:
    ~Property();
    Property();
    Property(const int &tp, const vector<double> &para, const int &n_para, const vector<AdditionalPert> &add_pert, const int &n_pert, const double &thk, const double &dep);
    int set(int &index, const vector<string> &v);
    friend class Group;
    Property operator+(const Property &obj){
        Property res;
        if(n_pert != obj.n_pert or tp != obj.tp or n_para!= obj.n_para) exit(0);
        res.para = para + obj.para;
        res.thk = thk + obj.thk;
        res.add_pert = add_pert + obj.add_pert;
        res.n_para = n_para;
        res.tp = tp;
        res.n_pert = n_pert;
        res.topdep = topdep + obj.topdep;
        return res;
    }
    Property operator/(const int &n){
        Property res;
        res.para = para/n;
        res.thk = thk/n;
        //cout<<"Property::/"<<res.thk<<endl;
        //res.add_pert = add_pert/n;
        res.add_pert.clear();
        for(int i=0;i<n_pert;i++) res.add_pert.push_back(add_pert[i]/n);
        res.n_para = n_para;
        res.n_pert = n_pert;
        res.tp = tp;
        res.topdep = topdep/n;
        return res;
    }
	vector<double> get_para() const;
    void print() const;
    void write_to_bin(std::ofstream &f) const;
    void read_from_bin(std::ifstream &f);
};

class Group{
    Property vs,vpvs,rho,Qs,Qp,T,P;
    double thk, dep, P0; // dep is the depth of the top of this group
    Lyr_model lyrm;
    static int count;
    int idx;
    int nlyr;
public:
    ~Group();
    Group();
    explicit Group(const int &n_group);
    Group(const int &n_group, const int &nlyr,  const Property &vs, const Property &vpvs, const Property &rho, const Property &Qs, const Property &Qp, const Property &T, const Property &P);
    Group operator+(const Group &obj){
        if(nlyr!=obj.nlyr or idx != obj.idx) exit(0);
        Group res(obj);
        res.vs = vs+obj.vs;
        res.vpvs = vpvs+obj.vpvs;
        res.rho = rho + obj.rho;
        res.Qp = Qp + obj.Qp;
        res.Qs = Qs + obj.Qs;
        res.P = P + obj.P;
        res.T = T + obj.T;
        res.thk = thk + obj.thk;
        res.dep = dep + obj.dep;
        res.P0 = P0 + obj.P0;
        return res;
    }
    Group operator/(const int &n){
        Group res;
        res.vs = vs/n;
        //cout<<"Group::operator /"<<endl;
        res.vpvs = vpvs/n;
        res.rho = rho/n;
        res.Qp = Qp/n;
        res.Qs = Qs/n;
        res.P = P/n;
        res.T = T/n;
        res.thk = thk/n;
        res.dep = dep/n;
        res.P0 = P0/n;
        res.nlyr = nlyr;
        res.idx = idx;
        return res;
    }
    void set_nlyr(const int &n);
    int readin(const vector<string> &v);
    void set_dep(); //set thk of group = thk of vs
    void set_dep(const double &dep);
    int update_from_sample(const Sample &s);
    int update_from_sample(const Sample &s,const double thickness);
    void update_lyrmodel();

    [[nodiscard]] vector<double> get_lyr_vs() const;
    [[nodiscard]] vector<double> get_lyr_vp() const;
    [[nodiscard]] vector<double> get_lyr_rho() const;
    [[nodiscard]] vector<double> get_lyr_dep() const;
    [[nodiscard]] vector<double> get_lyr_Qs() const;
    [[nodiscard]] vector<double> get_lyr_Qp() const;
    [[nodiscard]] vector<double> get_lyr_thk() const;
    [[nodiscard]] vector<double> get_lyr_vpvs() const;

    [[nodiscard]] double get_thk() const;
    void print_info() const;
    [[nodiscard]] double get_dep() const;
    [[nodiscard]] int get_nlyr() const;
    void write_lyr() const;

    void write_to_file(ofstream &f) const;
	void write_v_to_file(ofstream &f) const;
    void write_to_bin(ofstream &f) const;
    void read_from_bin(ifstream &f);
    void print() const;
};

class subSample{
    int n_para,n_pert; // ???type of property
    vector<double> para, left, right, step;
    vector<AdditionalPert> add_para, add_left, add_right, add_step;
public:
    subSample();
    ~subSample();
    double Initial(int &index, const vector<string> &v);
//    int readin(int &index, const vector<string> &v);

    void set(const int &para_index, const double &radius,const  double &step);
    void set(const char &f, const int &para_index, const double &radius, const double &step);
    void uniform_rand();
    void gaussian_rand();

    [[nodiscard]] vector<double> get_para() const;
    [[nodiscard]] vector<AdditionalPert> get_addpara() const;
    void print_sample() const;
    void print() const;
    friend class Sample;
};

class Sample {
    subSample vs,vpvs,rho,Qs,Qp,T,P;
    double thk,left_thk, right_thk, step_thk;
    public:
    Sample();
    ~Sample();
    Sample(subSample vs, subSample vpvs, subSample rho, subSample Qs, subSample Qp, subSample T, subSample P, double thk, double left_thk, double right_thk, double step_thk);

    void set(const subSample &vs, const subSample &vpvs, const subSample &rho, const subSample &Qs, const subSample &Qp, const subSample &T, const subSample &P, const double &thk, const double &left_thk, const double &right_thk, const double &step_thk);
    int Initial(const vector<string> &v);
    int readin(const vector<string> &v);
    int rand_start_sample(); // randomly initialize the value within model space
    int gen_new_sample(); // randomly generate new sample based on old one
    void set_thk(const double &thk);

    [[nodiscard]] vector<double> get_vs() const;
    [[nodiscard]] vector<AdditionalPert> get_addvs() const;
    [[nodiscard]] vector<double> get_vpvs() const;
    [[nodiscard]] vector<AdditionalPert> get_addvpvs() const;
    [[nodiscard]] vector<double> get_rho() const;
    [[nodiscard]] vector<AdditionalPert> get_addrho() const;
    [[nodiscard]] vector<double> get_P() const;
    [[nodiscard]] vector<AdditionalPert> get_addP() const;
    [[nodiscard]] vector<double> get_Qs() const;
    [[nodiscard]] vector<AdditionalPert> get_addQs() const;
    [[nodiscard]] vector<double> get_Qp() const;
    [[nodiscard]] vector<AdditionalPert> get_addQp() const;
    [[nodiscard]] vector<double> get_T() const;
    [[nodiscard]] vector<AdditionalPert> get_addT() const;
    [[nodiscard]] double get_thk() const ;
    void print() const;
    void print_sample() const;
    void write_to_file(ofstream &f) const;
};

class Model {
    vector<Group> groups; // all the groups the model have
    Lyr_model lyrm; // a layered model
    double Thk;

    int f_monot;
    vector<int> idx_monot;
    int n_group;//#groups
    int f_acpt; //  a flag: 1=accepted;0=rejected
    int f_E; // a flag: 1=okay;others: have a wrong phase energy
    SWData pre_sw;
    vector<RFData> pre_rf; // predicted RF waveform
    vector<vector<HkData> > pre_hk; // predicted Hk-related data

    SWData obs_sw; // observed data
    vector<RFData> obs_rf;
    vector<vector<RFData> > obs_hk;

    vector<double> sw_mis_one, sw_S_one; //ph_mis,gr_mis,hv_mis,la_mis; // individual misfit for individual SW data
    vector<double> rf_mis_one, rf_S_one;
    double sw_mis, rf_mis, Joint_mis, sw_L, rf_L; // misfit for SW, RF, and joint misfit for all data
    double w_sw, w_rf; // weighting
    vector<double> wPs, wPsPs, wPpPs;
    vector<vector<double> > E_Ps, E_PpPs, E_PsPs; // energy for individual phase
    vector<vector<double> > E, wE;
    double totalE,refE;
    double L_mis, L_E; // likelihood of misfit & energy, respectively
    int a_indx; // the index of this model
    static int  acount ; // number of model, accepted model
public:
    Model();
    explicit Model(int n_groups);
    Model(const vector<Group> &groups, const SWData &obs_sw, const vector<RFData>& obs_rf, const vector<vector<RFData> >& obs_hk);
    ~Model();
//    void reset_count(int num);
//    void reset_indx(int idx);
    void set_refE(const double &refE);
    void set_weight(const double &w_rf,const vector<double> &wPs,const vector<double> &wPsPs,const vector<double> &wPpPs,const vector<vector<double> > &wE);
    void set_monot(const vector<int> &monot);
    void set_group(const vector<Group> &group);
//    todo: test the read in
    int readin_groups(const string &fname); // ??? shouldn't be here
    int readin_swdata(const vector<string> &fname,const vector<char> &tp, const vector<char> &sw_tp); // shouldn't be here
    int readin_rfdata(const vector<string> &fname, const vector<double> &a, const vector<double> &rayp);
    int readin_hkdata(const vector<vector<string> > &fname, const vector<vector<int>> &idisco); //????

    int update_groups(const vector<Sample> &s);
    int update_lyrmod(); // update model part
    int update_SWdata(); // update SW data part
    int update_RFdata();
    int update_Hkdata(int flag);
    int update_data(int flag=0);

    bool check_prior(); //?????, need to be done check sample if follows the prior

    int cal_SWmisfit(); // calculate the SW misfit & Likelihood
    int cal_RFmisfit(); // calculate the RF misfit & Likelihood
    int cal_Jmisfit(); // calculate the joint misfit & Likelihood
    int cal_Energy(int ndisco); // calculate the energy & Likelihood
    int cal_mis();
    //double cal_Likelihood(); // calculate the likelihood, return Likelihood

    [[nodiscard]] int get_n_sw() const;
    [[nodiscard]] int get_aidx() const;
    int get_idx() const;
    [[nodiscard]] int get_ngroup() const;
    [[nodiscard]] double get_Lmis() const;
    [[nodiscard]] double get_LE() const;
    [[nodiscard]] double get_Thk() const;
    [[nodiscard]] vector<double> get_sw_one() const;
    [[nodiscard]] double get_Dispmis() const;
    [[nodiscard]] vector<double> get_rf_one() const;
    [[nodiscard]] double get_RFmis() const;
    [[nodiscard]] double get_Jointmis() const;
    [[nodiscard]] double get_Energy() const;
    [[nodiscard]] vector<vector<double>> get_Eps() const;
    [[nodiscard]] vector<vector<double>> get_Eppps() const;
    [[nodiscard]] vector<vector<double>> get_Epsps() const;
    [[nodiscard]] vector<double> get_Lyr_vs() const;
    [[nodiscard]] Lyr_model get_lyrm() const;
    [[nodiscard]] int get_flag() const;
    [[nodiscard]] int get_flagE() const;
    [[nodiscard]] vector<Group> get_groups() const;
    double get_refE() const;

    int update_accept_flag(const int &flag);

    void print_groups() const;
    void print_swdata() const;
    void print_rfdata() const;
    void write_lyr_to_file(const string &dir,const string &fnm) const;
    void write_data_to_file(const string &nm, const int flag=0) const;
    void write_groups_to_file(const string &dir,const string &fnm) const;
    void write_to_binary(ofstream &f) const;
    void write_info(const string &dir, const string &fnm) const;
};

class MC {
    vector<Sample> sample0, sample1;
    vector<double> L0,L1;
    vector<double> p;
    int N; // #models
    //vector<Model> models;
    Model m0,m1;
    Model refm;

    double w_rf,w_sw;
    static int mcount;
    double min_disp,min_rf,min_joint;
    int idx_minDisp,idx_minRF,idx_minJoint;
    double max_E;
    int idx_maxE;

    vector<int> a_idx;
public:
    explicit MC(int N);
    MC(vector<Sample> s, Model m, int n);
    ~MC();
    void set_weight(const double &w_rf);
    int Initial_m0(const int &n_group,const string &fgroupnm, const vector<string> &fdispnm, const vector<char> &disp_tp, const vector<char> &sw_tp, const vector<string> &frfnm, const vector<double> &a, const vector<double> &rayp, const vector<vector<string> > &fhknm, const vector<vector<int>> &idisco, const double &w_rf,const vector<double> &wPs, const vector<double> &wPsPs,const vector<double> &wPpPs,const vector<vector<double> > &wE, const double &refE, const vector<int> &idx_monot);
    int Initial_s0(const string &fnm, const int &n_group);
    int read_s0(const string &fnm);

    //int print();
    int rand_starting_point();
    //int set_sample0();//{sample0.initial();} // initial once, maybe moved to constructor
    //int update_model();
    int gen_good_sample();
    int update_m1();
    void set_L0();
    void set_L1();
    void cal_p(); // set p
    int accept_model(ofstream &f, vector<Model> &m);//{if(flag){sample0=sample1;L0=L1;}}

    [[nodiscard]] double get_minDisp() const;
    [[nodiscard]] int get_idxminDisp() const;
    [[nodiscard]] double get_minRF() const;
    [[nodiscard]] int get_idxminRF() const;
    [[nodiscard]] double get_minJoint() const;
    [[nodiscard]] int get_idxminJoint() const;
    [[nodiscard]] double get_maxE() const;
    [[nodiscard]] int get_idxmaxE() const;
    //[[nodiscard]] vector<Model> get_acc_models() const;
    //[[nodiscard]] vector<Model> get_models() const;
    [[nodiscard]] vector<int> get_adix() const;
    void write_samples_to_file(int idx_thread, int idx_search, int idx_model, ofstream &f) const;
    void print_s0() const;
};


#endif //MCMC_CPLUS_MC_H
