//
// Created by Hanxiao Wu on 7/11/23.
//
#include "function.h"
#include "MC.h"
//#include "function.h"
//#include "forward.h"
#include "GlobalParameter.h"
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <cmath>
#include <vector>
#include <string>
using namespace std;
/*
template<class T1, class T2>
class VecSum {
    int type_;
    const T1& op1_;
    const T2& op2_;
public:
    VecSum(int type, const T1& op1, const T2& op2): type_(type), op1_(op1), op2_(op2) {}

    int operator[](const int i) const {
        switch(type_) {
            case 0: return op1_[i] + op2_[i];
            case 1: return op1_[i] - op2_[i];
            case 2: return op1_[i] * op2_[i];
            case 3: return op1_[i] / op2_[i];
            default: throw "bad type";
        }
    }
};
template<class T1, class T2>
VecSum<T1, T2> operator+(const T1& t1, const T2& t2) {
    return VecSum<T1, T2>(0, t1, t2);
}

template<class T1, class T2>
VecSum<T1, T2> operator*(const T1& t1, const T2& t2) {
    return VecSum<T1, T2>(2, t1, t2);
}
template<class T1, class T2>
VecSum<T1, T2> operator-(const T1& t1, const T2& t2) {
    return VecSum<T1, T2>(1, t1, t2);
}
*/
extern"C"
{
void fast_surf_(int *n_layer0,int *kind0,float *a_ref0,float *b_ref0,float *rho_ref0,float *d_ref0,float *qs_ref0,float *cvper,int *ncvper,float *uR0,float *uL0,float *cR0,float *cL0,float *rR0, float *rL0, float *ampl0);
void theo_(int *n,float *fbeta,float *h,float *vps,float *qa,float *qb,float *fs,float *din,float *a0,float *c0,float *t0,int *nd,float *rx);
}
// ---------- Lyr_model Class -----------------------------
Lyr_model::Lyr_model(double Thickness, double DepthOfTop,int num_lyr) {
    nlyr = num_lyr; Thk = Thickness;dep_top = DepthOfTop;
    vs.resize(num_lyr,0.); vp.resize(num_lyr,0.0); vpvs.resize(num_lyr,0.0);
    rho.resize(num_lyr,0.0); Qs.resize(num_lyr,0.0); Qp.resize(num_lyr,0.0);
    thk.resize(num_lyr,0.0); depth.resize(num_lyr,0.0);
    cout << "a lyr_model object is constructed: "<< nlyr<<" layers" << endl;
}
Lyr_model::Lyr_model() {
    nlyr = DEFAULTi; Thk = DEFAULTd;dep_top =DEFAULTd;
    vs.clear();vp.clear();vpvs.clear();rho.clear();
    Qs.clear();Qp.clear();T.clear();P.clear();
    thk.clear();depth.clear();
}
Lyr_model::~Lyr_model() {
    //cout << "a lyr_model object is destructed." << endl;
}

int Lyr_model::set_Thk(double Thk) {this->Thk = Thk;return 0;}
int Lyr_model::set_topDep(double topdep) {this->dep_top=topdep;return 0;}
void Lyr_model::set_nlyr(int n) {
    nlyr = n;
    vs.resize(nlyr,DEFAULTd);
    vp.resize(nlyr,DEFAULTd);
    rho.resize(nlyr,DEFAULTd);
    vpvs.resize(nlyr,DEFAULTd);
    Qs.resize(nlyr,DEFAULTd);
    Qp.resize(nlyr,DEFAULTd);
    thk.resize(nlyr,DEFAULTd);
    T.resize(nlyr,DEFAULTd);
    P.resize(nlyr,DEFAULTd);
    depth.resize(nlyr,DEFAULTd);
}

vector<double> Lyr_model::bspline(const vector<double> &v, const vector<AdditionalPert> &add) {
    double dthk=0.,tmp_dep=0.,tmp_top=0.,tmp_bot=0.,tmp_dv=0.;
    vector<double> tmp(nlyr,DEFAULTd);
    //cout<<"bsplie:thk"<<Thk<<endl;
    dthk = Thk/double(nlyr);
    boost::math::cubic_b_spline<double> spline(v.begin(),v.end(),0.,Thk/(v.size()-1.));
    tmp_dep = dthk/2.0;
    for(int i=0;i<nlyr;i++){
        thk[i] = dthk;
        depth[i] = tmp_dep+dep_top;
        tmp[i] = spline(tmp_dep);
        tmp[i] = add_pertb(add,Thk,tmp_dep,tmp[i]);
        tmp_dep = tmp_dep+dthk;
    }
    return tmp;
}
vector<double> Lyr_model::gradient(const vector<double> &v, const vector<AdditionalPert> &add) {
    double dthk=0.,dv=0.,tmp_dep=0.,tmp_top=0.,tmp_bot=0.,tmp_dv=0.;
    vector<double> tmp(nlyr,DEFAULTd);
    dv = (v[1]-v[0])/(nlyr-1.0);
    dthk = Thk/double(nlyr);
    tmp_dep = dthk/2.0;
    for (int i=0;i<nlyr;i++){
        thk[i] = dthk;
        depth[i] = tmp_dep+dep_top; //????
        tmp[i] = v[0] + i * dv;
        tmp[i] = add_pertb(add, Thk, tmp_dep, tmp[i]);
        tmp_dep = tmp_dep + dthk;
    }
    return tmp;
}
vector<double> Lyr_model::layered(const vector<double> &v, const vector<AdditionalPert> &add) {
    double dthk=0.,dv=0.,tmp_vs=0.,tmp_dep=0.,tmp_top=0.,tmp_bot=0.,tmp_dv=0.;
    vector<double> tmp(nlyr,0.);
    dthk = Thk/double(nlyr);
    tmp_dep = dthk/2.0;
    for (int i=0;i<nlyr;i++){
        thk[i] = dthk;
        depth[i] = tmp_dep + dep_top;
        tmp[i] = v[i];
        tmp[i] = add_pertb(add,Thk,tmp_dep,tmp[i]);
        tmp_dep = tmp_dep + dthk;
    }
    return tmp;
}
vector<double> Lyr_model::bulk(const vector<double> &v, const vector<AdditionalPert> &add) {
    if(v.size()!=1){cout<<"update_vp: #vpvs_para is wrong."<<endl;exit(0);}
    double tmp_top=0., tmp_bot = 0., tmp_dv = 0.;
    double tmp_dep = 0.,dthk = 0.;
    vector<double> tmp(nlyr,0.);
    dthk = Thk/double(nlyr);
    tmp_dep = dthk/2.0;
    for(int i=0;i<nlyr;i++){
        thk[i] = dthk;
        depth[i] = tmp_dep+dep_top; //????
        tmp[i] = v[0];
        for (const auto & j : add) {
            tmp_top = j.get_top();
            tmp_bot = j.get_bot();
            tmp_dv = j.get_value();
            if (tmp_dep >= Thk*tmp_top && tmp_dep <= Thk*tmp_bot) { tmp[i] = tmp_dv; }
        }
        tmp_dep = tmp_dep + dthk;
    }
    return tmp;
}

int Lyr_model::update_vs(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_vs) {
    vector<double> tmp(1,0.);
    switch (flag) {
        case 1: // gradient type
            //cout << "group_tp=1, get lyred Vs from gradient group" << endl;
            vs = gradient(v,add_vs);
            break;
        case 2: // layered type
            //cout<<"Group_tp=2, get lyred Vs from lyred group"<<endl;
            vs = layered(v,add_vs);
            break;
        case 3:
            //cout<<"Group_tp=3, get lyred Vs from bs group"<<endl;
            vs = bspline(v,add_vs);
            break;
        case 4: // bulk value
            vs = bulk(v,add_vs);
            break;
        case -1: // water
            //cout<<"flag=0, get lyred Vs for water layer"<<endl;
            tmp[0] = water_vs;
            vs = bulk(tmp,add_vs);
            break;
        case -2: // ice
            //cout<<"flag=-1, get lyred Vs for ice layer"<<endl;
            tmp[0] = ice_vs;
            vs = bulk(tmp,add_vs);
            break;
        default:
            //cout<<"unrecognized flag for group_tp, couldn't update the Vs:"<<flag<<endl;
            return 0; // need to be changed.
    }
    return nlyr;
}
int Lyr_model::update_vp(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_vpvs) {
    vector<double> tmp(1,0.);
    double dthk=0.,dv=0.,tmp_vpvs=0.,tmp_dep=0.,tmp_top=0.,tmp_bot=0.,tmp_dv=0.,tmp_vp=0.;
    dthk = Thk/double(nlyr);
    tmp_dep = dthk/2.0;
    switch (flag) {
        case 1: // gradient type
            //cout << "group_tp=1, get lyred Vp from gradient group" << endl;
            vpvs = gradient(v,add_vpvs);
            break;
        case 2: // layered type
            //cout<<"Group_tp=2, get lyred Vp from lyred group"<<endl;
            vpvs = layered(v,add_vpvs);
            break;
        case 3:
            //cout<<"Group_tp=3, get lyred Vp from bs group"<<endl;
            vpvs = bspline(v,add_vpvs);
            break;
        case 4: //bulk vpvs
            //cout<<"vp_flag=4, update vp using bulk vp/vs"<<endl;
            vpvs = bulk(v,add_vpvs);
            break;
        case -1: // water layer
            //cout<<"vp_tp=0, water lyr"<<endl;
            tmp[0] = water_vp/water_vs;
            vpvs = bulk(tmp,add_vpvs);
            break;
        case -2: // ice layer
            //cout<<"vp_tp = -1, ice layer."<<endl;
            tmp[0] = ice_vp/ice_vs;
            vpvs = bulk(tmp,add_vpvs);
            break;
        case -3: // Brocher
            //cout<<"vp_flag=-4, update Vp using Brocher's"<<endl;
            for (int i=0;i<nlyr;i++){
                tmp_vp = get_vp_brocher(vs[i]);
                vp[i] = add_pertb(add_vpvs,Thk,tmp_dep,tmp_vp);
                this->vpvs[i] = vp[i]/vs[i];
            }break;
        default:
            cout<<"unrecognized flag for vp_tp, couldn't update the Vp"<<endl;
            return 0; // need to be changed
    }
    if(flag>-3){
        for (int i=0;i<nlyr;i++){
            vp[i] = vs[i]*vpvs[i];
        }
    }
    return 1;
}
int Lyr_model::update_rho(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_rho) {
    double tmp_dep,tmp_rho, tmp_top, tmp_bot, tmp_dv;
    switch (flag) {
        case 1:
            //cout<<"Group_tp=1, get lyred rho from gradient group"<<endl;
            rho = gradient(v,add_rho);
            break;
        case 2:
            //cout<<"Group_tp=2, get lyred rho from lyred group"<<endl;
            rho = layered(v,add_rho);
            break;
        case 3:
            //cout<<"Group_tp=3, get lyred rho from bs group"<<endl;
            rho = bspline(v,add_rho);
            break;
        case 4:
            //cout<<"Group_tp=4, get lyred rho from bulk group"<<endl;
            rho = bulk(v,add_rho);
            break;
        case -3: // Brocher?? scale from Vs
            //cout<<"Group_tp=-3, get lyred rho from brocher's group"<<endl;
            for (int i=0;i<nlyr;i++){
                tmp_rho = get_rho_brocher(vs[i],'s');
                tmp_dep=depth[i] - dep_top;
                rho[i] = add_pertb(add_rho,Thk,tmp_dep,tmp_rho);
            }break;
        case -4: //hacker, usually for mantle
            //cout<<"Group_tp=-4, get lyred rho from hacker's group"<<endl;
            for (int i=0;i<nlyr;i++){
                tmp_rho = get_rho_hacker(vs[i]);
                tmp_dep = depth[i] - dep_top;
                rho[i] = add_pertb(add_rho,Thk,tmp_dep,tmp_rho);
            }break;
		case -5: // Brocher scale from Vp
			for(int i=0;i<nlyr;i++){
				tmp_rho = get_rho_brocher(vp[i],'p');
				tmp_dep=depth[i] - dep_top;
                rho[i] = add_pertb(add_rho,Thk,tmp_dep,tmp_rho);
            }break;
        case -1: // water lyr
            //cout<<"Group_tp=-1, get lyred rho from water group"<<endl;
            for (int i=0;i<nlyr;i++) {
                rho[i] = water_rho;
            }break;
        case -2: // ice layer
            //cout<<"Group_tp=-2, get lyred rho from ice group"<<endl;
            for (int i=0;i<nlyr;i++){
                rho[i]= ice_rho;
            }break;
        default:
            //cout<<"rho_tp is wrong!"<<endl;
            return 0;
    }
    return 1;
}
int Lyr_model::update_Qs(const int &flag,const vector<double> &v,const vector<AdditionalPert> &add_Qs) {
    vector<double> tmp(1,0.);
    switch (flag) {
        case 1:
            Qs = gradient(v,add_Qs);
            break;
        case 2:
            Qs = layered(v,add_Qs);
            break;
        case 3:
            Qs = bspline(v,add_Qs);
            break;
        case 4:
            Qs = bulk(v,add_Qs);
            break;
        case -1:
            tmp[0] = water_Qs;
            Qs = bulk(tmp,add_Qs);
            break;
        case -2:
            tmp[0] = ice_Qs;
            Qs = bulk(tmp,add_Qs);
            break;
        case -3:
            tmp[0] = sedi_Qs;
            Qs = bulk(tmp,add_Qs);
            break;
        case -4:
            tmp[0] = cru_Qs;
            Qs = bulk(tmp,add_Qs);
            break;
        case -5:
            tmp[0] = man_Qs;
            Qs = bulk(tmp,add_Qs);
            break;
        default:
            cout<<"update Qs: wrong flag"<<endl;
            return 1;
    }
    return 0;
}
int Lyr_model::update_Qp(const int &flag,const vector<double> &v,const vector<AdditionalPert> &add_Qp) {
    vector<double> tmp(1,0.);
    switch (flag) {
        case 1:
            Qp = gradient(v,add_Qp);
            break;
        case 2:
            Qp = layered(v,add_Qp);
            break;
        case 3:
            Qp = bspline(v,add_Qp);
            break;
        case 4:
            Qp = bulk(v,add_Qp);
            break;
        case -1:
            tmp[0] = water_Qp;
            Qp = bulk(tmp,add_Qp);
            break;
        case -2:
            tmp[0] = ice_Qp;
            Qp = bulk(tmp,add_Qp);
            break;
        case -3: // sediment Qp
            tmp[0] = sedi_Qp;
            Qp = bulk(tmp,add_Qp);
            break;
        case -4:
            tmp[0] = cru_Qp;
            Qp = bulk(tmp,add_Qp);
            break;
        case -5:
            tmp[0] = man_Qp;
            Qp = bulk(tmp,add_Qp);
            break;
        default:
            cout<<"update Qp: wrong flag"<<endl;
            return 1;
    }
    return 0;
}
int Lyr_model::update_T(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_T) {

    switch (flag) {
        case 1:
            T = gradient(v,add_T);
            break;
        case 2:
            T = layered(v,add_T);
            break;
        case 3:
            T = bspline(v,add_T);
            break;
        case 4:
            T = bulk(v, add_T);
            break;
        default:
            cout<<"update T: wrong flag:"<<flag<<endl;
            return 1;
    }
    return 0;
}
int Lyr_model::update_P(const int &flag, const vector<double> &v, const vector<AdditionalPert> &add_P) {
    switch (flag) {
        case 1:
            P = gradient(v,add_P);
            break;
        case 2:
            P = layered(v,add_P);
            break;
        case 3:
            P = bspline(v,add_P);
            break;
        case 4:
            P = bulk(v, add_P);
            break;
        default:
            cout<<"update P: wrong flag:"<<flag<<endl;
            return 1;
    }
    return 0;
}

int Lyr_model::get_nlyr() const {return nlyr;}
double Lyr_model::get_Thk() const {return Thk;}
double Lyr_model::get_dep_top() const {return dep_top;}
vector<double> Lyr_model::get_vs() const {return vs;}
vector<double>Lyr_model::get_vp() const {return vp;}
vector<double> Lyr_model::get_vpvs() const {return vpvs;}
vector<double> Lyr_model::get_rho() const {return rho;}
vector<double> Lyr_model::get_Qs() const {return Qs;}
vector<double> Lyr_model::get_Qp() const {return Qp;}
vector<double> Lyr_model::get_thk() const {return thk;}
vector<double> Lyr_model::get_depth() const {return depth;}
vector<double> Lyr_model::get_T() const {return T;}
vector<double> Lyr_model::get_P() const {return P;}

int Lyr_model::set_thkdep(const vector<double> &thk, const vector<double> &dep) {
    this->thk = thk; this->depth = dep;return 0;
}
int Lyr_model::set_vs(const vector<double> &vs) {this->vs = vs;return 0;}
int Lyr_model::set_vp(const vector<double> &vp) {this->vp=vp;return 0;}
int Lyr_model::set_vpvs(const vector<double> &vpvs){this->vpvs=vpvs;return 0;}
int Lyr_model::set_rho(const vector<double> &rho) {this->rho = rho;return 0;}
int Lyr_model::set_Q(const vector<double> &Qs, const vector<double> &Qp) {
    this->Qs = Qs; this->Qp = Qp;return 0;
}

void Lyr_model::print_lyr() const { //to console
    cout<<depth.size()<<endl;
    for (int i=0;i<nlyr;i++){
        cout<<depth[i]<<" "<<vs[i]<<" "<<vp[i]<<" "<<rho[i]<<" ";
        cout<<Qs[i]<<" "<<Qp[i]<<" "<<thk[i]<<" "<<endl;
    }
}
void Lyr_model::write_to_txt(std::ofstream &f) const {
    for(int i=0;i<nlyr;i++){
        f<<depth[i]<<" "<<vs[i]<<" "<<vp[i]<<" "<<rho[i]<<" ";
        f<<Qs[i]<<" "<<Qp[i]<<" "<<thk[i]<<endl;
    }
	cout<<"Lyr_model::nlyr:"<<nlyr<<endl;
}
/*
vector<double> Lyr_model::get_T() {
    return T;
}
vector<double> Lyr_model::get_P() {
    return P;
}
vector<double> Lyr_model::get_porosity() {
    return porosity;
}
*/
/*
int Lyr_model::update_info() {
    nlyr = vs.size();
    dep_top=depth[0]-thk[0]/2;
    Thk = depth[-1] - dep_top + thk[-1]/2;
    return 0;
}*/
// ------------ END of Lyr_model Class ---------------

//-------------- Group Class ---------------------
AdditionalPert::AdditionalPert() {top=0.;bot=0.;value=0.;}
AdditionalPert::AdditionalPert(double top, double bot, double value) {
    this->top = top; this->bot = bot; this->value = value;
}
AdditionalPert::~AdditionalPert() = default;

int AdditionalPert::set(double top, double bot, double value) {
    this->top = top;this->bot = bot;this->value=value;
    return 0;
}
void AdditionalPert::set_top(double top) {this->top = top;}
void AdditionalPert::set_bot(double bot) {this->bot = bot;}
void AdditionalPert::set_value(double value) {this->value = value;}

double AdditionalPert::get_top() const {return top;}
double AdditionalPert::get_bot() const {return bot;}
double AdditionalPert::get_value() const {return value;}

AdditionalPert AdditionalPert::operator+(const AdditionalPert &a) const {
    return {top+a.top,bot+a.bot,value+a.value};
}
AdditionalPert AdditionalPert::operator-(const AdditionalPert &a) const {
    return {top-a.top,bot-a.bot,value-a.value};
}
AdditionalPert AdditionalPert::operator/(int n) {
    return {top/n,bot/n,value/n};
}

Property::Property() {
    para.clear();n_para = DEFAULTi;tp = DEFAULTi;
    thk = DEFAULTd;topdep = DEFAULTd;
    add_pert.clear();
    n_pert = DEFAULTi;
}
Property::~Property() = default;
Property::Property(const int &tp, const vector<double> &para, const int &n_para, const vector<AdditionalPert> &add_pert, const int &n_pert, const double &thk,
                   const double &dep) {
    this->tp = tp;
    this->n_para = n_para; this->para = para;
    this->n_pert = n_pert; this->add_pert = add_pert;
    this->thk = thk; this->topdep = dep;
}
vector<double> Property::get_para() const {return para;}
int Property::set(int &index, const vector<std::string> &v) {
    int nlyr;
    tp = atoi(v[index].c_str()); index++; // type of group(how to interp parameter)
    thk = atof(v[index].c_str()); index++;
    n_para = atoi(v[index].c_str()); index++;
    para.resize(n_para);
    for(int j=0;j<n_para;j++){
        para[j] = atof(v[index].c_str()); index++;
    }
    //cout<<"done:para"<<endl;
    n_pert = atoi(v[index].c_str()); index++;
    add_pert.resize(n_pert);
    double top,bot,value;
    for (int j=0;j<n_pert;j++){
        top = atof(v[index].c_str()); index++;
        bot = atof(v[index].c_str()); index++;
        value = atof(v[index].c_str()); index++;
        AdditionalPert tmp(top,bot,value);
        add_pert[j] = tmp;
    }
    //cout<<"done:add_para"<<endl;
    nlyr = atoi(v[index].c_str()); index++;
    if(index==v.size()-1) topdep = atof(v[index].c_str());
    else topdep = DEFAULTd;
    return nlyr;
}
void Property::print() const {
    cout<<tp<<" "<<n_para<<" ";
    for(int i=0;i<n_para;i++){
        cout<<para[i]<<" ";
    }
    cout<<n_pert<<" ";
    for (int i=0;i<n_pert;i++){
        cout<<add_pert[i].get_top()<<" "<<add_pert[i].get_bot()<<" "<<add_pert[i].get_value()<<" ";
    }
    cout<<thk<<" "<<topdep<<" ";
}
void Property::write_to_bin(std::ofstream &f) const {
    double tmp_d;
    int tmp_i;
    tmp_i = tp;
    f.write((char *)&tmp_i,SizeOfInt);
    //cout<<"Property::bin:tp"<<tmp_i<<endl;
    tmp_d = thk;
    f.write((char *)&tmp_d, SizeOfDouble);
    f.write((char *)&n_para,SizeOfInt);
    //cout<<"Property::bin:n_para"<<n_para<<endl;
    for(int i=0;i<n_para;i++){
        f.write((char *)&para[i],SizeOfDouble);
    }
    //cout<<"Property::bin:para"<<endl;
    f.write((char *)&n_pert,SizeOfInt);
    //cout<<"Property::bin:n_pert"<<n_pert<<endl;

    for(int i=0;i<n_pert;i++){
        tmp_d = add_pert[i].get_top();
        f.write((char *)&tmp_d, SizeOfDouble);
        tmp_d = add_pert[i].get_bot();
        f.write((char *)&tmp_d, SizeOfDouble);
        tmp_d = add_pert[i].get_value();
        f.write((char *)&tmp_d, SizeOfDouble);
    }
}
void Property::read_from_bin(std::ifstream &f) {
    double tmp_d;
    //int tmp_i;
    f.read((char *)&tp,SizeOfInt);
    //cout<<"Property:tp"<<tp<<endl;
    f.read((char *)&thk, SizeOfDouble);
    f.read((char *)&n_para,SizeOfInt);
    para.resize(n_para,DEFAULTd);
    for(int i=0;i<n_para;i++){
        f.read((char *)&para[i],SizeOfDouble);
    }
    f.read((char *)&n_pert,SizeOfInt);
    add_pert.resize(n_pert);
    for(int i=0;i<n_pert;i++){
        f.read((char *)&tmp_d,SizeOfDouble);
	add_pert[i].set_top(tmp_d);
        f.read((char *)&tmp_d,SizeOfDouble);
	add_pert[i].set_bot(tmp_d);
        f.read((char *)&tmp_d,SizeOfDouble);
	add_pert[i].set_value(tmp_d);
    }
    topdep = DEFAULTd; // todo: not used
}

int Group::count=0;
Group::~Group() = default;
Group::Group() {
    thk=DEFAULTd;dep=DEFAULTd;P0=DEFAULTd;
    nlyr = DEFAULTi;idx=DEFAULTi;
}; //todo :
Group::Group(const int &n_group) {
    idx = count; // index of group
    count++;
    P0 = DEFAULTd;
    if(count == n_group) count =0;
    nlyr = DEFAULTi; thk = DEFAULTd; dep = DEFAULTd;
}
Group::Group(const int &n_group, const int &nlyr, const Property &vs, const Property &vpvs, const Property &rho, const Property &Qs, const Property &Qp, const Property &T, const Property &P) {
    idx = count;
    count++;
    if(count == n_group) count = 0;
    this->nlyr = nlyr;
    this->vs = vs; this->vpvs = vpvs; this->rho = rho;
    this->Qs = Qs; this->Qp = Qp; this->T = T; this->P = P;
    thk = DEFAULTd;P0 = DEFAULTd;dep = DEFAULTd;

}

int Group::readin(const vector<string> &v) {
    int index = 0;// index of element in v
    int np;
    int idx = atoi(v[index].c_str()); index++;// index of group
    if (this->idx!=idx){cout<<"sth wrong with the index"<<endl;exit(0);} //????
    int property_tp = atoi(v[index].c_str()); index++; // type of property
    switch (property_tp) {
        case 1:
            nlyr = vs.set(index,v); //????
            thk = vs.thk;
            break;
        case 2:
            vpvs.set(index,v);
            break;
        case 3:
            rho.set(index, v);
            break;
        case 4:
            Qs.set(index,v);
            break;
        case 5:
            //cout<<"set Qp";
            Qp.set(index,v);
            //cout<<"done"<<endl;
            break;
        case 6:
            //cout<<"set T ";
            T.set(index, v);
            //cout<<"done"<<endl;
            break;
        case 7:
            P.set(index, v);
            break;
        default:
            cout<<"unrecognized property type"<<endl;
            return 1;
    }
    return 0;
}
void Group::set_nlyr(const int &n) {nlyr=n;}
void Group::set_dep() {
    dep = vs.topdep; // todo: to be further considered
    /*
    while(true){
        if(vs.topdep>=0) {dep = vs.topdep; break;}
        if(vpvs.topdep>=0) {dep = vpvs.topdep;break;}
        if(rho.topdep>=0){dep = rho.topdep;break;}
        if(Qs.topdep>=0){dep = Qs.topdep;break;}
        if(Qp.topdep>=0){dep = Qp.topdep;break;}
        if(T.topdep>=0){dep = T.topdep;break;}
        if(P.topdep>=0){dep = P.topdep;break;}
    }*/
}
void Group::set_dep(const double &dep) {this->dep = dep;}

int Group::update_from_sample(const Sample &s) {
    vs.para = s.get_vs(); vpvs.para = s.get_vpvs(); rho.para = s.get_rho();
    vs.thk = s.get_thk(); vpvs.thk = s.get_thk(); rho.thk = s.get_thk();
    Qs.para = s.get_Qs(); Qp.para = s.get_Qp(); T.para = s.get_T(); P.para = s.get_P();
    Qs.thk = s.get_thk(); Qp.thk = s.get_thk(); T.thk = s.get_thk();P.thk = s.get_thk();
    vs.add_pert = s.get_addvs(); vpvs.add_pert = s.get_addvpvs(); rho.add_pert = s.get_addrho();
    Qs.add_pert = s.get_addQs(); Qp.add_pert = s.get_addQp(); T.add_pert = s.get_addT(); P.add_pert = s.get_addP();
    thk = s.get_thk();
    return 1;
}
int Group::update_from_sample(const Sample &s, const double thickness) {
    vs.para = s.get_vs(); vpvs.para = s.get_vpvs(); rho.para = s.get_rho();
    vs.thk = s.get_thk(); vpvs.thk = s.get_thk(); rho.thk = s.get_thk();
    Qs.para = s.get_Qs(); Qp.para = s.get_Qp(); T.para = s.get_T(); P.para = s.get_P();
    Qs.thk = s.get_thk(); Qp.thk = s.get_thk(); T.thk = s.get_thk(); P.thk = s.get_thk();
    vs.add_pert = s.get_addvs(); vpvs.add_pert = s.get_addvpvs(); rho.add_pert = s.get_addrho();
    Qs.add_pert = s.get_addQs(); Qp.add_pert = s.get_addQp(); T.add_pert = s.get_addT(); P.add_pert = s.get_addP();
    thk = thickness;
    return 1;
}
void Group::update_lyrmodel(){
    lyrm.set_Thk(thk);
    lyrm.set_topDep(dep);
    lyrm.set_nlyr(nlyr);
    //???? sequence?
    lyrm.update_vs(vs.tp,vs.para,vs.add_pert);
    lyrm.update_vp(vpvs.tp,vpvs.para,vpvs.add_pert);
    lyrm.update_rho(rho.tp,rho.para,rho.add_pert);
    lyrm.update_Qs(Qs.tp,Qs.para,Qs.add_pert);
    lyrm.update_Qp(Qp.tp,Qp.para,Qp.add_pert);
    lyrm.update_T(T.tp,T.para,T.add_pert);
    lyrm.update_P(P.tp,P.para,P.add_pert);
}

vector<double> Group::get_lyr_vs() const {return lyrm.get_vs();}
vector<double> Group::get_lyr_vp() const {return lyrm.get_vp();}
vector<double> Group::get_lyr_rho() const {return lyrm.get_rho();}
vector<double> Group::get_lyr_dep() const {return lyrm.get_depth();}
vector<double> Group::get_lyr_Qs() const {return lyrm.get_Qs();}
vector<double> Group::get_lyr_Qp() const {return lyrm.get_Qp();}
vector<double> Group::get_lyr_thk() const {return lyrm.get_thk();}
vector<double> Group::get_lyr_vpvs() const {return lyrm.get_vpvs();}

double Group::get_thk() const {return thk;};
void Group::print_info() const {
    cout<<"index: "<<idx<<endl;
    cout<<"count: "<<count<<endl;
    cout<<"nlyr:" << nlyr<<endl;
    cout<<"thk:" << thk<<endl;
}
int Group::get_nlyr()  const {return nlyr;}
double Group::get_dep()  const {return dep;}
void Group::write_lyr()  const { lyrm.print_lyr();}
void Group::print() const {
    cout<<idx<<" ";vs.print();cout<<nlyr<<endl;
    cout<<idx<<" ";vpvs.print();cout<<nlyr<<endl;
    cout<<idx<<" ";rho.print();cout<<nlyr<<endl;
    cout<<idx<<" ";Qs.print();cout<<nlyr<<endl;
    cout<<idx<<" ";Qp.print();cout<<nlyr<<endl;
    cout<<idx<<" ";T.print();cout<<nlyr<<endl;
    cout<<idx<<" ";P.print();cout<<nlyr<<endl;
}
void Group::write_v_to_file(ofstream &f) const {
	f<<thk<<" ";
	for(int i=0;i<vs.n_para;i++) f<<vs.para[i]<<" ";
	for(int i=0;i<vpvs.n_para;i++) f<<vpvs.para[i]<<" ";
	for(int i=0;i<vpvs.n_pert;i++) f<<vpvs.add_pert[i].get_value()<<" ";
	f<<endl;
}
void Group::write_to_file(ofstream &f) const {
    f<<idx<<" 1 "<<vs.tp<<" "<<vs.thk<<" "<<vs.n_para<<" ";
    for(int i=0;i<vs.n_para;i++) f<<vs.para[i]<<" ";
    f<<vs.n_pert<<" ";
    for(int i=0;i<vs.n_pert;i++) f<<vs.add_pert[i].get_top()<<" "<<vs.add_pert[i].get_bot()<<" "<<vs.add_pert[i].get_value()<<" ";
    f<<nlyr<<" ";
    if(idx==0) f<<vs.topdep;
    f<<"\n";
    //vpvs
    f<<idx<<" 2 "<<vpvs.tp<<" "<<vpvs.thk<<" "<<vpvs.n_para<<" ";
    for(int i=0;i<vpvs.n_para;i++) f<<vpvs.para[i]<<" ";
    f<<vpvs.n_pert<<" ";
    for(int i=0;i<vpvs.n_pert;i++) f<<vpvs.add_pert[i].get_top()<<" "<<vpvs.add_pert[i].get_bot()<<" "<<vpvs.add_pert[i].get_value()<<" ";
    f<<nlyr<<"\n";
    //rho
    f<<idx<<" 3 "<<rho.tp<<" "<<rho.thk<<" "<<rho.n_para<<" ";
    for(int i=0;i<rho.n_para;i++) f<<rho.para[i]<<" ";
    f<<rho.n_pert<<" ";
    for(int i=0;i<rho.n_pert;i++) f<<rho.add_pert[i].get_top()<<" "<<rho.add_pert[i].get_bot()<<" "<<rho.add_pert[i].get_value()<<" ";
    f<<nlyr<<"\n";
    //Qs
    f<<idx<<" 4 "<<Qs.tp<<" "<<Qs.thk<<" "<<Qs.n_para<<" ";
    for(int i=0;i<Qs.n_para;i++) f<<Qs.para[i]<<" ";
    f<<Qs.n_pert<<" ";
    for(int i=0;i<Qs.n_pert;i++) f<<Qs.add_pert[i].get_top()<<" "<<Qs.add_pert[i].get_bot()<<" "<<Qs.add_pert[i].get_value()<<" ";
    f<<nlyr<<"\n";
    //Qp
    f<<idx<<" 5 "<<Qp.tp<<" "<<Qp.thk<<" "<<Qp.n_para<<" ";
    for(int i=0;i<Qp.n_para;i++) f<<Qp.para[i]<<" ";
    f<<Qp.n_pert<<" ";
    for(int i=0;i<Qp.n_pert;i++) f<<Qp.add_pert[i].get_top()<<" "<<Qp.add_pert[i].get_bot()<<" "<<Qp.add_pert[i].get_value()<<" ";
    f<<nlyr<<"\n";
    //T
    f<<idx<<" 6 "<<T.tp<<" "<<T.thk<<" "<<T.n_para<<" ";
    for(int i=0;i<T.n_para;i++) f<<T.para[i]<<" ";
    f<<T.n_pert<<" ";
    for(int i=0;i<T.n_pert;i++) f<<T.add_pert[i].get_top()<<" "<<T.add_pert[i].get_bot()<<" "<<T.add_pert[i].get_value()<<" ";
    f<<nlyr<<"\n";
    //P
    f<<idx<<" 7 "<<P.tp<<" "<<P.thk<<" "<<P.n_para<<" ";
    for(int i=0;i<P.n_para;i++) f<<P.para[i]<<" ";
    f<<P.n_pert<<" ";
    for(int i=0;i<P.n_pert;i++) f<<P.add_pert[i].get_top()<<" "<<P.add_pert[i].get_bot()<<" "<<P.add_pert[i].get_value()<<" ";
    f<<nlyr<<"\n";
}
void Group::write_to_bin(std::ofstream &f) const {
    f.write((char *)&idx,SizeOfInt);
    //cout<<"write idx"<<tmp_i<<endl;
    vs.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    f.write((char *)&dep,SizeOfDouble);
    vpvs.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    rho.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    Qs.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    Qp.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    T.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    P.write_to_bin(f);
    f.write((char *)&nlyr,SizeOfInt);
    //cout<<"done:Group::writeBin"<<endl;
}
void Group::read_from_bin(std::ifstream &f) {
    int tmp;
    f.read((char *)&idx,SizeOfInt);
    //cout<<"idx:"<<tmp<<endl;
    vs.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    f.read((char *)&dep,SizeOfDouble);
    vpvs.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    rho.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    Qs.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    Qp.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    T.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    P.read_from_bin(f);
    f.read((char *)&nlyr,SizeOfInt);
    thk = vs.thk; // todo: consistent with main code
    //todo: P0 need be added
}
// --------------end of Group Class --------------------------

//------------------ Sample Class ------------------------------
subSample::subSample() {
    para.clear();left.clear();right.clear();step.clear();
    add_para.clear();add_left.clear();add_right.clear();add_step.clear();
    n_pert = 0;n_para = 0;
}
subSample::~subSample() = default;
double subSample::Initial(int &index, const vector<string> &v) {
    int tp = atoi(v[index].c_str()); index++; // type of interp
    double thk = atof(v[index].c_str()); index++;
    n_para = atoi(v[index].c_str()); index++;
    if (n_para>0) {para.resize(n_para);step.resize(n_para,0.);}
    else {para.clear();step.clear();}
    for(int j=0;j<n_para;j++){
        para[j] = atof(v[index].c_str()); index++;
    }
    left = para;right = para;
    //left.resize(n_para,0.);
    //right.resize(n_para,0.);

    n_pert = atoi(v[index].c_str()); index++;
    if(n_pert>0) {add_para.resize(n_pert); add_step.resize(n_pert);}
    else {add_para.clear();add_step.clear();}
    double top,bot,value;
    for (int j=0;j<n_pert;j++){
        top = atof(v[index].c_str()); index++;
        bot = atof(v[index].c_str()); index++;
        value = atof(v[index].c_str()); index++;
        AdditionalPert tmp(top,bot,value);
        add_para[j] = tmp;
    }
    //add_left.resize(n_pert);add_right.resize(n_pert);
    add_left = add_para;add_right = add_para;
    return thk;
}
//int subSample::readin(int &index, const vector<std::string> &v) {}

void subSample::set(const int &para_index, const double &radius, const double &step) {
    left[para_index] = para[para_index] - radius;
    right[para_index] = para[para_index] + radius;
    this->step[para_index] = step;
}
void subSample::set(const char &f, const int &para_index, const double &radius, const double &step) {
    switch (f) {
        double top, bot, value;
        case 't':
            top = add_para[para_index].get_top();
            add_left[para_index].set_top(top-radius);
            add_right[para_index].set_top(top+radius);
            add_step[para_index].set_top(step);
            break;
        case 'b':
            bot = add_para[para_index].get_bot();
            add_left[para_index].set_bot(bot-radius);
            add_right[para_index].set_bot(bot+radius);
            add_step[para_index].set_bot(step);
            break;
        case 'v':
            value = add_para[para_index].get_value();
            add_left[para_index].set_value(value-radius);
            add_right[para_index].set_value(value+radius);
            add_step[para_index].set_value(step);
            break;
        default:
            cout<<"wrong flag for set addtional perturbation"<<endl;

    }
}
void subSample::uniform_rand() {
    for(int i=0;i<para.size();i++){
        if(step[i]>0.) para[i] = rand_uniform(left[i],right[i],para[i]);
    }
    for(int i=0;i<add_para.size();i++){
        add_para[i] = rand_uniform_add(add_left[i], add_right[i], add_para[i]);
    }
}
void subSample::gaussian_rand() {
    for(int i=0;i<para.size();i++){
        para[i] = rand_gaussian(para[i],step[i],left[i],right[i]);
        //cout<<"step: "<<step[i]<<endl;
    }
    for(int i=0;i<add_para.size();i++){
        add_para[i] = rand_gaussian_add(add_para[i], add_step[i], add_left[i], add_right[i]);
    }
    //cout<<"done:subSample gaussian rand"<<endl;
}
vector<double> subSample::get_para() const {return para;}
void subSample::print_sample() const {
    for(int i=0;i<n_para;i++) cout<<para[i]<<" +- "<<para[i]-left[i]<<"\n";
    for(int i=0;i<n_pert;i++) {
        cout<<add_para[i].get_top()<<" +- "<<add_para[i].get_top()-add_left[i].get_top()<<" ";
        cout<<add_para[i].get_bot()<<" +- "<<add_para[i].get_bot()-add_left[i].get_bot()<<" ";
        cout<<add_para[i].get_value()<<" +- "<<add_para[i].get_value()-add_left[i].get_value()<<"\n";
    }
    if(n_pert==0 and n_para==0) cout<<"\n";
}
void subSample::print() const {
    for (int i=0;i<para.size();i++){
        if(step[i]>0) cout<<para[i]<<" ";
    }
    for (int i=0;i<add_para.size();i++) {
        if(add_step[i].get_value()>0) cout<<add_para[i].get_top()<<" "<<add_para[i].get_bot()<<" "<<add_para[i].get_value()<<" ";
    }
}
vector<AdditionalPert> subSample::get_addpara() const {return add_para;}

Sample::Sample() {
    thk = DEFAULTd; step_thk = DEFAULTd;left_thk = DEFAULTd; right_thk = DEFAULTd;
}
Sample::~Sample() = default;
Sample::Sample(subSample vs, subSample vpvs, subSample rho, subSample Qs, subSample Qp, subSample T, subSample P, double thk, double left_thk, double right_thk, double step_thk){
    this->vs = vs; this->vpvs = vpvs; this->rho;
    this->Qs = Qs; this->Qp = Qp; this->T = T; this->P = P;
    this->thk = thk; this->left_thk = left_thk; this->right_thk = right_thk; this->step_thk = step_thk;
}
void Sample::set(const subSample &vs, const subSample &vpvs, const subSample &rho, const subSample &Qs, const subSample &Qp, const subSample &T, const subSample &P,
                 const double &thk, const double &left_thk, const double &right_thk, const double &step_thk) {
    this->vs = vs; this->vpvs = vpvs;this->rho = rho;
    this->Qs = Qs; this->Qs; this->Qp = Qp; this->T =T; this->P =P;
    this->thk = thk; this->left_thk = left_thk;this->right_thk = right_thk;this->step_thk = step_thk;
}
//initial from .mod file
int Sample::Initial(const vector<string> &v) {
    int index = 0;// index of element in v
    int idx = atoi(v[index].c_str());index++; // index of group
    int np;
    int property_tp = atoi(v[index].c_str()); index++; // type of property
    switch (property_tp) {
        case 1:
            cout<<"Sample"<<idx<<":initial vs"<<endl;
            thk = vs.Initial(index, v);
            break;
        case 2:
            cout<<"Sample"<<idx<<":initial vpvs"<<endl;
            vpvs.Initial(index, v);
            break;
        case 3:
            cout<<"Sample"<<idx<<":initial rho"<<endl;
            rho.Initial(index, v);
            break;
        case 4:
            cout<<"Sample"<<idx<<":initial Qs"<<endl;
            Qs.Initial(index, v);
            break;
        case 5:
            cout<<"Sample"<<idx<<":initial Qp"<<endl;
            Qp.Initial(index, v);
            break;
        case 6:
            cout<<"Sample"<<idx<<":initial T"<<endl;
            T.Initial(index, v);
            break;
        case 7:
            cout<<"Sample"<<idx<<":initial P"<<endl;
            P.Initial(index, v);
            break;
        default:
            cout<<"unrecognized property type"<<endl;
            return 1;
    }
    return 0;
}

//read in from in.para file
int Sample::readin(const vector<string> &v){
    int index = 0;
    int idx = atoi(v[0].c_str()); index++; // index of sample
    int property_tp = atoi(v[index].c_str()); index++;
    int radius_tp = atoi(v[index].c_str()); index++; // 0:percentage
    double radius = atof(v[index].c_str()); index++;
    double step = atof(v[index].c_str()); index++;
    int para_idx = DEFAULTi;
    double top=0.,bot=0.,value=0.;
    if(property_tp!=0) {para_idx = atoi(v[index].c_str()); index++;}
    switch (property_tp) {
        case 0:
            if(radius_tp==0) radius = abs(thk*radius);
            cout<<"thk radius:"<<radius<<endl;
            left_thk = thk - radius; right_thk = thk + radius;
            if(left_thk<0.) left_thk=0.;
			cout<<"thk left:"<<left_thk<<endl;
            step_thk = step;
            break;
        case 1:
            if(radius_tp==0) radius = abs(vs.para[para_idx]*radius);
            vs.set(para_idx,radius,step);
            break;
        case 2:
            if(radius_tp==0) radius = abs(vpvs.para[para_idx]*radius);
            vpvs.set(para_idx,radius,step);
            break;
        case 3:
            if(radius_tp==0) radius = abs(rho.para[para_idx]*radius);
            rho.set(para_idx,radius,step);
            break;
        case 4:
            if(radius_tp==0) radius = abs(Qs.para[para_idx]*radius);
            Qs.set(para_idx,radius,step);
            break;
        case 5:
            if(radius_tp==0) radius = abs(Qp.para[para_idx]*radius);
            Qp.set(para_idx,radius,step);
            break;
        case 6:
            if(radius_tp==0) radius = abs(T.para[para_idx]*radius);
            T.set(para_idx,radius,step);
            break;
        case 7:
            if(radius_tp==0) radius = abs(P.para[para_idx]*radius);
            P.set(para_idx,radius,step);
            break;
        case -10:
            if(radius_tp==0) radius = abs(vs.add_para[para_idx].get_top()*radius);
            vs.set('t',para_idx,radius,step);
            break;
        case -11:
            if(radius_tp==0) radius = abs(vs.add_para[para_idx].get_bot()*radius);
            vs.set('b',para_idx,radius,step);
            break;
        case -12:
            if(radius_tp==0) radius = abs(vs.add_para[para_idx].get_value()*radius);
            vs.set('v',para_idx,radius,step);
            break;
        case -20:
            if(radius_tp==0) radius = abs(vpvs.add_para[para_idx].get_top()*radius);
            vpvs.set('t',para_idx,radius,step);
            break;
        case -21:
            if(radius_tp==0) radius = abs(vpvs.add_para[para_idx].get_bot()*radius);
            vpvs.set('b',para_idx,radius,step);
            break;
        case -22:
            if(radius_tp==0) radius = abs(vpvs.add_para[para_idx].get_value()*radius);
            vpvs.set('v',para_idx,radius,step);
            break;
        case -30:
            if(radius_tp==0) radius = abs(rho.add_para[para_idx].get_top()*radius);
            rho.set('t',para_idx,radius,step);
            break;
        case -31:
            if(radius_tp==0) radius = abs(rho.add_para[para_idx].get_bot()*radius);
            rho.set('b',para_idx,radius,step);
            break;
        case -32:
            if(radius_tp==0) radius = abs(rho.add_para[para_idx].get_value()*radius);
            rho.set('v',para_idx,radius,step);
            break;
        case -40:
            if(radius_tp==0) radius = abs(Qs.add_para[para_idx].get_top()*radius);
            Qs.set('t',para_idx,radius,step);
            break;
        case -41:
            if(radius_tp==0) radius = abs(Qs.add_para[para_idx].get_bot()*radius);
            Qs.set('b',para_idx,radius,step);
            break;
        case -42:
            if(radius_tp==0) radius = abs(Qs.add_para[para_idx].get_value()*radius);
            Qs.set('v',para_idx,radius,step);
            break;
        case -50:
            if(radius_tp==0) radius = abs(Qp.add_para[para_idx].get_top()*radius);
            Qp.set('t',para_idx,radius,step);
            break;
        case -51:
            if(radius_tp==0) radius = abs(Qp.add_para[para_idx].get_bot()*radius);
            Qp.set('b',para_idx,radius,step);
            break;
        case -52:
            if(radius_tp==0) radius = abs(Qp.add_para[para_idx].get_value()*radius);
            Qp.set('v',para_idx,radius,step);
            break;
        case -60:
            if(radius_tp==0) radius = abs(T.add_para[para_idx].get_top()*radius);
            T.set('t',para_idx,radius,step);
            break;
        case -61:
            if(radius_tp==0) radius = abs(T.add_para[para_idx].get_bot()*radius);
            T.set('b',para_idx,radius,step);
            break;
        case -62:
            if(radius_tp==0) radius = abs(T.add_para[para_idx].get_value()*radius);
            T.set('v',para_idx,radius,step);
            break;
        case -70:
            if(radius_tp==0) radius = abs(P.add_para[para_idx].get_top()*radius);
            P.set('t',para_idx,radius,step);
            break;
        case -71:
            if(radius_tp==0) radius = abs(P.add_para[para_idx].get_bot()*radius);
            P.set('b',para_idx,radius,step);
            break;
        case -72:
            if(radius_tp==0) radius = abs(P.add_para[para_idx].get_value()*radius);
            P.set('v',para_idx,radius,step);
            break;
        default:
            cout<<"Sample:proper_type is wrong"<<endl;
            return 1;
        }
    return 0;
}

int Sample::rand_start_sample() {
    vs.uniform_rand();
    vpvs.uniform_rand();
    rho.uniform_rand();
    Qs.uniform_rand(); Qp.uniform_rand();
    T.uniform_rand();
    P.uniform_rand();
    thk = rand_uniform(left_thk,right_thk,thk);
    return 0;
}
int Sample::gen_new_sample() {
    //cout<<"Sample: gen_new_sample"<<endl;
    vs.gaussian_rand();vpvs.gaussian_rand();rho.gaussian_rand();
    Qs.gaussian_rand();Qp.gaussian_rand();T.gaussian_rand();P.gaussian_rand();
    thk = rand_gaussian(thk,step_thk,left_thk,right_thk);
    return 0;
}
void Sample::set_thk(const double &thk) {this->thk = thk;}

double Sample::get_thk() const {return thk;}
vector<double> Sample::get_vs() const {return vs.para;}
vector<AdditionalPert> Sample::get_addvs() const {return vs.add_para;}
vector<double> Sample::get_vpvs() const {return vpvs.para;}
vector<AdditionalPert> Sample::get_addvpvs() const {return vpvs.add_para;}
vector<double> Sample::get_rho() const {return rho.para;}
vector<AdditionalPert> Sample::get_addrho() const {return rho.add_para;}
vector<double> Sample::get_P() const {return P.para;}
vector<AdditionalPert> Sample::get_addP() const {return P.add_para;}
vector<double> Sample::get_Qs() const {return Qs.para;}
vector<AdditionalPert> Sample::get_addQs() const {return Qs.add_para;}
vector<double> Sample::get_Qp() const {return Qp.para;}
vector<AdditionalPert> Sample::get_addQp() const {return Qp.add_para;}
vector<double> Sample::get_T() const {return T.para;}
vector<AdditionalPert> Sample::get_addT() const { return T.add_para;}
void Sample::print_sample() const {
    cout<<"Thk:"<<thk<<" +- "<<thk-left_thk<<endl;
    cout<<"Vs:";vs.print_sample();
    cout<<"Vp/Vs:";vpvs.print_sample();
    cout<<"rho:";rho.print_sample();
    cout<<"Qs:";Qs.print_sample();
    cout<<"Qp:";Qp.print_sample();
    cout<<"T:";T.print_sample();
    cout<<"P:";P.print_sample();
}
void Sample::print() const {
    vs.print();vpvs.print();rho.print();Qs.print();T.print();P.print();
    if(step_thk>0.) cout<<thk;
    cout<<"\n";
}
void Sample::write_to_file(std::ofstream &f) const {
    if(step_thk>0) f<<thk<<" ";
    //vs
    for(int k=0;k<vs.n_para;k++){
        if(vs.step[k]>0) f<<vs.para[k]<<" ";
    }
    for(int k=0;k<vs.n_pert;k++){
        if(vs.add_step[k].get_top()>0) f<<vs.add_para[k].get_top()<<" ";
        if(vs.add_step[k].get_bot()>0) f<<vs.add_para[k].get_bot()<<" ";
        if(vs.add_step[k].get_value()>0) f<<vs.add_para[k].get_value()<<" ";
    }
    //vpvs
    for(int k=0;k<vpvs.n_para;k++){
        if(vpvs.step[k]>0) f<<vpvs.para[k]<<" ";
    }
    for(int k=0;k<vpvs.n_pert;k++){
        if(vpvs.add_step[k].get_top()>0) f<<vpvs.add_para[k].get_top()<<" ";
        if(vpvs.add_step[k].get_bot()>0) f<<vpvs.add_para[k].get_bot()<<" ";
        if(vpvs.add_step[k].get_value()>0) f<<vpvs.add_para[k].get_value()<<" ";
    }
    //rho
    for(int k=0;k<rho.n_para;k++){
        if(rho.step[k]>0) f<<rho.para[k]<<" ";
    }
    for(int k=0;k<rho.n_pert;k++){
        if(rho.add_step[k].get_top()>0) f<<rho.add_para[k].get_top()<<" ";
        if(rho.add_step[k].get_bot()>0) f<<rho.add_para[k].get_bot()<<" ";
        if(rho.add_step[k].get_value()>0) f<<rho.add_para[k].get_value()<<" ";
    }
    //Qs
    for(int k=0;k<Qs.n_para;k++){
        if(Qs.step[k]>0) f<<Qs.para[k]<<" ";
    }
    for(int k=0;k<Qs.n_pert;k++){
        if(Qs.add_step[k].get_top()>0) f<<Qs.add_para[k].get_top()<<" ";
        if(Qs.add_step[k].get_bot()>0) f<<Qs.add_para[k].get_bot()<<" ";
        if(Qs.add_step[k].get_value()>0) f<<Qs.add_para[k].get_value()<<" ";
    }
    //Qp
    for(int k=0;k<Qp.n_para;k++){
        if(Qp.step[k]>0) f<<Qp.para[k]<<" ";
    }
    for(int k=0;k<Qp.n_pert;k++){
        if(Qp.add_step[k].get_top()>0) f<<Qp.add_para[k].get_top()<<" ";
        if(Qp.add_step[k].get_bot()>0) f<<Qp.add_para[k].get_bot()<<" ";
        if(Qp.add_step[k].get_value()>0) f<<Qp.add_para[k].get_value()<<" ";
    }
    //T
    for(int k=0;k<T.n_para;k++){
        if(T.step[k]>0) f<<T.para[k]<<" ";
    }
    for(int k=0;k<T.n_pert;k++){
        if(T.add_step[k].get_top()>0) f<<T.add_para[k].get_top()<<" ";
        if(T.add_step[k].get_bot()>0) f<<T.add_para[k].get_bot()<<" ";
        if(T.add_step[k].get_value()>0) f<<T.add_para[k].get_value()<<" ";
    }
    //P
    for(int k=0;k<P.n_para;k++){
        if(P.step[k]>0) f<<P.para[k]<<" ";
    }
    for(int k=0;k<P.n_pert;k++){
        if(P.add_step[k].get_top()>0) f<<P.add_para[k].get_top()<<" ";
        if(P.add_step[k].get_bot()>0) f<<P.add_para[k].get_bot()<<" ";
        if(P.add_step[k].get_value()>0) f<<P.add_para[k].get_value()<<" ";
    }
    //f<<"\n";
}
//------------------ end of Sample Class --------------------------------

// ------------------ Class -----------------------------
swData::swData(const char &sw_tp, const char &tp) {
    this->sw_tp = sw_tp;
    this->data_tp = tp;
    T.clear();data.clear();err.clear();ref_data.clear();
    n=(int)data.size();
}
swData::~swData() = default;

int swData::read_from_file(const string &fname) {
    ifstream fin;
    int rows, cols;
    double tmp;
    fin.open(fname);
    if(!fin) {cout<<"cannot open the file "<<fname<<endl;exit(0);}
    fin >> rows >> cols;
    T.resize(rows,DEFAULTd);data.resize(rows,DEFAULTd);err.resize(rows,DEFAULTd);ref_data.resize(rows,DEFAULTd);
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            switch (j) {
                case 0:
                    fin >> T[i];break;
                case 1:
                    fin >> data[i];break;
                case 2:
                    fin >> err[i];break;
                case 3:
                    fin >> ref_data[i];break;
                default:
                    cout<<"the #column in sw file is wrong."<<endl;
            }
        }
    }
    n = rows;
    fin.close();
    return 1;
}
char swData::get_tp() const {return sw_tp;}
vector<double> swData::get_data() const {return data;}
vector<double> swData::get_err() const {return err;}
int swData::print() const {
    for(int i=0;i<n;i++){
        cout<<T[i]<<" "<<data[i]<<" "<<err[i]<<" "<<ref_data[i]<<endl;
    }
    return 0;
}
void swData::Clear() {T.clear();data.clear();err.clear();ref_data.clear();n=DEFAULTi;}
void swData::write_to_txt(const std::string &nm) const {
    ofstream f;
    string fnm= nm + "_"+sw_tp+data_tp+".dat";
    f.open(fnm);
    //f<<n<<" 4"<<endl;
    for (int i=0;i<n;i++){
        f<<T[i]<<" "<<data[i]<<" "<<err[i]<<" "<<ref_data[i]<<endl;
    }
    f.close();
}

void SWData::add_sw(const swData &sw) {this->sw.push_back(sw);n = n++;}
int SWData::Disp(const Lyr_model &lyrm) {
    int newnlayer,nn,k,nper,nper1,N=200; //lyrm.get_nlyr()+10;
    float *tvs, *tvp, *tvpvs, *tqs, *tqp, *trho, *tthick; // can use a 2d array for simplicity. but to make code clear, use 1d here.
    float *uR0,*uL0,*cR0,*cL0,*rR0,*rL0,*period,*ampl0;
    float temp_amp;
    vector<double>::iterator id;
    vector<double> period1,tvel;
    tvs = new float[N]; tvp = new float[N]; tvpvs = new float [N];
    tqs= new float[N]; tqp = new float [N]; period = new float [N];
    trho= new float[N]; tthick= new float[N];
    uR0 = new float[N]; uL0 = new float[N]; cR0= new float[N]; cL0= new float[N];
    rR0 = new float[N]; rL0 = new float[N]; ampl0 = new float[N];
    for(int i=0;i<N;i++){
        *(tvs+i)=*(tvp+i)=*(tqs+i)=*(trho+i)=*(tthick+i)=0.;
        *(uR0+i)=*(uL0+i)=*(cR0+i)=*(cL0+i)=*(period+i)=0.;
        *(rR0+i)=*(rL0+i)=0.;
    }
    nn=0;
    newnlayer = lyrm.get_nlyr();
    //cout<<lyrm.get_nlyr()<<endl;
    //lyrm.get_nlyr();
    for(int i=0;i<newnlayer;i++){
        //cout<<nn<<endl;
        tvs[i]=(float)lyrm.get_vs()[i];
        tvp[i]=(float)lyrm.get_vp()[i];
        tvpvs[i] = tvp[i]/tvs[i];
        trho[i]=(float)lyrm.get_rho()[i];
        tthick[i]=(float)lyrm.get_thk()[i];//depth
        tqs[i] = (float)1./lyrm.get_Qs()[i]; // todo to be changed
        tqp[i] = (float)1./lyrm.get_Qp()[i];
        nn=nn+1;
    }//for i
    tvs[nn]=tvs[nn-1]+0.01;
    tvp[nn]=tvp[nn-1]+0.01;
    tqs[nn]=tqs[nn-1];//tqs[nn-2];
    trho[nn]=trho[nn-1]+0.01;
    tthick[nn]=0.;
    nn=nn+1;//#of element inside the arrays
    //cout<<"add half space:done"<<endl;
    /***
    for(int i=0;i<nn;i++){
        cout<<tthick[i]<<" "<<tvs[i]<<" ";
        cout<<tvp[i]<<" "<<trho[i]<<" "<<tqs[i]<<endl;
    }***/
    nper = 0;
    period1.clear();
    for(int i=0;i<n;i++){
        period1.insert(period1.end(),sw[i].T.begin(),sw[i].T.end());
    }
    //cout<<"combine all periods"<<endl;
    sort(period1.begin(),period1.end());
    //cout<<"sort period"<<endl;
    id = unique(period1.begin(),period1.end());
    period1.erase(id,period1.end());
    //cout<<"erase:"<<period1.size()<<endl;
    nper = period1.size();
    /***
    for(int i=0;i<nper;i++){
        cout<<period1[i]<<endl;
    }***/
    //cout<<"new a period,size"<<nper<<endl;
    for(int j=0;j<nper;j++){
        period[j] = (float)period1[j];
    }
    //cout<<"done: period"<<endl;
    int cflag = 2; // todo : add Love wave : 1?
    //cout<<"nn:"<<nn<<endl;
    fast_surf_(&nn,&cflag,tvp,tvs,trho,tthick, tqs,period,&nper,uR0,uL0,cR0,cL0,rR0,rL0,ampl0);
    /***
    for(int i=0;i<nper;i++){
        cout<<period[i]<<" "<<uR0[i]<<" "<<uL0[i]<<" "<<cR0[i]<<" "<<cL0[i]<<" "<<rR0[i]<<" "<<rL0[i]<<" "<<ampl0[i]<<endl;
    }***/
    //cout<<"done: fast_surf"<<endl;
    for(int j=0;j<n;j++){
        switch (sw[j].data_tp) {
            case 'p':
                if (sw[j].T.size()>0){
                    for(int i=0;i<sw[j].T.size();i++){
                        id=find(period1.begin(),period1.end(),sw[j].T[i]);
                        if(id==period1.end()){cout<<"#######!!!canot find pper"<<sw[j].T[i]<<"in perio1\n";exit(0);}
                        if(sw[j].sw_tp == 'R' or sw[j].sw_tp == 'r') sw[j].data[i] = cR0[id - period1.begin()];
                        else if(sw[j].sw_tp == 'L' or sw[j].sw_tp == 'l') sw[j].data[i] = cL0[id - period1.begin()];
                    }
                }break;
            case 'g':
                if(sw[j].T.size()>0){
                    for(int i=0;i<sw[j].T.size();i++){
                        id=find(period1.begin(),period1.end(),sw[j].T[i]);
                        if(id==period1.end()){cout<<"#######!!!canot find pper"<<sw[j].T[i]<<"in perio1\n";exit(0);}
                        if(sw[j].sw_tp == 'R' or sw[j].sw_tp == 'r') sw[j].data[i] = uR0[id - period1.begin()];
                        else if(sw[j].sw_tp == 'L' or sw[j].sw_tp == 'l') sw[j].data[i] = uL0[id - period1.begin()];
                    }
                }break;
            case 'e':
                if(sw[j].T.size()>0){
                    for(int i=0;i<sw[j].T.size();i++){
                        id=find(period1.begin(),period1.end(),sw[j].T[i]);
                        if(id==period1.end()){cout<<"#######!!!canot find pper"<<sw[j].T[i]<<"in perio1\n";exit(0);}
                        if(sw[j].sw_tp == 'R' or sw[j].sw_tp == 'r') sw[j].data[i] = rR0[id - period1.begin()];
                        else if(sw[j].sw_tp == 'L' or sw[j].sw_tp == 'l') sw[j].data[i] = rL0[id - period1.begin()];
                    }
                }break;
            case 'a':
                if(sw[j].T.size()>0){
                    for(int i=0;i<sw[j].T.size();i++){
                        id=find(period1.begin(),period1.end(),sw[j].T[i]);
                        if(id==period1.end()){cout<<"#######!!!canot find pper"<<sw[j].T[i]<<"in perio1\n";exit(0);}
                        if(sw[j].sw_tp == 'R' or sw[j].sw_tp == 'r') sw[j].data[i] = 100. + 100. * (ampl0[id - period1.begin()] - sw[j].ref_data[i]) / sw[j].ref_data[i];
                        //else if(la.tp=='L' or la.tp=='l') la.data[i] = ampl0[id-period1.begin()];
                    }
                }break;
            default: cout<<"wrong data type:"<<sw[j].data_tp<<endl;
        }
    }
    delete[] tvs;
    delete[] tvp;
    delete[] tvpvs;
    delete[] tqs;
    delete[] tqp;
    delete[] period;
    delete[] trho;
    delete[] tthick;
    delete[] uR0;
    delete[] uL0;
    delete[] cR0;
    delete[] cL0;
    delete[] rR0;
    delete[] rL0;
    delete[] ampl0;
}
void SWData::Clear() {sw.clear();}
SWData::SWData() {sw.clear();n=0;}
SWData::~SWData() = default;
void SWData::set_n(const int &n) {this->n = n;}
int SWData::get_n() const {return n;}
vector<double> SWData::get_data(const int &idx) {return sw[idx].data;}
vector<double> SWData::get_err(const int &idx) {return sw[idx].err;}

void SWData::print_swdata() const {
    for(int i=0;i<sw.size();i++){
        cout<<"data type:"<<sw[i].data_tp<<" of "<<sw[i].sw_tp<< "wave"<<endl;
        sw[i].print();
    }
}
void SWData::write_to_txt(const std::string &nm) const {
    for (int i=0;i<sw.size();i++){
        sw[i].write_to_txt(nm);
    }
}
// ------------------ end of swData Class -------------------

// ------------------ RFData Class --------------------------
RFData::RFData() {
    a=DEFAULTd;rayp=DEFAULTd;dist = DEFAULTd; baz = DEFAULTd; snr = DEFAULTd; fitness = DEFAULTd;
    t.clear();rf.clear();err.clear();
    E.clear();Eppps.clear();Eps.clear();Epsps.clear();
    TPsPs.clear();TPs.clear();TPpPs.clear();disco.clear();
}
RFData::RFData(const double &a, const double &rayp) {
    this->a=a;this->rayp=rayp;dist = DEFAULTd; baz = DEFAULTd;
    t.clear();rf.clear();err.clear();
    E.clear();Eppps.clear();Eps.clear();Epsps.clear();
    TPsPs.clear();TPs.clear();TPpPs.clear();disco.clear();
}
RFData::RFData(const vector<double> &time, const vector<double> &amplitude, const vector<double> &err, const double &a, const double &rayp, const vector<int> &disco ){
    n = time.size();
    t = time;
    rf = amplitude;
    this->a = a; this->rayp = rayp;
    this->disco = disco;
    this->err = err;
}
RFData::~RFData() {}

int RFData::read_from_file(const string &fname) {
    ifstream fin;
    int rows, cols;

    fin.open(fname);
    if(!fin) {cout<<"cannot open the file:"<<fname<<endl;exit(0);}
    fin >> rows >> cols;
    t.resize(rows);rf.resize(rows);err.resize(rows);
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            switch (j) {
                case 0:
                    fin >> t[i];break;
                case 1:
                    fin >> rf[i];break;
                case 2:
                    fin >> err[i];break;
                default:
                    cout<<"the #column in RF file is wrong."<<endl;
            }
        }
    }
    n = rows;
    rayp = 0.06;dist = DEFAULTd; baz = DEFAULTd;
    fin.close();
    return 1;
}
int RFData::read_from_sac(const string &fname) {
    SAC::SacStream sac(fname); //todo close file?
    rf = sac.data1;
    //cout<<"nvhdr"<<sac.nvhdr<<endl;
    double b,delta,tmp=0.;
    b = sac.b;
    delta=sac.delta; // todo:: sth wrong
    if(delta<0.00001) delta = sac.user4;
    //if(abs(b)>0.00001) {
    //    cout<<"b time is wrong"<<endl;return 1;
    //}
    n = sac.npts;
    //cout<<"RFData::readfromsac: delta"<<delta<<endl;
    //cout<<"RFData::readfromsac:npts "<<n<<endl;
    if(n != rf.size()) {
        cout << "NPTS is wrong" << endl;
        exit(0);
    };
    t.clear();
    t.push_back(b);
    for (int i=1;i<n;i++){
        tmp = b + i*delta;
        t.push_back(tmp);
        err.push_back(DEFAULTd);
    }
    a = sac.user0; // gaussian parameter
    snr = sac.user2; // SNR ratio
    fitness = sac.user1; //fitness of deconvolution
    rayp = sac.user3/earth_radius; // ray parameter
    // todo: unit of rayp
    dist = sac.dist;
    baz = sac.baz;
    // todo:mag
    //cout<<delta<<" "<<a<<" "<<rayp<<" "<<n<<endl;
    return 1;
}
int RFData::set_disco(const vector<int> &disco) {
    this->disco = disco;
    return 0;
}

int RFData::RF(const Lyr_model &lyrm) { // ????, to be done
    rf.resize(n,0.);
    int N = lyrm.get_nlyr()+1;
    float *tvs, *tvpvs, *tqs, *tqp, *trho, *tthick;
    tvs = new float[N];  tvpvs = new float [N];
    tqs= new float[N]; tqp= new float[N];
    trho= new float[N]; tthick= new float[N];
    for(int i=0;i<lyrm.get_nlyr();i++){
        tvs[i] = (float)lyrm.get_vs()[i];
        tvpvs[i] = (float)lyrm.get_vpvs()[i];
        tqs[i] = (float)lyrm.get_Qs()[i];
        tqp[i] = (float)lyrm.get_Qp()[i];
        trho[i] = (float)lyrm.get_rho()[i];
        tthick[i] = (float)lyrm.get_thk()[i];
    }
    int nl = N-1;
    tvs[nl]=tvs[nl-1];
    tvpvs[nl]=tvpvs[nl-1];
    tqs[nl]=tqs[nl-1];
    tqp[nl]=tqp[nl-1];
    tthick[nl]=0.;
    nl=nl+1; // total number of lyr, including infinite half-space
    float slow,din,rt;
    float *rx;
    rx=new float[n];
    int nn = n; //#points in rf waveform
    slow=rayp;
//  cout<<"===== pi="<<pi<<endl;
    din=180.*asin(tvs[nl-1]*tvpvs[nl-1]*slow)/pi;
//    cout<<tvs[nl-1]*tvpvs[nl-1]<<" "<<slow<<" "<<pi<<endl;
    float gau = a;
    float dt = 0.005; // ???parameter of a minimum amplitude level
    float t0 = 0; //
    //cout<<"RFData::RF: t0"<<t0<<endl;
    int nn1 = nl-1; //???
    //cout<<"RFdata::RF-"<<"slow:"<<slow<<"; gau:"<<gau;
    rt = int(1./(t[1]-t[0])+0.5); // todo: to be disccussed
    //cout<<t[1]<<" "<<t[0]<<endl;
    //cout<<nn1<<" "<<rt<<" "<<din<<" "<<gau<<" "<<dt<<" "<<t0<<" "<<rx<<endl;
    theo_(&nn1,tvs,tthick,tvpvs,tqp,tqs,&rt,&din,&gau,&dt,&t0,&nn,rx);
    for(int i=0;i<nn;i++){
        rf[i] = rx[i];
        t[i] = t0+i*1/rt;
    }
    delete[] tvs;
    delete[] tvpvs;
    delete[] tqs;
    delete[] tqp;
    delete[] trho;
    delete[] tthick;
    delete[] rx;
    return 1;
}
double RFData::get_a() const {return a;}
double RFData::get_rayp() const {return rayp;}
vector<double> RFData::get_t() const {return t;}
vector<double> RFData::get_rf() const {return rf;}
vector<double> RFData::get_err() const {return err;}
vector<double> RFData::get_E() const {return E;}
vector<double> RFData::get_Eppps() const {return Eppps;}
vector<double> RFData::get_Eps() const {return Eps;}
vector<double> RFData::get_Epsps() const {return Epsps;}
vector<double> RFData::get_Tppps() const {return TPpPs;}
vector<double> RFData::get_Tps() const {return TPs;}
vector<double> RFData::get_Tpsps() const {return TPsPs;}
int RFData::cal_Energy(const double &wPs, const double &wPsPs, const double &wPpPs) {
    int idx;
    Eps.clear();Epsps.clear();Eppps.clear();E.clear();
    //cout<<"RFData::calEnergy:weight "<<wPsPs<<endl;
    //cout<<"RFdata:cal_Energy-disco.size"<<disco.size()<<endl;
    for (int i=0;i<disco.size();i++){
        idx = find_closest_value(t, TPs[i]);
        Eps.push_back(rf[idx]);
        idx = find_closest_value(t,TPsPs[i]);
        Epsps.push_back(rf[idx]);
        idx = find_closest_value(t, TPpPs[i]);
        Eppps.push_back(rf[idx]);
        E.push_back(wPs*Eps[i]+wPpPs*Eppps[i]-wPsPs*Epsps[i]);
        //cout<<"RFdata:cal_Energy-"<<i<<" "<<E[i]<<endl;
        //cout<<"RFData::calEnergy:Epsps "<<Epsps[i]<<endl;
    }
    return 1;
}
int RFData::cal_arrival(const Lyr_model &lyrm) {
    vector<double> vs, vp, dep, thk;
    double sqrt_s,sqrt_p;
    TPs.clear();TPpPs.clear();TPsPs.clear();
    TPs.resize(disco.size(),0.);TPpPs.resize(disco.size(),0.);TPsPs.resize(disco.size(),0.);
    vs = lyrm.get_vs();vp = lyrm.get_vp();
    thk = lyrm.get_thk(); dep = lyrm.get_depth();
    //cout<<"RFData::dsco.size"<<disco.size()<<endl;
    //cout<<"RFData::Tps.size "<<TPs.size()<<endl;
    for(int j=0;j<disco.size();j++){
        for(int i=0;i<disco[j];i++){
            sqrt_s = sqrt(1/(vs[i]*vs[i])-rayp*rayp);
            sqrt_p = sqrt(1/(vp[i]*vp[i])-rayp*rayp);
            TPs[j] = TPs[j] + thk[i]*(sqrt_s-sqrt_p);
            TPpPs[j] = TPpPs[j] + thk[i]*(sqrt_s+sqrt_p);
            TPsPs[j] = TPsPs[j] + thk[i]*sqrt_s*2.;
            //cout<<"RFData::cal_arrival:"<<TPs[j]<<" "<<TPpPs[j]<<" "<<TPsPs[j]<<endl;
        }
    }
    return 0;
}
void RFData::print() const {
    for(int i=0;i<n;i++){
        cout<<t[i]<<" "<<rf[i]<<" "<<err[i]<<endl;
    }
}
void RFData::write_to_txt(const string &nm) const {
    ofstream f;
    string fnm= nm + "."+ to_string(a)+"_"+ to_string(rayp)+".rf.dat";
    f.open(fnm);
    //f<<n<<" 3"<<endl;
    for(int i=0;i<n;i++){
        f<<t[i]<<" "<<rf[i]<<" "<<err[i]<<endl;
    }
    f.close();
}
void RFData::write_to_sac(const std::string &nm) const {
    SAC::SacStream sac; //("random_1.sac"); //todo:: change template sac
    sac.data1 = rf;
    sac.b = t[0];
    sac.delta = t[1]-t[0];
    sac.npts = n;
    sac.dist = dist;
    sac.baz = baz;
    sac.user0 = a;
    sac.user1 = fitness;
    sac.user2 = snr;
    sac.user3 = rayp*earth_radius;
    sac.nvhdr = 7;
    filesystem::path fnm= "./"+nm + "."+ to_string(a)+"_"+ to_string(rayp)+".rf.sac";
    sac.write(fnm);
}
// ------------------- end of RF Class ----------------------

// ------------------- HkData Class --------------------------
HkData::HkData() {
    a=DEFAULTd;rayp=DEFAULTd;dist=DEFAULTd;baz=DEFAULTd;
    E.clear();Eps.clear();Eppps.clear();Epsps.clear();
    TPs.clear();TPpPs.clear();TPsPs.clear();
    disco.clear();
}
HkData::~HkData() {}

void HkData::set_idsco(const vector<int> &disco) {
    this->disco=disco;
}
void HkData::set(const double &a, const double &rayp, const double &dist, const double &baz) {
    this->a = a; this->rayp = rayp; this->dist = dist; this->baz = baz;
}
void HkData::set_E(vector<double> E, vector<double> Eps, vector<double> Eppps, vector<double> Epsps) {
    this->E = E; this->Eps = Eps; this->Eppps = Eppps; this->Epsps = Epsps;
}
void HkData::set_T(vector<double> Tps, vector<double> Tppps, vector<double> Tpsps) {
    this->TPs = Tps; this->TPsPs = Tpsps; this->TPpPs = Tppps;
}
vector<double> HkData::get_Eps() const {return Eps;}
vector<double> HkData::get_Epsps() const {return Epsps;}
vector<double> HkData::get_Eppps() const {return Eppps;}
vector<double> HkData::get_E() const {return E;}
double HkData::get_a() const {return a;}
void HkData::write_to_txt(ofstream &f) const {
    for(int i=0;i<TPs.size();i++){
        f<<a<<" "<<rayp<<" "<<TPs[i]<<" "<<TPpPs[i]<<" "<<TPsPs[i]<<" ";
    }
    f<<"\n";
}
// ------------------- end of HkData Class -------------------

// ------------------- Model Class ---------------------------
int Model::acount=0;
Model::Model() {
    groups.clear(); n_group = groups.size();
    f_acpt = DEFAULTi; f_E = DEFAULTi;f_monot = DEFAULTi; idx_monot.clear();
    pre_rf.clear();obs_rf.clear();
    pre_hk.clear();obs_hk.clear();
    pre_sw.Clear();obs_sw.Clear();
    sw_mis_one.clear();rf_mis_one.clear();
    sw_S_one.clear();rf_S_one.clear();
    sw_mis = DEFAULTd; rf_mis = DEFAULTd; Joint_mis = DEFAULTd;sw_L = DEFAULTd; rf_L = DEFAULTd;
    E_PpPs.clear();E_PsPs.clear();E_Ps.clear();E.clear();wE.clear();
    totalE = DEFAULTd;L_mis = DEFAULTd; L_E = DEFAULTd;refE = DEFAULTd;
    a_indx = DEFAULTi;
}
Model::Model(int n_groups) {
    groups.clear();
    for(int i=0;i<n_groups;i++) {
        Group tmp(n_groups);
        groups.push_back(tmp);
    }
    this->n_group = n_groups;
    f_acpt = DEFAULTi; f_E = DEFAULTi; f_monot = DEFAULTi; idx_monot.clear();
    pre_sw.Clear();obs_sw.Clear();pre_rf.clear();obs_rf.clear();
    pre_hk.clear();obs_hk.clear();
    sw_mis_one.clear();rf_mis_one.clear();
    sw_S_one.clear();rf_S_one.clear();
    sw_mis = DEFAULTd; rf_mis = DEFAULTd; Joint_mis = DEFAULTd;sw_L = DEFAULTd; rf_L = DEFAULTd;
    E_PpPs.clear();E_PsPs.clear();E_Ps.clear();E.clear();wE.clear();
    totalE = DEFAULTd;L_mis = DEFAULTd; L_E = DEFAULTd;refE = DEFAULTd;
    a_indx = DEFAULTi;
}
Model::Model(const vector<Group>& groups, const SWData &obs_sw, const vector<RFData>& obs_rf, const vector<vector<RFData> >& obs_hk) {
    this->groups = groups;
    this->obs_sw = obs_sw;
    this->obs_rf = obs_rf;
    this->obs_hk = obs_hk;
    n_group = this->groups.size();
    f_acpt = DEFAULTi; f_E = DEFAULTi; f_monot = DEFAULTi; idx_monot.clear();
    pre_sw.Clear();pre_rf.clear();pre_hk.clear();
    sw_mis_one.clear();sw_S_one.clear();
    rf_mis_one.clear();rf_S_one.clear();
    sw_mis = DEFAULTd; rf_mis = DEFAULTd; Joint_mis = DEFAULTd;sw_L = DEFAULTd;rf_L = DEFAULTd;
    w_sw = DEFAULTd; w_rf = DEFAULTd;
    wPs.clear();wPsPs.clear();wPpPs.clear();
    //wPs = DEFAULTd; wPsPs = DEFAULTd; wPpPs = DEFAULTd;
    E_PpPs.clear();E_PsPs.clear();E_Ps.clear();
    E.clear();wE.clear();
    totalE = DEFAULTd;L_mis = DEFAULTd; L_E = DEFAULTd;refE = DEFAULTd;
    a_indx = DEFAULTi; //??
}
//void Model::reset_count(int num) {count=num;}
//void Model::reset_indx(int idx) {indx = idx;}

Model::~Model() {}
void Model::set_refE(const double &refE) {
    this->refE = refE;
}
void Model::set_weight(const double &w_rf, const vector<double> &wPs, const vector<double> &wPsPs, const vector<double> &wPpPs, const vector<vector<double> > &wE) {
    this->w_rf=w_rf;
    w_sw = 1. - w_rf;
    this->wPs = wPs;
    this->wPpPs = wPpPs;
    this->wPsPs = wPsPs;
    this->wE = wE;
}
void Model::set_monot(const vector<int> &monot) {this->idx_monot = monot;}
void Model::set_group(const vector<Group> &group) {
    groups = group;
    /***
    for(int i=0;i<n_group;i++){
        cout<<"thk:"<<groups[i].get_thk()<<endl;
    }
    ***/
}
int Model::readin_groups(const string &fname) {
    double dep;
    ifstream f;
    f.open(fname);
    if(f.is_open()){
        string line; vector<string> v;
        while(getline(f,line)){
            Split(line, v, " ");
            int idx = atoi(v[0].c_str()); // index of group
            //cout<<"Model::readin_groups-"<<idx<<endl;
            groups[idx].readin(v);
        }
    }
    else {
        cout<<"couldn't open group file:"<<fname<<endl;
        return 0;
    }
    f.close();
    cout<<"Model::readin_groups:f.close"<<endl;
    Thk =0.0;
    for(int i=0;i<n_group;i++){
        if(i==0) {groups[i].set_dep();}
        else {
            dep = groups[i-1].get_dep() + groups[i-1].get_thk();
            groups[i].set_dep(dep);
        }
        Thk = Thk + groups[i].get_thk();
    }
    cout<<"done: read in groups file."<<endl;
    return 1;
}
int Model::readin_swdata(const vector<string> &fname, const vector<char> &tp, const vector<char> &sw_tp) {
    obs_sw.set_n(tp.size());
    for(int i=0;i<sw_tp.size();i++) {
        swData tmp(sw_tp[i], tp[i]);
        tmp.read_from_file(fname[i]);
        obs_sw.add_sw(tmp);
    }
   // cout<<"done: read in swdata"<<endl;
    return 1;
}
int Model::readin_rfdata(const vector<string> &fname, const vector<double> &a, const vector<double> &rayp) {
	obs_rf.clear();
    for(int i=0;i<a.size();i++){
        RFData tmp(a[i],rayp[i]); //???? rayp ?
        tmp.read_from_file(fname[i]);
        obs_rf.push_back(tmp);
    }
    //cout<<"done: read in rf data."<<endl;
    return 1;
}
int Model::readin_hkdata(const vector<vector<std::string> > &fname, const vector<vector<int>> &idisco) {
    obs_hk.clear();pre_hk.clear();
    for(int i=0;i<fname.size();i++){
        vector<RFData> ttmp;
        vector<HkData> hk_ttmp;
        hk_ttmp.clear();
        ttmp.clear();
        vector<int> tmp_disco; tmp_disco.resize(idisco[i].size(),0);
        for(int j=0;j<idisco[i].size();j++){
            for(int k=0;k<=idisco[i][j];k++){
                tmp_disco[j] = groups[k].get_nlyr()+tmp_disco[j];
            }
            cout<<"tmp_disco"<<j<<" "<<tmp_disco[j]<<endl;
        }
        for(int j=0;j<fname[i].size();j++) {
            RFData tmp;
            tmp.read_from_sac(fname[i][j]);
            tmp.set_disco(tmp_disco);
            ttmp.push_back(tmp);
            HkData hk_tmp;
            hk_tmp.set_idsco(tmp_disco);
	    hk_tmp.set(tmp.get_a(),tmp.get_rayp(),DEFAULTd,DEFAULTd);
            //cout<<"Model:readin_hkdata:idisco.size"<<idisco[i].size()<<endl;
            hk_ttmp.push_back(hk_tmp);
        }
        obs_hk.push_back(ttmp);
        pre_hk.push_back(hk_ttmp);
    }
    cout<<"done: read in hk data"<<endl;
    return 1;
}

bool Model::check_prior() {
    bool flag = true;
    // check monotonically increase (Vs)
    for(int i=0;i<idx_monot.size();i++){
        //cout<<"checking mono"<<endl;
        flag = check_monot(groups[i].get_lyr_vs());
        if(!flag) {
            cout<<"mono check failed"<<endl;
            return flag;
        }
    }
    // check Vs < 4.9 km/s
    //cout<<"checking max Vs"<<lyrm.get_vs().size()<<endl;
    vector<double> tmp = lyrm.get_vs();
    double max_Vs = *max_element(tmp.begin(),tmp.end());
    //cout<<"max_vs "<<max_Vs;
    flag = (max_Vs < MAX_VS); if(!flag) {
        //cout<<"maxVs failed:"<<max_Vs<<endl;
        return flag;
    }
    //cout<<"check maxVs:"<<flag<<endl;
    // check Vs jump
    for(int i=0;i<n_group-1;i++){
        //cout<<"checking Vs jump"<<endl;
        flag = (groups[i].get_lyr_vs().back()<groups[i+1].get_lyr_vs().front());
        if(!flag) {
            //cout<<"Vs jump failed:"<<groups[i].get_lyr_vs().back()<<" "<<groups[i+1].get_lyr_vs().front()<<endl;
            return flag;
        }
    }
    // todo: velocity gradient
    //cout<<"pass prior check"<<endl;
    return flag;
}
int Model::update_groups(const vector<Sample> &s) {
    double dep;
    double tmp=0.;
    for (int i=0;i<n_group-1;i++){
        groups[i].update_from_sample(s[i]);
        if(i>0) {
            dep = groups[i-1].get_dep() + groups[i-1].get_thk();
            groups[i].set_dep(dep);
        }
        tmp = tmp + groups[i].get_thk();
        //cout<<"Model::updataGroups:thk "<<groups[i].get_thk()<<endl;
    }
    groups[n_group-1].update_from_sample(s[n_group-1],Thk-tmp);
    //cout<<"Model::updateGroups:thk"<<groups[n_group-1].get_thk()<<endl;
    dep = groups[n_group-2].get_dep() + groups[n_group-2].get_thk();
    groups[n_group-1].set_dep(dep);
    //cout<<"done: update groups"<<endl;
    return 1;
}
int Model::update_lyrmod() { // todo add T&P
    int nlyr=0;
    vector<double> vs,vp,thk,rho,vpvs,Qs,Qp,dep;
    vs.clear();vp.clear();vpvs.clear();thk.clear();rho.clear();Qs.clear();Qp.clear();dep.clear();
	//cout<<n_group<<endl;
    for(int i=0;i<n_group;i++){
        groups[i].update_lyrmodel();
        //cout<<"update lyr:"<<i<<"th group thk "<<groups[i].get_thk()<<endl;
        vs = joint_lyr(vs,groups[i].get_lyr_vs());
        vp = joint_lyr(vp,groups[i].get_lyr_vp());
        vpvs = joint_lyr(vpvs,groups[i].get_lyr_vpvs());
        rho = joint_lyr(rho,groups[i].get_lyr_rho());
        thk = joint_lyr(thk,groups[i].get_lyr_thk());
        dep = joint_lyr(dep, groups[i].get_lyr_dep());
        Qs = joint_lyr(Qs, groups[i].get_lyr_Qs());
        Qp = joint_lyr(Qp, groups[i].get_lyr_Qp());
        nlyr = nlyr+groups[i].get_nlyr();
		//cout<<"Model::update_lyr:nlyr=="<<nlyr<<endl;
    }
    lyrm.set_vs(vs);lyrm.set_vp(vp);lyrm.set_vpvs(vpvs);lyrm.set_rho(rho);
    lyrm.set_Q(Qs,Qp);lyrm.set_thkdep(thk,dep);//lyrm.update_info();
    lyrm.set_nlyr(nlyr);
    //cout<<"done: update lyr model"<<endl;
    return 1;
}
int Model::update_SWdata() {
    SWData tmp(obs_sw);
    tmp.Disp(lyrm);
    pre_sw = tmp;
    //cout<<"done: update swdata"<<endl;
    return 0;
}
int Model::update_RFdata() {
    pre_rf.clear();
    for (auto tmp : obs_rf){ //????
        tmp.RF(lyrm);
        pre_rf.push_back(tmp);
    }
    //cout<<"done: update rfdata"<<endl;
    return 0;
}
int Model::update_Hkdata(int flag) { // flag=1, generate RFwaveform; flag = 0, generate only arrival time
    f_E = 1;
    switch (flag) {
        case 0: // only update travel time & Energy
            //cout<<"Model::updata_Hkdata:pre_hk.size"<<pre_hk.size()<<endl;
		//lyrm.print_lyr();
            for (int i = 0; i < pre_hk.size(); i++) {
                for (int j = 0; j < pre_hk[i].size(); j++) {
                    obs_hk[i][j].cal_arrival(lyrm);
                    //lyrm.print_lyr();
                    vector<double> tps,tpsps, tppps;
                    tps = obs_hk[i][j].get_Tps();
                    tpsps = obs_hk[i][j].get_Tpsps();
                    tppps = obs_hk[i][j].get_Tppps();
                    //cout<<"Model::update_Hkdata:arrival time:"<<tps[0]<<" "<<tppps[0]<<" "<<tpsps[0]<<endl;
                    pre_hk[i][j].set_T(tps,tppps,tpsps);
                    obs_hk[i][j].cal_Energy(wPs[i],wPsPs[i],wPpPs[i]);
                    vector<double> E,Eps,Eppps,Epsps;
                    E = obs_hk[i][j].get_E();
                    Eps = obs_hk[i][j].get_Eps();
                    Epsps = obs_hk[i][j].get_Epsps();
                    Eppps = obs_hk[i][j].get_Eppps();
                    for(int k=0;k<Eps.size();k++){
                        if(Eps[k]<0 || Eppps[k]<0 || Epsps[k]>0) f_E = 0;
                    }
                    //cout<<"update_Hkdata0-Esize:"<<E.size()<<endl;
                    //cout<<"Model::update_Hkdata-E.size"<<E.size()<<endl;
                    pre_hk[i][j].set_E(E,Eps,Eppps,Epsps);
                }
            }
            break;
        case 1: //update waveform
            for (int i = 0; i < obs_hk.size(); i++) {
                for (int j = 0; j < obs_hk[i].size(); j++) {
                    obs_hk[i][j].RF(lyrm);
                    obs_hk[i][j].cal_arrival(lyrm);
                }
            }
            break;
        default:
            cout<<"update Hk data: wrong flag: "<<flag<<endl;
            return 0;
    }
    //cout<<"done: updata hk data"<<endl;
    return 1;
}  // ???
int Model::update_data(int flag){
    if(w_sw>0.0001){
        //cout<<"begin update SWdata"<<endl;
        update_SWdata();
		//cout<<"update SWdata"<<endl;
    } // todo: need to be fixed
    if(w_rf>0.0001) {
        //cout<<"begin update RFdata"<<endl;
        update_RFdata();
		//cout<<" update RFdata"<<endl;
    }
    //cout<<"begin update Hkdata"<<endl;
    update_Hkdata(flag);
	//cout<<"update Hkdata"<<endl;
    return 0;
}

int Model::cal_SWmisfit() {
    sw_mis = 0.0;
    sw_L = 1.0;
    int N=0;
    double mis,S;
    vector<double> tmp(obs_sw.get_n(),0.);
    sw_mis_one.clear();
    sw_S_one.clear();
    for (int i=0;i<obs_sw.get_n();i++){
        vector<double> pre = pre_sw.get_data(i);
        vector<double> obs = obs_sw.get_data(i);
        vector<double> err = obs_sw.get_err(i);
        tmp[i] = Gm_d(pre,obs,err);
        mis = sqrt(tmp[i]/pre.size()); // sqrt((Gm-d)/N)
        N = N + (int)pre.size(); // total number of all sw data points
        sw_mis_one.push_back(mis);
        sw_S_one.push_back(tmp[i]);
        S = reduceLbyhand(tmp[i]);
        sw_L = sw_L*exp(-factor_M*S); // todo: to be discuss
    }
    sw_mis = sqrt(accumulate(tmp.begin(), tmp.end(),0.0)/N);
    //mis = reduceLbyhand(sw_mis);
    //sw_L = exp(-factor_M*mis);

    return 0;
}
int Model::cal_RFmisfit() {
    rf_mis = 0.;rf_L = 1.;
    int N = 0;
    vector<double> tmp(obs_rf.size());
    double mis,S;
    rf_mis_one.clear();rf_S_one.clear();
    for (int i=0;i<obs_rf.size();i++){
        vector<double> pre = pre_rf[i].get_rf();
        vector<double> obs = obs_rf[i].get_rf();
        vector<double> err = obs_rf[i].get_err();
        N = N + pre.size();
        tmp[i] = Gm_d(pre,obs,err);
        mis = sqrt(tmp[i]/pre.size());
        rf_mis_one.push_back(mis);
        rf_S_one.push_back(tmp[i]);
        //rf_mis = rf_mis + rf_mis_one[i]*w_rf_one[i]; //??
        S = reduceLbyhand(tmp[i]);
        rf_L = rf_L*exp(-factor_M*S);
    }
    rf_mis = sqrt(accumulate(tmp.begin(), tmp.end(),0.0)/N);
    return 0;
}
int Model::cal_Jmisfit() {
    double tmp_rf,tmp_sw;
    double S;
    if(w_rf>0.0001 and w_sw>0.0001) {
        tmp_rf = accumulate(rf_S_one.begin(),rf_S_one.end(),0.0);
        tmp_rf = reduceLbyhand(tmp_rf);
        tmp_sw = accumulate(sw_S_one.begin(),sw_S_one.end(),0.0);
        tmp_sw = reduceLbyhand(tmp_sw);
        //cout<<"S_rf "<<tmp_rf;
        //cout<<" S_sw "<<tmp_sw<<endl;
        S = w_rf * tmp_rf + w_sw * tmp_sw;
        Joint_mis = w_rf*rf_mis + w_sw*sw_mis;
        //cout<<"Joint mis: "<<Joint_mis;
        //cout<<" = "<<w_sw<<"*"<<sw_mis<<"+"<<w_rf<<"*"<<rf_mis<<endl;
    }
    else if(w_rf<=0.0001){
        tmp_sw = accumulate(sw_S_one.begin(),sw_S_one.end(),0.0);
        tmp_sw = reduceLbyhand(tmp_sw);
        S =  tmp_sw;
        Joint_mis = sw_mis;
        //cout<<"S: "<<S<<endl;
    }
    else if(w_sw<=0.0001){
        tmp_rf = accumulate(rf_S_one.begin(),rf_S_one.end(),0.0);
        //cout<<"before,after "<<tmp_rf;
        tmp_rf = reduceLbyhand(tmp_rf);
        S = tmp_rf;
        Joint_mis = rf_mis;
        //cout<<" "<<S<<endl;
        //cout<<"Model::cal_Jmisfit:"<<Joint_mis<<endl;
    }
    S = reduceLbyhand(S);
    L_mis = exp(-factor_M*S);
    //cout<<"Model::cal_Jmisfit: "<<L_mis<<" "<<S<<" "<<Joint_mis<<endl;
    return 0;
}
int Model::cal_Energy(int ndisco=1) {
    if(*max_element(wE[0].begin(),wE[0].end())<0){
        totalE = DEFAULTd;
        L_E = 1;
        return 0;
    }
    totalE = 0.;
    E_Ps.clear();E_PsPs.clear();E_PpPs.clear();E.clear();
    vector<double> tmp_ps,tmp_ppps, tmp_psps,tmp;
    //cout<<"cal_Energy-prehk.size:"<<pre_hk.size()<<endl;
    //cout<<"cal_Energy-getE.size:"<<pre_hk[0][0].get_E().size()<<endl;
    for (auto & i : pre_hk){
        tmp_ps.clear();tmp_ppps.clear();tmp_psps.clear();
        //todo: problem with the size
        tmp_ps.resize(i[0].get_Eps().size(),0.0);
        tmp_ppps.resize(i[0].get_Eppps().size(),0.0);
        tmp_psps.resize(i[0].get_Epsps().size(),0.0);
        tmp.resize(i[0].get_E().size(),0.0);
        //cout<<"Model:cal_Energy-E.size"<<i[0].get_E().size()<<endl;
        for (const auto & j : i){
            tmp_ps = tmp_ps + j.get_Eps();
            tmp_ppps = tmp_ppps + j.get_Eppps();
            tmp_psps = tmp_psps + j.get_Epsps();
            tmp = tmp + j.get_E();
        }
        //cout<<"cal_Energy-E.size:"<<i.size()<<endl;
        //cout<<"Model::cal_Energy-tmp.size:"<<tmp.size()<<endl;
        E_Ps.push_back(tmp_ps/i.size());
        E_PpPs.push_back(tmp_ppps/i.size());
        E_PsPs.push_back(tmp_psps/i.size());
        E.push_back(tmp/i.size());
        //cout<<"Model::E.size"<<E.size()<<endl;

    }
    //cout<<"Model:cal_Energy:"<<wE.size()<<" "<<wE[0].size()<<endl;
    //cout<<"Model:cal_Energy:"<<E.size()<<" "<<E[0].size()<<endl;
    for (int i=0;i<E.size();i++){
        for (int j=0;j<E[i].size();j++){
            //cout<<"Model:cal_Energy-"<<"i "<<i<<";j "<<j<<endl;
            totalE = totalE + wE[i][j]*E[i][j];
			//cout<<wE[i][j]<<" "<<E[i][j]<<endl;
            //cout<<"Model:cal_Energy-"<<"i "<<i<<";j "<<j<<endl;
        }
    }
    //cout<<"Model::cal_Energy-totoalE"<<totalE<<endl;
    if(totalE>0) {
        L_E = pow(exp(totalE/refE), factor_E);
        //cout<<"Model::cal_Energy+:total E"<<totalE<<endl;
    }
    else { //?????
        L_E = pow(exp(totalE), factor_E);
        //cout<<"Model::cal_Energy-: total E"<<totalE<<endl;
    }
    //cout<<"Model::cal_Energy:"<<totalE<<" "<<refE<<" "<<L_E<<endl;
    return 0;
}
int Model::cal_mis() {
    if(w_sw>0.0001) {
        cal_SWmisfit();
        //cout << "done: sw misfit" << endl;
    }
    if(w_rf>0.0001) {
        cal_RFmisfit();
        //cout<<"done:rfmisfit"<<endl;
    }
    cal_Jmisfit();
    //cout<<"Model:cal_mis-Lmis "<<L_mis<<endl;
    return 0;
}

int Model::get_n_sw() const {return obs_sw.get_n();}
int Model::get_aidx() const {return a_indx;}
int Model::get_idx() const {return acount;}
int Model::get_ngroup() const {return n_group;}
double Model::get_Lmis() const {return L_mis;}
double Model::get_LE() const {return L_E;}
double Model::get_Thk() const {return Thk;}
vector<double> Model::get_sw_one() const {return sw_mis_one;}
double Model::get_Dispmis() const {return sw_mis;}
vector<double> Model::get_rf_one() const {return rf_mis_one;}
double Model::get_RFmis() const {return rf_mis;}
double Model::get_Jointmis() const {return Joint_mis;}
double Model::get_Energy() const {return totalE;}
vector<vector<double>> Model::get_Eps() const {return E_Ps;}
vector<vector<double>> Model::get_Epsps() const {return E_PsPs;}
vector<vector<double>> Model::get_Eppps() const {return E_PpPs;}
vector<double> Model::get_Lyr_vs() const {return lyrm.get_vs();}
Lyr_model Model::get_lyrm() const {return lyrm;}
int Model::get_flag() const {return f_acpt;}
int Model::get_flagE() const {return f_E;}
vector<Group> Model::get_groups() const {return groups;}
double Model::get_refE() const {return refE;}
int Model::update_accept_flag(const int &flag) {
    f_acpt = flag;
    if(f_acpt==1) {a_indx=acount;acount++;}
    return 0;
}


void Model::print_swdata() const {
    cout<<"Observed SWData:"<<endl;
    obs_sw.print_swdata();
    cout<<"Predicted SWData:"<<endl;
    pre_sw.print_swdata();
    //vector<double> T, pre, obs, err, ref;
}
void Model::print_rfdata() const {
    //cout<<"Obserbed SWData:"<<endl;
    vector<double> t,pre,obs,err;
    for(int i=0;i<obs_rf.size();i++) {
        t = obs_rf[i].get_t();
        pre = pre_rf[i].get_rf();
        obs = obs_rf[i].get_rf();
        err = obs_rf[i].get_err();
        for(int j=0;j<t.size();j++){
            cout<<t[j]<<" "<<pre[j]<<" "<<obs[j]<<" "<<err[j]<<endl;
        }
    }

}
void Model::write_lyr_to_file(const string &dir,const string &fnm) const {
    ofstream f;
    string file="./"+dir+"/"+"MC."+fnm+".lyr";
    f.open(file);
    lyrm.write_to_txt(f);
    f.close();
}
void Model::write_groups_to_file(const string &dir,const string &fnm) const {
    ofstream f;
    string file="./"+dir+"/MC."+fnm+".group";
    f.open(file);
    for(int i=0;i<n_group;i++) {
        //cout<<"Model::writeGroupToFile:"<<groups[i].get_thk()<<endl;
        groups[i].write_to_file(f);
    }
    f.close();
}
void Model::write_data_to_file(const string &nm, const int flag) const {
    pre_sw.write_to_txt(nm);
    for(int i=0;i<pre_rf.size();i++){
        pre_rf[i].write_to_txt(nm);
    }
    if(flag==1){
    	for(int i=0;i<obs_hk.size();i++){
        	for(int j=0;j<obs_hk[i].size();j++)
            	obs_hk[i][j].write_to_sac(nm);
    	}
        for(int i=0;i<obs_hk.size();i++){
            for(int j=0;j<obs_hk[i].size();j++){
                obs_hk[i][j].write_to_txt(nm);
                //obs_hk[i][j].print();
            }
        }
    }
    // write travel time
    for(int i=0;i<pre_hk.size();i++) {
        ofstream f;
        string fnm=nm+ ".hk."+to_string(pre_hk[i][0].get_a())+".hk.T";
        f.open(fnm);
        for(int j=0;j<pre_hk[i].size();j++){
            pre_hk[i][j].write_to_txt(f);
        }
        f.close();
    }
}
void Model::print_groups() const {
    for(int i=0;i<n_group;i++){
        groups[i].print();
    }
}
void Model::write_to_binary(std::ofstream &f) const {
    double tmp;
    int tp;
    //cout<<"n_group:"<<n_group<<endl;
    for(int i=0;i<n_group;i++) {
        groups[i].write_to_bin(f);
       //cout<<"Model::writeBin:group"<<i<<" "<<sizeof(groups[i])<<endl;
    }
    for(int i=0;i<sw_mis_one.size();i++){
        f.write((char *) &sw_mis_one[i],SizeOfDouble);
        //f.write((char *) &sw_S_one[i],SizeOfDouble);
    }
    f.write((char*) &sw_mis,SizeOfDouble);
    f.write((char *) &rf_mis,SizeOfDouble);
    f.write((char *) &Joint_mis, SizeOfDouble);
/***
    for(int i=0;i<E_Ps.size();i++){
        for(int j=0;j<E_Ps[0].size();j++){
            f.write((char *) &E_Ps[i][j],SizeOfDouble);
            f.write((char *) &E_PpPs[i][j],SizeOfDouble);
            f.write((char *) &E_PsPs[i][j],SizeOfDouble);
        }
    }
***/
    f.write((char *) &totalE,SizeOfDouble);
    //f.write((char *) &L_E,SizeOfDouble);
}
void Model::write_info(const std::string &dir, const std::string &fnm) const {
    ofstream f;
    string file = "./" + dir + "/MC."+fnm+".info";
    f.open(file);
    f<<"#Total Energy: "<< totalE<<endl;
	for(int i =0;i<E.size();i++){
		f<<"--Energy for "<<i<<"th hk set: "<<endl;
		for(int j=0;j<E[i].size();j++){
			f<<"----"<<j<<"th disco:"<<E[i][j]<<" "<<E_Ps[i][j]<<" "<<E_PpPs[i][j]<<" "<<E_PsPs[i][j]<<endl;
		}
	}
    f << "#Joint misfit: " << Joint_mis << endl;
    f<<"--SW   misfit: "<<sw_mis<<endl;
    for(int i=0;i<sw_mis_one.size();i++){
        f<<"----"<<i<<"th SW misfit:"<<sw_mis_one[i]<<endl;
    }
    f<<"--RF   misfit: "<<rf_mis<<endl;
    for(int i=0;i<rf_mis_one.size();i++){
        f<<"----"<<i<<"th RF misfit:"<<rf_mis_one[i]<<endl;
    }
}
// ------------------ end of Model Class ----------------------------


// ------------------- MC Class --------------------------------------
MC::MC(int N){
    sample0.clear();sample1.clear();
    L0.resize(N_L,-DEFAULTd);L1.resize(N_L,-DEFAULTd);
    p.resize(N_L,DEFAULTd);
    //models.resize(N+1);
    this->N = N;
    min_disp = 9999.;min_rf = 999.;min_joint = 999.;
    max_E = -999.;
    a_idx.clear();
}

MC::MC(vector<Sample> s, Model m, int n) {
    //models.clear();
    // L0.clear();L1.clear();p.clear();
    //models.reserve(n+1);
    N = n;
    sample0 = s;
    refm = m;
    //models.push_back(m);
}
MC::~MC(){cout<<"an MC machine is destructed."<<endl;}
int MC::mcount = 0;
void MC::set_weight(const double &w_rf) {
    this->w_rf = w_rf;
    w_sw = 1.0 - w_rf;
}
int MC::Initial_m0(const int &n_group,const string &fgroupnm, const vector<string> &fdispnm, const vector<char> &disp_tp, const vector<char> &sw_tp, const vector<string> &frfnm, const vector<double> &a, const vector<double> &rayp, const vector<vector<string> > &fhknm, const vector<vector<int>> &idisco, const double &w_rf,const vector<double> &wPs, const vector<double> &wPsPs,const vector<double> &wPpPs,const vector<vector<double> > &wE, const double &refE, const vector<int> &idx_monot){
    mcount = 0;
    Model m0(n_group);
    cout<<"MC::initial_m0-construct m0"<<endl;
//    m0.reset_indx(0);
//    m0.reset_count(1);
    m0.readin_groups(fgroupnm);
    cout<<"MC_initialm0:read in groups"<<endl;
    m0.update_lyrmod();
    cout<<"MC_initialm0:update lyr"<<endl;
    //read in data
    m0.readin_swdata(fdispnm,disp_tp,sw_tp);
    m0.readin_rfdata(frfnm,a,rayp);
    m0.readin_hkdata(fhknm,idisco);
    cout<<"MC_initialm0:read in data"<<endl;
    m0.set_weight(w_rf,wPs,wPsPs,wPpPs,wE);
    m0.set_refE(refE);
    m0.set_monot(idx_monot);
    //m0.update_data();
    //models.push_back(m0);
    this->m0 = m0;
    refm = m0;
    //models[mcount] = m0;
    mcount++;
    return 1;
}
int MC::Initial_s0(const string &fnm, const int &n_group) {
    sample0.clear();
    sample0.resize(n_group);
    ifstream f;
    f.open(fnm);
    if(f.is_open()){
        string line;vector<string> v;
        while(getline(f,line)){
            Split(line,v," ");
            int idx = atoi(v[0].c_str());
            sample0[idx].Initial(v);
        }
    }
    else{
        cout<<"sth wrong with initial s0"<<endl;
        return 0;
    }
    f.close();
    return 1;
}
int MC::read_s0(const string &fnm) {
    ifstream f;
    f.open(fnm);
    if(f.is_open()){
        string line;vector<string> v;
        while(getline(f,line)){
            Split(line,v," ");
            int idx = atoi(v[0].c_str());
            sample0[idx].readin(v);
        }
    }
    else{
        cout<<"sth wrong with read_s0"<<endl;
        return 0;
    }
    f.close();
    return 1;
}

int MC::rand_starting_point() {
    int nbad=-1;
    double Thk=0.;
    Model tmp(refm);
    //cout<<"#model in MC"<<models.size()<<endl;
    do {
        nbad++;
        if (nbad > Nbad_threshold) {
            cout << "couldn't find a model satisfying prior constrains" << endl;
            return 0;
        }
        Thk =0.;
        //generate a random sample point
        for(int i = 0; i < sample0.size(); i++) {
            sample0[i].rand_start_sample();
            if(i==sample0.size()-1) sample0[i].set_thk(refm.get_Thk()-Thk);
            Thk = Thk + sample0[i].get_thk();
            //sample0[i].print();
        }
        //transfer the sample to a model
        tmp.update_groups(sample0); // update groups from sample
        tmp.update_lyrmod(); // update lyred model from groups
        //tmp.write_lyr_to_file("./","start_"+to_string(nbad)+".lyr"); //todo: to be removed
    }while(!tmp.check_prior());
	//tmp.write_lyr_to_file("./","start.lyr");
    cout<<"rand sample"<<endl;
    //generate data based on this model
    tmp.update_data(0);
    cout<<"random_starting: update data"<<endl;
    tmp.cal_mis(); //??
    //cout<<"done: "<<mcount<<"calculate Lmis:"<<tmp.get_Lmis()<<endl;
    tmp.cal_Energy(1); // need to be changed
    //cout<<"done:"<<mcount<<" calculate energy:"<<tmp.get_LE()<<endl;
    tmp.update_accept_flag(1); //??, accept the first model
    //cout<<"done:"<<mcount<<"update accept_flag"<<endl;
    //models.push_back(tmp);
    //models[mcount] = tmp;
    mcount++;
    m0 = tmp;
    m1 = tmp;
    //cout<<"MC::random_starting-misfit&eneryg "<<m0.get_Jointmis()<<" "<<m0.get_Energy()<<endl;
    set_L0();
    //cout<<"MC::random_starting-L0 "<<L0[0]<<" "<<L0[1]<<endl;
    //min_disp = m0.get_Dispmis(); idx_minDisp = 0;
    //min_rf = m0.get_RFmis(); idx_minRF = 0;
    //min_joint = m0.get_Jointmis(); idx_minJoint = 0;
    //max_E = m0.get_Energy(); idx_maxE = 0;
    sample1 = sample0;
    return 1;
}
void MC::set_L0() {
    //cout<<"MC::setL0-"<<m0.get_Lmis()<<" "<<m0.get_LE()<<endl;
    L0[0]=m0.get_Lmis();L0[1]=m0.get_LE();
    //cout<<"MC::setL0-"<<L0[0]<<" "<<L0[1]<<endl;
}
int MC::gen_good_sample() {
    int nbad=-1;
    double Thk=0.;
    Model tmp(refm);
    do {
        nbad++;
        if (nbad > Nbad_threshold) {
            cout << "couldn't find a model satisfying prior constrains" << endl;
            return 0;
        }
        //generate a random sample point
        sample1 = sample0;
        Thk = 0.;
        for (int i = 0; i < sample1.size(); i++) {
            sample1[i].gen_new_sample();
            if(i==sample1.size()-1) sample1[i].set_thk(refm.get_Thk()-Thk);
            Thk = Thk + sample1[i].get_thk();
            //cout<<"MC::genGoodSample: sample"<<i<<":thk"<<sample1[i].get_thk()<<endl;
        }
        //transfer the sample to a model
        tmp.update_groups(sample1); // update groups from sample
        //cout<<"genGood:n_group"<<tmp.get_ngroup()<<endl;
        //cout<<"MC_gen_goodsample: update groups"<<endl;
        tmp.update_lyrmod(); // update lyred model from groups
        //cout<<"MC_gen_goodsample: update lyr"<<endl;
        //tmp.write_lyr_to_file("./", to_string(nbad)+ to_string(tmp.get_aidx())); // todo: need to be removed
    }while(!tmp.check_prior());
	m1 = tmp;
	return 1;
}
int MC::update_m1(){
    //tmp.write_lyr_to_file("./", to_string(tmp.get_idx()));
    //cout<<"MC:good sample: write to file"<<endl;
    m1.update_data(0);
    //cout<<"MC_gen_goodsample: update data"<<endl;
    m1.cal_mis();
    //cout<<"MC_gen_goodsample: update misfit"<<tmp.get_RFmis()<<endl;
    m1.cal_Energy(1); // todo: need to be changed, not necessarilly calculate energy.
    //cout<<"MC_gen_goodsample: update energy."<<endl;
    //m1 = tmp;
    set_L1();
    return 1;
}
void MC::set_L1() {L1[0]=m1.get_Lmis();L1[1]=m1.get_LE();}
void MC::cal_p() {
    for (int i=0;i<N_L;i++){
        p[i] = L1[i]/L0[i];
		//cout<<L1[i]<<" "<<L0[i]<<" ";
    }
	//cout<<"\n";
}
int MC::accept_model(ofstream &f, vector<Model> &m) {
    random_device generator;
    uniform_real_distribution<double> uni_rand01(0.0,1.0);
    vector<double> randp(N_L);
    vector<int> flag(N_L,0);
    int Flag=1;
    for (int i=0;i<N_L;i++){
        randp[i] = uni_rand01(generator);
        uni_rand01.reset();
        if (p[i]>=randp[i]) flag[i]=1;
        Flag = Flag*flag[i];
    }
    //cout<<"Flag:"<<Flag<<endl;
    if(Flag==1) {
        m1.update_accept_flag(1);
        //cout<<"flag1:m1:n_group"<<m1.get_ngroup()<<endl;
        //if(m1.get_flagE()==1) {
		if(true){
            m1.write_to_binary(f);
            //cout<<"write to bin"<<endl;
        }
        //a_idx.push_back(m1.get_aidx());
        a_idx.push_back(mcount);
        //cout<<"push back mcount"<<endl;
        sample0 = sample1;
        L0 = L1;
        m0 = m1;
        //m.push_back(m1); //todo: to be done
        //cout<<"model accepted"<<endl;
    }
    //models[mcount]=m1;
    mcount++;
    //cout<<m1.get_Dispmis()<<" "<<m1.get_RFmis()<<endl;
    return 0;
}

double MC::get_minDisp() const {return min_disp;}
int MC::get_idxminDisp() const {return idx_minDisp;}
double MC::get_minRF() const {return min_rf;}
int MC::get_idxminRF() const {return idx_minRF;}
double MC::get_minJoint() const {return min_joint;}
int MC::get_idxminJoint() const {return idx_minJoint;}
double MC::get_maxE() const {return max_E;}
int MC::get_idxmaxE() const {return idx_maxE;}
//vector<Model> MC::get_models() const {return models;}
vector<int> MC::get_adix() const {return a_idx;}

void MC::write_samples_to_file(int idx_thread, int idx_search, int idx_model, ofstream &f) const {
    f<<idx_thread<<" "<<idx_search<<" "<<idx_model<<" "<<m1.get_flag()<<" "<<m1.get_flagE()<<" ";
    for(const auto & i : sample1) i.write_to_file(f);
    f<<m1.get_Lmis()<<" ";
    f<<m1.get_LE()<<" ";
    for(int i=0;i<m1.get_sw_one().size();i++){
        f<<m1.get_sw_one()[i]<<" ";
    }
    f<<m1.get_Dispmis()<<" ";
    for(int i=0;i<m1.get_rf_one().size();i++){
        f<<m1.get_rf_one()[i]<<" ";
    }
    f<<m1.get_RFmis()<<" ";
    f<<m1.get_Jointmis()<<" "<<m1.get_Energy()<<" ";
    for(int i=0;i<m1.get_Eps().size();i++){
        for(int j=0;j<m1.get_Eps()[i].size();j++){
            f<<m1.get_Eps()[i][j]<<" "<<m1.get_Eppps()[i][j]<<" "<<m1.get_Epsps()[i][j]<<" ";
        }
    }
    f<<p[0]<<" "<<p[1];
    f<<"\n";
}
void MC::print_s0() const {
    for(int i=0;i<sample0.size();i++){
        sample0[i].print_sample();
    }
}
// ------------------ end of MC Class ----------------------------------
