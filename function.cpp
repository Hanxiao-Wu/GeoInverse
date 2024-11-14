//
// Created by Hanxiao Wu on 7/15/23.
//
#include "function.h"

using namespace std;

double get_vp_brocher(const double &vs){
    double vp;
    vp = 0.9409 + 2.0947*vs-0.8206*pow(vs,2)+0.2683*pow(vs,3)-0.0251*pow(vs,4);
    return vp;
}
double get_rho_brocher(const double &v,const char &flag){
    double rho;
	switch(flag){
		case 's':
			rho = 1.22679+1.53201*v-0.83668*pow(v,2)+0.20673*pow(v,3)-0.01656*pow(v,4);
			break;
		case 'p':
			rho = 1.6612*v-0.4721*v*v+0.0671*v*v*v-0.0043*v*v*v*v+0.000106*v*v*v*v*v;
			break;
	}
    return rho;
}
double get_rho_hacker(const double &vs){
    double rho;
    rho = 3.42+0.01*100*(vs-4.5)/4.5;
    return rho;
}

double rand_uniform(const double &left, const double &right,const double &value){
    //default_random_engine generator;
    random_device generator;
    uniform_real_distribution<double> uni_rand01(0.0,1.0);
    double range = right - left;
    double tmp = value;
    if (range>0) tmp = uni_rand01(generator)*range + left;
    return tmp;
}
AdditionalPert rand_uniform_add(const AdditionalPert &left, const AdditionalPert &right, const AdditionalPert &v){
    double tmp_top, tmp_bot, tmp_v;
    tmp_top = rand_uniform(left.get_top(),right.get_top(),v.get_top());
    tmp_bot = rand_uniform(left.get_bot(),right.get_bot(),v.get_bot());
    tmp_v = rand_uniform(left.get_value(),right.get_value(),v.get_value());
    AdditionalPert tmp(tmp_top,tmp_bot,tmp_v);
    return tmp;
}

double rand_gaussian(const double &mean, const double &sigma, const double &left, const double &right){
    //default_random_engine generator;
    random_device generator;
    normal_distribution<double> rand_norm(0.0,sigma);
    double tmp = mean;
    int nbad=0;
    if(right-left>0) {
        while (true) {
            tmp = mean + rand_norm(generator);
            rand_norm.reset();
            if(tmp > right) tmp = right - (tmp-right);
            else if (tmp < left)  tmp = left + (left-tmp);
            if (tmp <= right && tmp >= left) break;
            nbad++; //
            if (nbad > 1000) {cout << "couldn't generate a model within model space" << endl;return 1;}
        }
        if(nbad>0) cout<<"#bad samples:"<<nbad<<endl;
    }
    return tmp;
}
AdditionalPert rand_gaussian_add(const AdditionalPert &mean, const AdditionalPert &sigma, const AdditionalPert &left, const AdditionalPert &right) {
    double tmp_top, tmp_bot, tmp_value;
    tmp_top = rand_gaussian(mean.get_top(),sigma.get_top(),left.get_top(),right.get_top());
    tmp_bot = rand_gaussian(mean.get_bot(),sigma.get_bot(),left.get_bot(),right.get_bot());
    tmp_value = rand_gaussian(mean.get_value(),sigma.get_value(),left.get_value(),right.get_value());
    return {tmp_top, tmp_bot, tmp_value};
}

double add_pertb(const vector<AdditionalPert> &add_pert, const double &Thk, const double &dep, double value){
    double tmp_top, tmp_bot, tmp_dv;
    for (const auto & j : add_pert){
        tmp_top = j.get_top();
        tmp_bot = j.get_bot();
        tmp_dv = j.get_value();
        if(dep>=Thk*tmp_top && dep<=Thk*tmp_bot) value = value +tmp_dv;

    }
    return value;
}

int find_closest_value(const vector<double> &v, const double &x) {
    auto it = lower_bound(v.begin(),v.end(),x);
    int idx = it - v.begin();
    if(v.empty()){cout<<"find_closest_value: size of vector is wrong"<<endl;exit(0);}
    if(idx == 0) return idx;
    else if(idx == v.size()) return idx-1;
    else if(abs(x-v[idx-1]) < abs(x-v[idx])) return idx-1;
    else if(abs(x-v[idx-1]) >= abs(x-v[idx])) return idx;
    else {cout<<"couldn't find a closest value"<<endl;return -1;}
}

vector<double> joint_lyr(vector<double> v1, vector<double> v2){
    v1.insert(v1.end(),v2.begin(),v2.end());
    return v1;
}

double Gm_d(const vector<double> &pre, const vector<double> &obs, const vector<double> &err){
    double tmp=0.;
    for(int i=0;i<pre.size();i++){
        tmp = tmp + pow((pre[i]-obs[i])/err[i],2);
    }
    //if(tmp>0.){tmp = sqrt(tmp/pre.size());}
    return tmp;
}
double reduceLbyhand(double mis){
    //cout<<"reduceL:mis"<<mis;
    if(mis>30) mis = sqrt(mis*30);
    //if(mis>30) mis = sqrt(mis*30);
    //if(mis>30) mis = sqrt(mis*30);
    //if(mis>30) mis = sqrt(mis*30);
    //cout<<" "<<mis<<endl;
    return mis;
}

void get_bound(const vector<double> &v, const vector<double> &radius, vector<double> &left, vector<double> &right){
// radius must be positive
    if(radius.empty()){left = v;right=v;}
    else{
        left = v - radius;
        right = v + radius;
    }
}
void get_bound(const vector<AdditionalPert> &v, const vector<AdditionalPert> &radius, vector<AdditionalPert> &left, vector<AdditionalPert> &right){
    if(radius.empty()){left = v;right=v;}
    else{
        left = v - radius;
        right = v + radius;
    }
}

bool check_monot(const vector<double> &v){
    bool inc = true;
    for(int i=1;i<v.size();i++){
        //cout<<"check_monot: "<<v.size()<<endl;
        // To check if array is not increasing
        if (v[i] < v[i - 1]) {
            inc = false;
        }
    }
    return inc;
}

void Split(const string &s, vector<string> &v, const string& del = " ")
{
    int start, end = -1*del.size();
    v.clear();
    do {
        start = end + del.size();
        end = s.find(del, start);
        v.push_back(s.substr(start,end-start));
    } while (end != -1);
}

Model average_m(const Model &m, const vector<vector<Group> > &groups){
    vector<Group> tmp_group = groups[0];
    Model tmp(m);
    for(int i=1;i<groups.size();i++){
        for(int j=0;j<tmp_group.size();j++){
            tmp_group[j] = tmp_group[j] + groups[i][j];
            //cout<<i<<"th groups is added"<<endl;
            //cout<<i<<"th group's thk"<<tmp_group[j].get_thk()<<endl;
        }
    }
    cout<<"done: sum up"<<endl;
    for(auto & j : tmp_group){
        j = j/groups.size();
        //cout<<"divided group thk: "<<j.get_thk()<<endl;
        j.update_lyrmodel();
    }
    cout<<"done: divided"<<endl;
    tmp.set_group(tmp_group);
    return tmp;
}
/***
template<typename T>
T variance(const std::vector<T> &vec){
	const size_T sz = vec.size();
	if (sz<=1){
		return 0.;
	}

	const T mean = std::accumulate(vec.begin(),vec.end(),0.0)/sz;

	auto variance_func = [&mean,&sz] (T accumulator, const T& val){
		return accumulator + ((val-mean)*(val-mean)/(sz-1));
	};
	return std::accumulate(vec.begin(),vec.end(),0.0,variance_func);
}

Model std_m(const Model &m, const vector<vector<Group> > &groups){
	vector<Group> tmp_group = groups[0];
	Model tmp(m);
	for(int i=0;i<groups.size();i++){
		tmp_group[i] = variance(groups[i]);
	tmp.set_group(tmp_group);
	return tmp;
}
***/	
