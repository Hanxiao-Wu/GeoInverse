//
// Created by Hanxiao Wu on 7/14/23.
//
#include "MC.h"
#include <iostream>
#include <fstream>
//#include <cmath>
#include "GlobalParameter.h"
#include "function.h"
#include <random>
#include <vector>
#include "SAC_IO/sac_stream.hpp"
using namespace std;
int main() {
    double dep;
    int n_group = 3;
    vector<Group> groups;
    vector<Sample> s;
    for(int i=0;i<n_group;i++) {
        Group tmp(n_group);
        Sample tmps;
        groups.push_back(tmp);
        s.push_back(tmps);
    }
    ifstream f;
    f.open("test.mod");
    if(f.is_open()){
        string line; vector<string> v;
        while(getline(f,line)){
            Split(line, v, " ");
            int idx = atoi(v[0].c_str()); // index of group
            groups[idx].readin(v);
            s[idx].Initial(v);
        }
    }
    f.close();

    for(int i=0;i<n_group;i++){
        if(i==0) {groups[i].set_dep();}
        else {
            dep = groups[i-1].get_dep() + groups[i-1].get_thk();
            groups[i].set_dep(dep);
        }
    }
//update lyr model
    for(int i=0;i<n_group;i++){
        groups[i].update_lyrmodel();
        groups[i].write_lyr();
    }

    f.open("in.para");
    if(f.is_open()){
        string line; vector<string> v;
        while(getline(f,line)){
            Split(line, v, " ");
            int idx = atoi(v[0].c_str()); // index of group
            s[idx].readin(v);
        }
    }
    for(int i=0;i<n_group;i++){
        cout<<i<<"th sample:"<<endl;
        s[i].print();
    }
    // generate starting point
    int nbad=-1;
    double Thk=0.;
    nbad++;
    //if(nbad>Nbad_threshold){cout<<"couldn't find a model satisfying"<<endl;}
    for(int i=0;i<s.size();i++) {
        s[i].rand_start_sample();
        if(i==n_group-1) s[i].set_thk(176.-Thk);
        Thk = Thk + s[i].get_thk();
        groups[i].update_from_sample(s[i]);
        if(i>0){
            dep = groups[i-1].get_dep() + groups[i-1].get_thk();
            groups[i].set_dep(dep);
        }
        groups[i].update_lyrmodel();
        groups[i].write_lyr();
    }




}
