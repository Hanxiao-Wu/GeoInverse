//
// Created by Hanxiao Wu on 7/21/23.
//
#include <iostream>
#include "MC.h"
#include <chrono>

using namespace chrono;
using namespace std;
int main(){/*
    auto start = chrono::system_clock::now();
    vector<double> para;
    vector<double> radius,step;
    double thk,vpvs,Rho;
    double thk_radius,vpvs_radius,Rho_radius;
    double thk_step, vpvs_step, Rho_step;
    vector<AdditionalPert> para_add_vs, para_add_vp, para_add_rho;
    vector<AdditionalPert> radius_add_vs, radius_add_vp, radius_add_rho;
    vector<AdditionalPert> step_add_vs, step_add_vp, step_add_rho;

    para.clear();para.push_back(1.5);para.push_back(2.5);
    radius.clear();radius.push_back(0.5);radius.push_back(0.5);
    step.clear();step.push_back(0.05);step.push_back(0.05);

    thk = 30.; vpvs = 1.75; Rho = 3.9;
    thk_radius = 5.; vpvs_radius = 0.15; Rho_radius = 0.03;
    thk_step = 0.5; vpvs_step = 0.015; Rho_step = 0.003;

    AdditionalPert tmp(0.5,1.0,0.02);
    para_add_vs.clear(); para_add_vs.push_back(tmp);
    tmp.set(0.,0.,0.01);
    radius_add_vs.clear(); radius_add_vs.push_back(tmp);
    tmp.set(0.,0.,0.001);
    step_add_vs.clear();step_add_vs.push_back(tmp);

    tmp.set(0.,0.5,1.71);
    para_add_vp.clear(); para_add_vp.push_back(tmp);
    tmp.set(0.,0.,0.15);
    radius_add_vp.clear();radius_add_vp.push_back(tmp);
    tmp.set(0.,0.,0.015);
    step_add_vp.push_back(tmp);

    para_add_rho.clear();radius_add_rho.clear();step_add_rho.clear();

    vector<Sample> ss;ss.clear();
    Sample s(para,radius,step,thk,thk_radius,thk_step,vpvs,vpvs_radius,vpvs_step,Rho,Rho_radius,Rho_step,para_add_vs,radius_add_vs,step_add_vs,para_add_vp,radius_add_vp,step_add_vp,para_add_rho,radius_add_rho,step_add_rho);
    //s.rand_start_sample();
    int N = 50*5000; // the radius should be within 10 times of the step


    for(int i=0;i<N;i++){
        s.gen_new_sample();
        //s.rand_start_sample();
        ss.push_back(s);
        s.print();
    }
    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end-start);
    cout<<double(duration.count())*microseconds::period::num / microseconds::period::den<<"s"<<endl;
*/
    cout<<"hello world"<<endl;
    return 0;
}