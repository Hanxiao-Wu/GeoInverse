//
// Created by Hanxiao Wu on 8/4/23.
//
#include "MC.h"
//#include "MC.cpp"
#include "GlobalParameter.h"
//#include "SAC_IO/sac_stream.hpp"
//#include "SAC_IO/sac_io.hpp"
#include "SAC_IO/sac_stream.cpp"
#include "SAC_IO/sac_io.cpp"
#include "function.h"
//#include "function.cpp"
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
using namespace std;
int main(){
    int n_group = 3;
    //pi=atan(1)*4;
    cout<<pi<<endl;
    Model m(n_group);
    //set up models
    m.readin_groups("test.mod");
    m.update_lyrmod();
    Lyr_model tmp_lyrm = m.get_lyrm();
    //cout<<"#lyr:"<<tmp_lyrm.get_nlyr()<<endl;
    //tmp_lyrm.print_lyr();
    // read in swdata
    vector<string> fname;
    vector<char> tp, sw_tp;
    fname.clear();tp.clear();sw_tp.clear();
    fname.push_back("tar_g.dat");tp.push_back('g');sw_tp.push_back('R');
    fname.push_back("tar_p.dat");tp.push_back('p');sw_tp.push_back('R');
    fname.push_back("tar_e.dat");tp.push_back('e');sw_tp.push_back('R');
    fname.push_back("tar_a.dat");tp.push_back('a');sw_tp.push_back('R');
    m.readin_swdata(fname,tp,sw_tp);
    vector<string> rfname;
    vector<double> a,rayp;
    rfname.clear();rfname.push_back("tar_rf.dat");
    a.clear();a.push_back(2.5);
    rayp.clear();rayp.push_back(0.06);
    m.readin_rfdata(rfname,a,rayp);
    cout<<"read in data:done"<<endl;
    //m.print_swdata();
    //forward
    m.update_SWdata();
    cout<<"update sw data done"<<endl;
    m.print_swdata();

    m.update_RFdata();
    cout<<"update rf data done"<<endl;
    m.print_rfdata();
    /*
    //sample
    vector<Sample> s0;s0.clear();
    for (int i=0;i<n_group;i++){
        Sample tmp;
        s0.push_back(tmp);
    }

    ifstream f;
    f.open("test.mod");
    if(f.is_open()){
        string line; vector<string> v;
        while(getline(f,line)){
            Split(line, v, " ");
            int idx = atoi(v[0].c_str()); // index of group
            s0[idx].Initial(v);
        }
    }
    f.close();
    cout<<"Initial sample"<<endl;
    f.open("in.para");
    if(f.is_open()){
        string line; vector<string> v;
        while(getline(f,line)){
            Split(line, v, " ");
            int idx = atoi(v[0].c_str()); // index of group
            s0[idx].readin(v);
        }
    }
    cout<<"read in in.para"<<endl;



    int N = 5000;
    MC mc(s0,m,N);
    cout<<"construct mc"<<endl;

    mc.rand_starting_point();
    cout<<"get random starting point"<<endl;
    Model cm = mc.get_currentmodel();

    Lyr_model lym = cm.get_lyrm();
    lym.print_lyr();
    //mc.gen_sample1();
*/
    cout<<"hello world"<<endl;
}