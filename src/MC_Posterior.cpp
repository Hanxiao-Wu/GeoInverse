//
// Created by Hanxiao Wu on 10/17/23.
//
#include "MC.cpp"
#include "SAC_IO/sac_stream.cpp"
#include "SAC_IO/sac_io.cpp"
#include "function.cpp"
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "read_write.cpp"
#include <map>
#include <filesystem>
#include <thread>
#include <cstdio>
#include <cstdlib>
using namespace std;
int main(){
    clock_t start,end;
    start = clock();
    string fnm_bin,fnm_control;
    int Nacc;
    cin>>fnm_control>>Nacc;
    int n_group,n_sw, n_hk;
    vector<int> n_disco;
    vector<string> outdir;

    string fgroupnm,finparanm;
    vector<string> fdispnm,frfnm,fhknm0;
    vector<vector<string> > fhknm;
    vector<char> sw_tp,disp_tp;
    vector<double> a, rayp;
    vector<int> idx_monotlyr;
    vector<vector<int>> idisco; // index of layer that monotonically increase
    int n_model,n_search;//#models in one MC search, #MCsearch
    double w_rf,w_sw; //weighting of rf and sw misfit
    vector<double> wPs, wPsPs, wPpPs;//weighting of energy of different phases
    vector<vector<double>> wE;
    double refE=0.;
    int err_flag=0;
    //read in control file
    err_flag = read1(fnm_control,fgroupnm,finparanm,n_group,fdispnm,sw_tp,disp_tp,frfnm,a,rayp,fhknm0,idisco,w_rf,wPs,wPsPs,wPpPs,wE,idx_monotlyr,n_model,n_search,outdir,refE);
    //cout<<"fdispnm:"<<fdispnm.size()<<endl;
    if(err_flag != 1) {
        cout<<"Fail: read in control file!"<<endl;
        exit(0);
    }
    else {
        cout<<"Done: read in control file." << endl;
    }
    w_sw = 1. - w_rf;
    err_flag = read2(fhknm0,fhknm);
    if(err_flag !=1 ){
        cout<<"Fail: read in hk.lst file!"<<endl;
        exit(0);
    }
    else {
        cout<<"Done: read in hk.lst file."<<endl;
    }
    n_hk = fhknm0.size();
    n_sw = fdispnm.size();
    n_disco.clear();
    for(int i=0;i<idisco.size();i++){
        n_disco.push_back(idisco[i].size());
    }
    vector<Group> ref_groups; ref_groups.resize(n_group);
    Model ref_m(n_group);
    ref_m.readin_groups(fgroupnm);
    ref_groups = ref_m.get_groups();
    //cout<<ref_groups[0].get_dep()<<endl;
    ref_m.update_lyrmod();
    ref_m.readin_swdata(fdispnm,disp_tp,sw_tp);
    ref_m.readin_rfdata(frfnm,a,rayp);
    ref_m.readin_hkdata(fhknm,idisco);
    ref_m.set_monot(idx_monotlyr);
    ref_m.set_weight(w_rf,wPs,wPsPs,wPpPs,wE);
    ref_m.set_refE(refE);
    cout<<"set done"<<endl;
    ref_m.update_data(0);
    cout<<"update data"<<endl;
    ref_m.cal_mis();
    cout<<"cal_mis"<<endl;
    ref_m.cal_Energy(1);
    cout<<"cal Energy:"<<ref_m.get_Energy()<<endl;
    fnm_bin = "./"+outdir[0]+"/MC."+outdir[1]+".bin";
    ifstream fbin;
    fbin.open(fnm_bin,ios::in | ios::binary);
    if(!fbin){
        cout<<"error: cannot open "<<fnm_bin<<endl;
        exit(0);
    }
    filesystem::path dir_path("./"+outdir[0]+"/"+outdir[0]);
    try{
        filesystem::create_directory(dir_path);
        cout<<"Successfully created directory:"<<dir_path<<endl;
    }catch (filesystem::filesystem_error &e){
        cerr<<e.what()<<endl;
    }

    Model ave_m(ref_m),min_m(ref_m),minsw_m(ref_m),minrf_m(ref_m),maxE_m(ref_m);
    vector<vector<Group>> Groups;
    Groups.clear();
    vector<Group> tmp_groups;
    tmp_groups.resize(n_group);
    int idx_maxE,idx_minjoint,idx_minrf,idx_mindisp;
    vector<vector<double>> sw_mis_ones;sw_mis_ones.clear();
    sw_mis_ones.resize(n_sw);
    //vector<double> tmp_sw_mis_one; tmp_sw_mis_one.clear();
    vector<double> sw_miss;sw_miss.clear();
    vector<double> rf_miss;rf_miss.clear();
    vector<double> Joint_miss;Joint_miss.clear();
    vector<double> E;
	string dir;
    dir = "./"+outdir[0]+"/"+outdir[0];
	vector<ofstream> facc_group;facc_group.clear();
	ofstream facc_mis;
	string fnm_acc_group;
	for(int i=0;i<n_group;i++){
		fnm_acc_group = dir + "/ACC.group."+to_string(i)+".txt";
		facc_group.emplace_back(std::ofstream{fnm_acc_group});
	}
	facc_mis.open(dir+"/acc.misNenergy.txt");
    double tmp;
	cout<<"begin read in bin:"<<Nacc<<endl;
    for(int k=0;k<Nacc;k++){
        //cout<<k<<endl;
        for(int i=0;i<n_group;i++){
            tmp_groups[i].read_from_bin(fbin);
	    if(i==0) tmp_groups[i].set_dep(ref_groups[i].get_dep());
            //tmp_groups[i].print();
        }

        Groups.push_back(tmp_groups);
        //cout<<"n_sw"<<n_sw<<endl;
        for(int i=0;i<n_sw;i++){
            fbin.read((char *) &tmp,sizeof(double));
            sw_mis_ones[i].push_back(tmp);
            //tmp_sw_mis_one.push_back(tmp);
        }
        //sw_mis_ones.push_back(tmp_sw_mis_one);

        fbin.read((char *) &tmp,sizeof(double));
        //cout<<"sw_mis"<<tmp<<endl;
        sw_miss.push_back(tmp);

        fbin.read((char *) &tmp,sizeof(double));
        rf_miss.push_back(tmp);
        //cout<<"rf_miss:"<<tmp<<endl;
        fbin.read((char *) &tmp, sizeof(double ));
        Joint_miss.push_back(tmp);
/***
        for(int i=0;i<n_hk;i++){
            for(int j=0;j<n_disco[i];j++){
                fbin.read((char *) &tmp,sizeof(double));
                tmp_EPs.push_back(tmp);
                fbin.read((char *) &tmp,sizeof(double));
                tmp_EPpPs.push_back(tmp);
                fbin.read((char *) &tmp, sizeof(double));
                tmp_EPsPs.push_back(tmp);
            }
        }
        EPs.push_back(tmp_EPs);
        EPpPs.push_back(tmp_EPpPs);
        EPsPs.push_back(tmp_EPsPs);
***/
        fbin.read((char *) &tmp,sizeof(double));
        E.push_back(tmp);

    }
    double min_disp=999.,min_rf=999.,min_joint=999.,max_E=-999.;
    vector<double> min_sw_one;
    min_sw_one.resize(n_sw,999.);
    double kai_cri,kai_min=999.;
    int idx_kaimin;
    vector<double> kai_joint;
    //vector<int> a_idx; // store all the index of accepted models
    kai_joint.clear();//a_idx.clear();
    // cal the joint kai & minimum joint kai
    //cout<<"size:"<<sw_miss.size()<<" "<<rf_miss.size()<<" "<<E.size()<<endl;
    min_disp = *min_element(sw_miss.begin(),sw_miss.end());
    min_rf = *min_element(rf_miss.begin(),rf_miss.end());
    max_E = *max_element(E.begin(),E.end());
    for(int i=0;i<n_sw;i++){
        tmp = *min_element(sw_mis_ones[i].begin(),sw_mis_ones[i].end());
        min_sw_one.push_back(tmp);
        cout<<"min "<<i<<"th sw misfit:"<<min_sw_one.back()<<endl;
    }
    cout<<"min mis"<<min_disp<<endl;
    cout<<"min RF "<<min_rf<<endl;
    cout<<"max E "<<max_E<<endl;
    for(int i=0;i<sw_miss.size();i++){
        if(w_sw<0.00001) kai_joint.push_back(rf_miss[i]/min_rf);
        else if(w_rf<0.00001) kai_joint.push_back(sw_miss[i]/min_disp);
        else kai_joint.push_back((w_sw*sw_miss[i]/min_disp+w_rf*rf_miss[i]/min_rf)/(w_sw/min_disp+w_rf/min_rf));
        //if(kai_joint[i]<kai_min) {kai_min = kai_joint[i];}
    }
    kai_min = *min_element(kai_joint.begin(),kai_joint.end());
    cout<<"done:calculate kai_joint"<<kai_joint.size()<<endl;
    // find the critical joint kai
    if(w_sw<0.00001) kai_cri = kai_min + threshold_mis/min_rf;
    else if(w_rf<0.00001) kai_cri = kai_min + threshold_mis/min_disp;
    else {
		if(kai_min<threshold_mis) kai_cri = 2*kai_min;
		else kai_cri = kai_min + threshold_mis;
	}
	cout<<"minDisp:"<<min_disp<<endl;
	cout<<"minRF:"<<min_rf<<endl;
    cout<<"threshold: kai_cri "<<kai_min<<" "<<kai_cri<<endl;
    cout<<"threshold: Energy "<<threshold_E*max_E<<endl;
    // posterior model space
	double E_cri = threshold_E*max_E;
    min_disp = 999.;min_rf = 999.;min_joint = 999.;max_E = -999.;
    vector<vector<Group> > acc_groups;acc_groups.clear();
    vector<int> flag;
    flag.resize(n_sw,0.0);
    int Flag = 1.;
    for(int i=0;i<sw_miss.size();i++){
        Flag = 1.;flag.clear();flag.resize(n_sw,0.0);
        if(kai_joint[i]<kai_cri && E[i]>E_cri ){
            for(int k=0;k<n_sw;k++){
                if(sw_mis_ones[k][i]<min_sw_one[k]+threshold_sw_one) flag[k]=1;
                //Flag = Flag * flag[k];
				Flag = 1; // to be fixed
            }
            if(Flag==1) {
				facc_mis<<kai_cri<<" "<<kai_joint[i]<<" "<<E_cri<<" "<<E[i]<<" ";
				facc_mis<<rf_miss[i]<<" ";
				facc_mis<<sw_miss[i]<<" ";
				for(int k=0;k<n_sw;k++){ facc_mis<<sw_mis_ones[k][i]<<" ";}
				facc_mis<<endl;
                acc_groups.push_back(Groups[i]);
				for(int j=0;j<Groups[i].size();j++) {
					Groups[i][j].write_v_to_file(facc_group[j]);
				}
                if (sw_miss[i] < min_disp) {
                    min_disp = sw_miss[i];
                    idx_mindisp = i;
                }
                if (rf_miss[i] < min_rf) {
                    min_rf = rf_miss[i];
                    idx_minrf = i;
                }
                if (Joint_miss[i] < min_joint) {
                    min_joint = Joint_miss[i];
                    idx_minjoint = i;
                }
                if (E[i] > max_E) {
                    max_E = E[i];
                    idx_maxE = i;
                }
            }

        }
    }
	facc_mis.close();
	for(int i=0;i<n_group;i++){facc_group[i].close();}
    cout<<"final accepted:"<<acc_groups.size()<<endl;
    cout<<"maxE: "<<max_E<<endl;
    maxE_m.set_group(Groups[idx_maxE]);
    maxE_m.update_lyrmod();
    maxE_m.update_data(0);
    maxE_m.cal_Energy(1);
    maxE_m.cal_mis();

    cout<<"min joint misfit: "<<min_joint<<endl;
    min_m.set_group(Groups[idx_minjoint]);
    min_m.update_lyrmod();
    min_m.update_data(0);
    min_m.cal_Energy(1);
    min_m.cal_mis();


    cout<<"min RF misfit: "<<min_rf<<endl;
    minrf_m.set_group(Groups[idx_minrf]);
    minrf_m.update_lyrmod();
    minrf_m.update_data(0);
    minrf_m.cal_Energy(1);
    minrf_m.cal_mis();

    cout<<"min disp misfit: "<<min_disp<<endl;
    minsw_m.set_group(Groups[idx_mindisp]);
    minsw_m.update_lyrmod();
    minsw_m.update_data(0);
    minsw_m.cal_Energy(1);
    minsw_m.cal_mis();

    //get average model
    if(acc_groups.size()>1) {
        ave_m = average_m(ref_m, acc_groups);
		//std_m = std_m(ref_m,acc_groups);
    }
    else if (acc_groups.size()==1){
        ave_m.set_group(acc_groups[0]);
		//std_m.set_group(acc_groups[0]);
    }
    else {
        cout<<"no accepted models left"<<endl;
        return 0;
    }
    cout<<"done: cal ave_m"<<endl;
    ave_m.update_lyrmod();
    ave_m.update_data(0);
    ave_m.cal_Energy(1);
    ave_m.cal_mis();
    cout<<"ave_m"<<ave_m.get_Energy()<<endl;
	//std_m.update_lyrmod();
	//std_m.update_data(0);
	//std_m.cal_Energy(1);
	//std_m.cal_mis();
    //write model to file
    //string dir;
    //dir = "./"+outdir[0]+"/"+outdir[0];
    maxE_m.write_groups_to_file(dir,outdir[1]+".maxE");
    maxE_m.write_lyr_to_file(dir,outdir[1]+".maxE");
    maxE_m.write_data_to_file(dir+"/MC."+outdir[1]+".maxE",0);
    maxE_m.write_info(dir,outdir[1]+".maxE");

    minrf_m.write_groups_to_file(dir,outdir[1]+".minrf");
    minrf_m.write_lyr_to_file(dir,outdir[1]+".minrf");
    minrf_m.write_data_to_file(dir+"/MC."+outdir[1]+".minrf",0);
    minrf_m.write_info(dir,outdir[1]+".minrf");

    minsw_m.write_groups_to_file(dir,outdir[1]+".minsw");
    minsw_m.write_lyr_to_file(dir,outdir[1]+".minsw");
    minsw_m.write_data_to_file(dir+"/MC."+outdir[1]+".minsw",0);
    minsw_m.write_info(dir,outdir[1]+".minsw");

    min_m.write_groups_to_file(dir,outdir[1]+".min");
    min_m.write_lyr_to_file(dir,outdir[1]+".min");
    min_m.write_data_to_file(dir+"/MC."+outdir[1]+".min",0);
    min_m.write_info(dir,outdir[1]+".min");

    ave_m.write_groups_to_file(dir,outdir[1] + ".ave");
    ave_m.write_lyr_to_file(dir,outdir[1]+".ave");
    ave_m.write_data_to_file(dir+"/MC."+outdir[1]+".ave",0);
    ave_m.write_info(dir,outdir[1]+".ave");

    //std_m.write_groups_to_file(dir,outdir[1] + "std");
	//std_m.write_lyr_to_file(dir,outdir[1]+".std");
	//std_m.write_data_to_file(dir+"/MC."+outdir[1]+".std",0);
	//std_m.write_info(dir,outdir[1]+".std");
	end = clock();
    cout<<"posterior process time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<endl;
    return 0;
}
