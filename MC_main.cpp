//
// Created by Hanxiao Wu on 8/14/23.
//
#include <time.h>
#include "MC.cpp"
//#include "MC.cpp"
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
using namespace std;
int main(){
    clock_t start,end1,end2;
    start = clock();
    string f_control;
    int n_thread=1; // default
    cin>>f_control>>n_thread; // input name and n_thread

    vector<string> outdir;
    string fgroupnm,finparanm;
    vector<string> fdispnm,frfnm,fhknm0;
    vector<vector<string> > fhknm;
    vector<char> sw_tp,disp_tp;
    vector<double> a, rayp;
    vector<int> idx_monotlyr;
    vector<vector<int>> idisco; // index of layer that monotonically increase
    int n_group; //#groups
    int n_model,n_search;//#models in one MC search, #MCsearch
    double w_rf,w_sw; //weighting of rf and sw misfit
    vector<double> wPs, wPsPs, wPpPs;//weighting of energy of different phases
    vector<vector<double>> wE;
    double refE=0.;
    int err_flag=0;
	int flag_forward=0;
    //read in control file
    err_flag = read1(f_control,fgroupnm,finparanm,n_group,fdispnm,sw_tp,disp_tp,frfnm,a,rayp,fhknm0,idisco,w_rf,wPs,wPsPs,wPpPs,wE,idx_monotlyr,n_model,n_search,outdir,refE);
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
    //create dir
    filesystem::path dir_path("./"+outdir[0]);
    try{
        filesystem::create_directory(dir_path);
        cout<<"Successfully created directory:"<<dir_path<<endl;
    } catch (filesystem::filesystem_error &e){
        cerr<<e.what()<<endl;
    }
    //
    Model ref_m(n_group);
    vector<Model> models; //store all the models
    models.clear();
    //int k=0,idx_mc;
    int idx_mindisp=0, idx_minrf=0, idx_minjoint=0, idx_maxE=0;
    //reference model
    ref_m.readin_groups(fgroupnm);
    //ref_m.print_groups();
    ref_m.update_lyrmod();
    //read in data
    cout<<"sizeOfModel"<<sizeof(Model)<<endl;
    ref_m.readin_swdata(fdispnm,disp_tp,sw_tp);
    ref_m.readin_rfdata(frfnm,a,rayp);
    ref_m.readin_hkdata(fhknm,idisco);
    cout<<"done: read in data"<<endl;
    ref_m.set_weight(w_rf,wPs,wPsPs,wPpPs,wE);
    //cout<<"wE:"<<wE.size()<<" "<<wE[0].size()<<endl;
    //cout<<"refE:"<<refE<<endl;
    ref_m.set_refE(refE);
    //cout<<"refE:"<<ref_m.get_refE()<<endl;
    ref_m.set_monot(idx_monotlyr);
    ref_m.write_groups_to_file(outdir[0],outdir[1]);
    ref_m.write_lyr_to_file(outdir[0],outdir[1]);
    cout<<"done: write lyr to file"<<endl;
    //models.push_back(ref_m);
    //if(n_model<0) {
    cout<<"begin forward calculation:"<<endl;
    cout<<"update data"<<endl;
    ref_m.update_data(0);
    cout<<"update Hkdata"<<endl;
    if(n_model<0) {flag_forward=1;}
    ref_m.cal_Energy(flag_forward);
	cout<<"cal_energy"<<endl;
    ref_m.cal_mis();
	ref_m.write_info(outdir[0],outdir[1]);
    cout<<"cal_mis"<<endl;
    ref_m.cal_RFmisfit();
    ref_m.cal_Jmisfit();
    ref_m.print_swdata();
    ref_m.print_rfdata();
    ref_m.update_data(flag_forward);
	ref_m.write_data_to_file(outdir[0]+"/forward",flag_forward); //todo: need to be fixed
    cout<<"done: forward calculation."<<endl;
    cout<<"Energy:"<<ref_m.get_Energy()<<endl;
    cout<<"Misfit:"<<ref_m.get_Jointmis()<<endl;
    if(n_model<0) return 1;
    

    models.clear();
    //models.reserve(n_search*(n_model+1));
    //start inversion
    ofstream fout_s;
    ofstream fout_bin;
    fout_bin.open("./"+outdir[0]+"/MC."+outdir[1]+".bin",ios::out|ios::binary);
    fout_s.open("./"+outdir[0]+"/MC."+outdir[1]+".para",ios::out);
    cout<<"done: open test/MC.*.para file"<<endl;
    int idx_thread;
    if(n_search<0){
		for (int i = 0; i < -n_search; i++) {
	         cout<<"begin prior sampling"<<endl;
		     MC mc(n_model);
		     mc.set_weight(w_rf);
		     mc.Initial_m0(n_group, fgroupnm, fdispnm, disp_tp, sw_tp, frfnm, a, rayp, fhknm        , idisco,w_rf, wPs, wPsPs, wPpPs, wE, refE,idx_monotlyr);
			mc.Initial_s0(fgroupnm, n_group);
			if(mc.read_s0(finparanm)==0){exit(0);}
			mc.rand_starting_point();
			for (int j = 1; j < n_model; j++) {
			    mc.gen_good_sample();
                mc.cal_p();
				mc.update_m1();
                mc.accept_model(fout_bin,models);
                idx_thread = DEFAULTi;
                mc.write_samples_to_file(idx_thread,0,j,fout_s);
			}
		}   
         fout_s.close();
         fout_bin.close();
         return 1;
     }   

    //begin inversion
    {
        for (int i = 0; i < n_search; i++) {
            //construct a MC machine
            MC mc(n_model);
            cout<<i<<"th mc is constructed"<<endl;
            //cout<<"#model is "<<n_model<<" in "<<i<<"th search"<<endl;
            mc.set_weight(w_rf);
            //cout<<"done:set w_rf"<<endl;
            //initial model(read in reference model)
            mc.Initial_m0(n_group, fgroupnm, fdispnm, disp_tp, sw_tp, frfnm, a, rayp, fhknm, idisco,w_rf, wPs, wPsPs, wPpPs, wE, refE,
                          idx_monotlyr);
            //cout<<"refE:"<<refE<<endl;
            //cout<<i<<"th search: initial m0 done"<<endl;
            //initial sample
            mc.Initial_s0(fgroupnm, n_group);
            //cout<<i<<"th search: initial s0 done:"<<endl;
            if(mc.read_s0(finparanm)==0){exit(0);}
            //mc.print_s0();
            cout<<"done:read in para"<<endl;
            //get a random starting point
            if(mc.rand_starting_point()==0){i--;continue;}
            cout<<"done: get a starting point."<<endl;
            //mc.write_samples_to_file(idx_thread,i,0,fout_s);
            //begin search
            for (int j = 1; j < n_model; j++) {
                //cout<<"begin gen_good_sample:"<<j<<endl;
                if(mc.gen_good_sample()==0){
                    i--;
                    continue;
                } //set m1
				mc.update_m1();
                //write sample to file
                // calculate chance of being accepted
                //cout<<"cal_p"<<j<<endl;
                mc.cal_p(); //p=L1/L0
                //decide if accept model
                //cout<<"done:cal_p"<<j<<endl;
                mc.accept_model(fout_bin,models);
                //Model tmp_m =mc.get_models().back();
                //models.push_back(tmp_m);
                idx_thread = DEFAULTi; //todo
                //cout<<"write sample to file"<<j<<endl;
                mc.write_samples_to_file(idx_thread,i,j,fout_s);
                cout<<i<<"th search,model:"<<j<<endl;
                //cout<<"done:write sample to file"<<j<<endl;
            }
            //store all the models & info
            //cout<<mc.get_models().size()<<endl;
            //vector<Model> tmp_m = mc.get_models();
            //models.insert(models.end(),tmp_m.begin()+1,tmp_m.end());
             cout<<i<<"th search is done, models inserted"<<endl;
            //a_idx.insert(a_idx.end(),mc.get_adix().begin(),mc.get_adix().end());
        }
    }
    fout_s.close();
    fout_bin.close();
    cout<<"close fouts"<<endl;
    end1 = clock();

/***

    //posterior process
    Model minrf_m(ref_m),minsw_m(ref_m),min_m(ref_m),maxE_m(ref_m), ave_m(ref_m);
    double min_disp=999.,min_rf=999.,min_joint=999.,max_E=-999.;
    double kai_cri,kai_min=999.;
    vector<double> sw_mis,rf_mis,joint_mis,E;
    vector<vector<Group>> Groups;Groups.clear();
    vector<double> kai_joint;
    sw_mis.clear();rf_mis.clear();joint_mis.clear();E.clear();
    // cal the joint kai & minimum joint kai
    cout<<"begin posterior process:"<<models.size()<<endl;
    if(wE[0][0]<0){
        for(int i=0;i<models.size();i++){
            sw_mis.push_back(models[i].get_Dispmis());
            rf_mis.push_back(models[i].get_RFmis());
            joint_mis.push_back(models[i].get_Jointmis());
            E.push_back(999.);
            Groups.push_back(models[i].get_groups());
        }
    }
    for(int i=0;i<models.size();i++) {
        if(models[i].get_flag()==1 && models[i].get_flagE()==1) {
            sw_mis.push_back(models[i].get_Dispmis());
            rf_mis.push_back(models[i].get_RFmis());
            joint_mis.push_back(models[i].get_Jointmis());
            E.push_back(models[i].get_Energy());
            Groups.push_back(models[i].get_groups());
            //cout<<sw_mis.back()<<" "<<rf_mis.back()<<" "<<E.back()<<endl;
        }
        //cout<<"mis:"<<models[i].get_Dispmis()<<" "<<models[i].get_RFmis()<<endl;
    }
    min_disp = *min_element(sw_mis.begin(),sw_mis.end());
    min_rf = *min_element(rf_mis.begin(),rf_mis.end());
    max_E = *max_element(E.begin(),E.end());
    cout<<"min mis"<<min_disp<<endl;
    cout<<"min RF "<<min_rf<<endl;
    cout<<"max E "<<max_E<<endl;
    kai_joint.clear();
    for(int i=0;i<sw_mis.size();i++){
        if(w_sw<0.00001) kai_joint.push_back(rf_mis[i]/min_rf);
        else if(w_rf<0.00001) kai_joint.push_back(sw_mis[i]/min_disp);
        else kai_joint.push_back((w_sw*sw_mis[i]/min_disp+w_rf*rf_mis[i]/min_rf)/(w_sw/min_disp+w_rf/min_rf));
        //if(kai_joint[i]<kai_min) {kai_min = kai_joint[i];}
        //cout<<"kai_joint:"<<kai_joint[i]<<endl;
    }
    kai_min = *min_element(kai_joint.begin(),kai_joint.end());
    cout<<"done:calculate kai_joint"<<kai_joint.size()<<endl;
    // find the critical joint kai
    if(w_sw<0.00001) kai_cri = kai_min + threshold_mis/min_rf;
    else if(w_rf<0.00001) kai_cri = kai_min + threshold_mis/min_disp;
    else kai_cri = kai_min + threshold_mis;
    cout<<"done: kai_cri "<<kai_min<<" "<<kai_cri<<endl;
    // posterior model space
    min_disp = 999.;min_rf = 999.;min_joint = 999.;//max_E = -999.;
    vector<vector<Group>> acc_groups;acc_groups.clear();
    for(int i=0;i<sw_mis.size();i++){
        //cout<<kai_joint[i]<<" "<<kai_cri<<endl;
        //cout<<models[idx].get_Energy()<<" "<<0.9*max_E<<endl;
        if(kai_joint[i]<kai_cri && E[i]>threshold_E*max_E){
            //for(int j=0;j<models[i].get_Eps().size();j++){
              //  for(int k=0;k<models[i].get_Eps()[j].size();k++){
                //    if(models[i].get_Eps()[j][k]>0 && models[i].get_Eppps()[j][k]>0 && models[i].get_Epsps()[j][k]<0){
                        acc_groups.push_back(Groups[i]);
                        if(sw_mis[i]<min_disp){
                            min_disp = sw_mis[i];
                            idx_mindisp = i;
                        }
                        if(rf_mis[i]<min_rf){
                            min_rf = rf_mis[i];
                            idx_minrf = i;
                        }
                        if(joint_mis[i]<min_joint){
                            min_joint = joint_mis[i];
                            idx_minjoint = i;
                        }
                        if(E[i]>max_E){
                            max_E = E[i];
                            idx_maxE = i;
                        }
                    //}
            //    }
          //  }
        }
    }
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
    minsw_m.update_data(1);
    minsw_m.cal_Energy(1);
    minsw_m.cal_mis();

    //get average model
    if(acc_groups.size()>1) {
        ave_m = average_m(ref_m, acc_groups);
    }
    else if (acc_groups.size()==1){
        ave_m.set_group(acc_groups[0]);
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

    //write model to file
    maxE_m.write_groups_to_file(outdir[0],outdir[1]+".maxE");
    maxE_m.write_lyr_to_file(outdir[0],outdir[1]+".maxE");
    maxE_m.write_data_to_file(outdir[0]+"/MC."+outdir[1]+".maxE",0);
    maxE_m.write_info(outdir[0],outdir[1]+".maxE");

    minrf_m.write_groups_to_file(outdir[0],outdir[1]+".minrf");
    minrf_m.write_lyr_to_file(outdir[0],outdir[1]+".minrf");
    minrf_m.write_data_to_file(outdir[0]+"/MC."+outdir[1]+".minrf",0);
    minrf_m.write_info(outdir[0],outdir[1]+".minrf");

    minsw_m.write_groups_to_file(outdir[0],outdir[1]+".minsw");
    minsw_m.write_lyr_to_file(outdir[0],outdir[1]+".minsw");
    minsw_m.write_data_to_file(outdir[0]+"/MC."+outdir[1]+".minsw",0);
    minsw_m.write_info(outdir[0],outdir[1]+".minsw");

    min_m.write_groups_to_file(outdir[0],outdir[1]+".min");
    min_m.write_lyr_to_file(outdir[0],outdir[1]+".min");
    min_m.write_data_to_file(outdir[0]+"/MC."+outdir[1]+".min",0);
    min_m.write_info(outdir[0],outdir[1]+".min");

    ave_m.write_groups_to_file(outdir[0],outdir[1] + ".ave");
    ave_m.write_lyr_to_file(outdir[0],outdir[1]+".ave");
    ave_m.write_data_to_file(outdir[0]+"/MC."+outdir[1]+".ave",0);
    ave_m.write_info(outdir[0],outdir[1]+".ave");
***/
    //end2 = clock();
    cout<<"inversion time = "<<double(end1-start)/CLOCKS_PER_SEC<<"s"<<endl;
    //cout<<"posterior process time = "<<double(end2-end1)/CLOCKS_PER_SEC<<"s"<<endl;
    //cout<<"total time = "<<double(end2-start)/CLOCKS_PER_SEC<<"s"<<endl;
    return 0;
}
