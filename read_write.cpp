//
// Created by Hanxiao Wu on 8/14/23.
//

#include "read_write.h"
#include "function.h"
int read1(const string& f_control, string& fgroupnm, string &finparanm, int &n_group, vector<string> &fdispnm, vector<char> &sw_tp, vector<char> &disp_tp, vector<string> &frfnm,vector<double> & a, vector<double> &rayp, vector<string> &fhknm,vector<vector<int>> & idisco, double &w_rf, vector<double> &wPs,vector<double> &wPsPs,vector<double> &wPpPs, vector<vector<double>> &wE,vector<int> &idx_monotlyr, int &n_model, int &n_search, vector<string> &outdir, double &refE){
    ifstream f;
    f.open(f_control);
    if (not f.is_open()) {
        cout << "#####control file does not exist:" << f_control << endl;
        exit(0);
    }
// initialize numbers for switch
    map<string, int> numbers;

    numbers["model"] = 1;
    numbers["para"] = 2;
    numbers["disp"] = 3;
    numbers["rf"] = 4;
    numbers["rfweight"] = 5;
    numbers["hk"] = 6;
    numbers["hkweight"] = 7;
    numbers["monol"] = 9;
    numbers["Eweight"] = 8;
    numbers["#model"] = 10;
    numbers["#search"] = 11;
    numbers["outdir"] = 12;
    numbers["end"] = 0;

    vector<string> vargv;
    string line;
    string tmp_str;
    int tmp_int;
    vector<double> tmp_wE;
    vector<int> tmp_idisco;
    char scheck[300];
    while (getline(f, line)) {
        cout << line<<endl;
        vargv.clear();
        Split(line, vargv, " ");
        tmp_str = vargv[0];
        snprintf(scheck, tmp_str.size()+1, "%s", tmp_str.c_str());
        //cout<<scheck<<endl;
        switch (numbers[scheck]) {
            case 1: //model information(groups)
                if (vargv.size() != 3) {
                    cout << "wrong model information!" << endl;
                    return 0;
                }
                n_group = atoi(vargv[1].c_str());
                fgroupnm = vargv[2];
                break;
            case 2: //para information
                if (vargv.size() != 2) {
                    cout << "wrong in.para information!" << endl;
                    return 0;
                }
                finparanm = vargv[1];
                break;
            case 3: //sw data
                if(vargv.size()<3){
                    cout << "wrong disp information!" << endl;
                    return 0;
                }
                tmp_str = vargv[1]; //R or L
                tmp_int = atoi(vargv[2].c_str());
                if(tmp_str=="R" or tmp_str=="r"){
                    for(int k=0;k<tmp_int;k++){
                        sw_tp.push_back('R');
                        assert(vargv[3+k*2].size()==1);
                        disp_tp.push_back(vargv[3+k*2][0]);
                        fdispnm.push_back(vargv[4+k*2]);
                    }
                }
                else if(tmp_str=="L" or tmp_str=="l"){
                    for(int k=0;k<tmp_int;k++){
                        sw_tp.push_back('L');
                        assert(vargv[3+k*2].size()==1);
                        disp_tp.push_back(vargv[3+k*2][0]);
                        fdispnm.push_back(vargv[4+k*2]);
                    }
                }
                break;
            case 4: //rf data
                if(vargv.size()!=4){
                    cout << "wrong rf information!" << endl;
                    return 0;
                }
                a.push_back(atof(vargv[1].c_str()));
                rayp.push_back(atof(vargv[2].c_str()));
                frfnm.push_back(vargv[3]);
                break;
            case 5: //weighting of rf data
                if(vargv.size()!=2){
                    cout << "wrong rf weighting information!" << endl;
                    return 0;
                }
                w_rf = atof(vargv[1].c_str());
                break;
            case 6: //hk data
                if(vargv.size()<4){
                    cout << "wrong hk information!" << endl;
                    return 0;
                }
                fhknm.push_back(vargv[1]);
                tmp_int = atoi(vargv[2].c_str());
				tmp_idisco.clear();
                for(int k=0;k<tmp_int;k++) {
                    tmp_idisco.push_back(atoi(vargv[k+3].c_str()));
                }
                idisco.push_back(tmp_idisco);
                break;
            case 7: // weighting of different phases in hk
                if(vargv.size()==1+2*idisco.size()){
                    for(int k=0;k<idisco.size();k++) {
                        wPs.push_back(atof(vargv[1+k*2].c_str()));
                        wPpPs.push_back(atof(vargv[2+k*2].c_str()));
                        wPsPs.push_back(1. - atof(vargv[1+k*2].c_str()) - atof(vargv[2+k*2].c_str()));
                    }
                }
                else if(vargv.size()==1+3*idisco.size()){
                    for(int k=0;k<idisco.size();k++) {
                        wPs.push_back(atof(vargv[1+k*3].c_str()));
                        wPpPs.push_back(atof(vargv[2+k*3].c_str()));
                        wPsPs.push_back(atof(vargv[3+k*3].c_str()));
                    }
                }
                else{
                    cout << "wrong hk weighting information!" << endl;
                    return 0;
                }
                break;
            case 8:// weighting of E from different discontinuity
                if(vargv.size()!=2+idisco.size()){
                    cout << "wrong E weighting information!" << endl;
                    return 0;
                }
				tmp_wE.clear();
                for(int k=0;k<idisco.size();k++){
					for(int kk=0;kk<idisco[k].size();kk++){
                    		tmp_wE.push_back(atof(vargv[1+k].c_str()));
					}
                
					wE.push_back(tmp_wE);
					tmp_wE.clear();
				}
                refE = atof(vargv[idisco.size()+1].c_str());
                cout<<"read&write:"<<refE<<endl;
                break;
            case 9: //monotonically increasing layer
                for(int k=1;k<vargv.size();k++){
                    idx_monotlyr.push_back(atoi(vargv[k].c_str()));
                }
                break;
            case 10:
                if(vargv.size()!=2){
                    cout << "wrong #model information!" << endl;
                    return 0;
                }
                n_model = atoi(vargv[1].c_str());
                break;
            case 11:
                if(vargv.size()!=2){
                    cout << "wrong #search information!" << endl;
                    return 0;
                }
                n_search = atoi(vargv[1].c_str());
                break;
            case 12:
                if(vargv.size()!=3){
                    cout << "wrong #model information!" << endl;
                    return 0;
                }
                outdir.push_back(vargv[1]);
                outdir.push_back(vargv[2]);
                break;
            case 0:
                cout<<"input is over!"<<endl;
                break;
            default:
                cout<<"I don't know what it is, exit"<<endl;
                return 0;
        }
    }
    f.close();
    return 1;
}

int read2(const vector<string> &fhklst,vector<vector<string> > &fhknm){
    ifstream f;
    string line;
    vector<string> tmp;
    for(const auto & i : fhklst) {
        cout<<i<<endl;
        f.open(i);
        if (not f.is_open()) {
            cout << "#####hk.lst file does not exist!\n" << i << endl;
            exit(0);
        }
	tmp.clear();
        while (getline(f, line)) {
            tmp.push_back(line);
        }
        fhknm.push_back(tmp);
        f.close();
    }
    return 1;
}

