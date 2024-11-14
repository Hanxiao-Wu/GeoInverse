//
// Created by Hanxiao Wu on 8/14/23.
//

#ifndef MCMC_CPLUS_READ_WRITE_H
#define MCMC_CPLUS_READ_WRITE_H
#include "MC.h"
#include <map>
int read1(const string& f_control, string& fgroupnm, string &finparanm, int &n_group, vector<string> &fdispnm, vector<char> &sw_tp, vector<char> &disp_tp, vector<string> &frfnm,vector<double> & a,vector<double> &rayp, vector<string> &fhknm,vector<int> & idisco, double &w_rf, vector<double> &wPs,vector<double> &wPsPs,vector<double> &wPpPs, vector<vector<double>> &wE, vector<int> &idx_monotlyr, int &n_model, int &n_search, vector<string> &outdir, double &refE);
int read2(const vector<string> &fhklst,vector<vector<string> > &fhknm);

#endif //MCMC_CPLUS_READ_WRITE_H
