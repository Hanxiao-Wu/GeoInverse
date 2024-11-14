//
// Created by Hanxiao Wu on 8/8/23.
//

#ifndef MCMC_CPLUS_FORWARVAR_H
#define MCMC_CPLUS_FORWARVAR_H
namespace Forward {
    const int LIN = 5, LOT = 6, NL = 200, NL2 = NL + NL;
    const int NPERIOD=2049;
    //double fval[NPERIOD];
    //int nfval;
    int mmax, mode;
    vector<int> iwat;
    double twopi,displ,dispr;
    double mvts[2]
    vector<vector<double> > vts;
    vector<double> od,oa,ob,orho,qa,qb,etap,etas,frefp,frefs;
    double refdep;
    bool allfluid;
}
#endif //MCMC_CPLUS_FORWARVAR_H
