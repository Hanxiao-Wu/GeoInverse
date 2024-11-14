//
// Created by Hanxiao Wu on 8/9/23.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "MC.cpp"
#include "function.cpp"
#include "SAC_IO/sac_stream.cpp"
#include "SAC_IO/sac_io.cpp"
#include <string>
using namespace std;
/*
extern"C"
{
void fast_surf_(int *n_layer0,int *kind0,float *a_ref0,float *b_ref0,float *rho_ref0,float *d_ref0,float *qs_ref0,float *cvper,int *ncvper,float *uR0,float *uL0,float *cR0,float *cL0,float *rR0, float *rL0, float *ampl0);
void theo_(int *n,float *fbeta,float *h,float *vps,float *qa,float *qb,float *fs,float *din,float *a0,float *c0,float *t0,int *nd,float *rx);
}
*/
int main(){
    int newnlayer,nn,i,k,nper,nper1,N=200;
    float *tvs, *tvp, *tvpvs, *tqs1, *tqs2, *tqp, *trho, *tthick; // can use a 2d array for simplicity. but to make code clear, use 1d here.
    float *uR0,*uL0,*cR0,*cL0,*rR0,*rL0,*period,*ampl0;
    double ttt;
    float temp_amp;
    ifstream fin;
    vector<double>::iterator id;
    vector<double> tvel;
    tvs = new float[N]; tvp = new float[N]; tvpvs = new float [N];
    tqs1= new float[N]; tqs2= new float[N]; tqp = new float [N];
    trho= new float[N]; tthick= new float[N];
    uR0 = new float[N]; uL0 = new float[N]; cR0= new float[N]; cL0= new float[N];
    period = new float[N];
    rR0 = new float[N]; rL0 = new float[N]; ampl0 = new float[N];
    for(i=0;i<N;i++)
    {
        *(tvs+i)=*(tvp+i)=*(tqs1+i)=*(tqs2+i)=*(trho+i)=*(tthick+i)=0.;
        *(uR0+i)=*(uL0+i)=*(cR0+i)=*(cL0+i)=*(period+i)=0.;
        *(rR0+i)=*(rL0+i)=0.;
    }
    cout<<"initialization:done"<<endl;
    nn=0;
    //cout<<"#layer is :";
    //cin>>newnlayer;
    cin>>newnlayer;
    fin.open("./test/MC.test.lyr");
    float tmp;
    for(i=0;i<newnlayer;i++){
        for(int j=0;j<7;j++){
            switch (j) {
                case 0:
                    fin >> tmp;break;
                case 1:
                    fin >> tvs[i];break;
                case 2:
                    fin >> tvp[i];break;
                case 3:
                    fin >> trho[i];break;
                case 4:
                    fin >> tqs1[i];
                    tqs2[i] = 1./tqs1[i];
                    break;
                case 5:
                    fin >> tqp[i];break;
                case 6:
                    fin >> tthick[i];break;
            }
        }
        //cout<<"Vs is:";
        //cin>>tvs[i];
        //trho[i] = 1.22679+1.53201*tvs[i]-0.83668*pow(tvs[i],2)+0.20673*pow(tvs[i],3)-0.01656*pow(tvs[i],4);
        //tvp[i] = 0.9409 + 2.0947*tvs[i]-0.8206*pow(tvs[i],2)+0.2683*pow(tvs[i],3)-0.0251*pow(tvs[i],4);
        cout<<tthick[i]<<" "<<tvs[i]<<" "<<tvp[i]<<" "<<trho[i]<<" "<<tqs1[i]<<" "<<tqs2[i]<<endl;
        tvpvs[i] = tvp[i]/tvs[i];
        //cout<<"thk is:";
        //cin>>tthick[i];
        //cout<<"Qs is:";
        //cin>>tqs[i];
        //tqs[i] = 1./tqs[i];
        //cout<<"Qp is:";
        //cin>>tqp[i];
        //tqp[i] = 1./tqp[i];
        nn=nn+1;
    }//for i
    int nl = nn;
    cout<<"#period: ";
    cin>>nper;
    for(i=0;i<nper;i++){
        cout<<"period: ";
        cin>>period[i];
    }
    tvs[nn]=tvs[nn-1]+(float)0.01;
    tvp[nn]=tvp[nn-1]+(float)0.01;
    tqs1[nn]=tqs1[nn-1];//tqs[nn-2];
    tqs2[nn]=tqs2[nn-1];
    trho[nn]=trho[nn-1]+(float)0.01;
    tthick[nn]=0.;
    cout<<tthick[nn]<<" "<<tvs[nn]<<" "<<tvp[nn]<<" "<<trho[nn]<<" "<<tqs1[nn]<<" "<<tqs2[nn]<<endl;
    nn=nn+1;//#of element inside the arrays

    int cflag=2;
    cout<<"done:model"<<endl;
    cout<<"nn:"<<nn<<endl;
    cout<<"period:"<<nper<<" "<<period[0]<<endl;
    fast_surf_(&nn,&cflag,tvp,tvs,trho,tthick,tqs2,period,&nper,uR0,uL0,cR0,cL0,rR0,rL0,ampl0);
    for(i=0;i<nper;i++){
        cout<<period[i]<<" "<<cR0[i]<<endl;
    }

    float slow,pi,din,rt;
    int NN=2000;
    float *rx;
    rx=new float[NN];
    nn=150;
    slow=0.06;
    pi=atan(1.)*4.;
    //cout<<"===== pi="<<pi<<endl;
    din=180.*asin(tvs[nl-1]*tvpvs[nl-1]*slow)/pi;
    cout<<"din:"<<din<<endl;
    float gau = 2.5;
    float dt = 0.005; // ???parameter of a minimum amplitude level
    float t0 = 0.;
    int nn1 = (int)newnlayer;
    rt = (int)1./0.1;
    theo_(&nn1,tvs,tthick,tvpvs,tqp,tqs1,&rt,&din,&gau,&dt,&t0,&nn,rx);
    for(i=0;i<nn;i++){cout<<i*1./rt<<" "<<rx[i]<<endl;}
    return 0;
}