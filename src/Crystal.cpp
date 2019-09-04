#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include<iomanip>
#include"Crystal.h"


Crystal::Crystal(const char *fname){
    std::string line;
    std::ifstream params (fname);
    while(std::getline(params,line)){
        if((line.at(0)=='#')) continue;
        else{
            std::stringstream ss(line);
            ss >> nCrystals >> cStar >> gamma >> c0 >> k; 
            
            c0 *= cStar;
            mu = new double[nCrystals];
            xStar = new double[nCrystals]; 

            for(int i=0;i<nCrystals;i++){
                ss >> mu[i] >> xStar[i];
            }
        }
    }
}

Crystal::~Crystal(){
    if(this->nCrystals > 1){
        delete[] mu;
        delete[] xStar;
    }
}

Crystal& Crystal::operator=(const Crystal& that){
   nCrystals = that.nCrystals;        
   gamma = that.gamma;
   c0 = that.c0;
   cStar = that.cStar;
   k = that.k;
   double *localMu = new double[that.nCrystals];
   double *localxStar = new double[that.nCrystals];
   for(int i=0; i<that.nCrystals; i++){
       localMu[i] = that.mu[i];
       localxStar[i] = that.xStar[i];
   }
   delete[] mu;
   delete[] xStar;
   for(int i=0; i<that.nCrystals; i++){
       mu[i] = localMu[i];
       xStar[i] = localxStar[i];
    }
    return *this;
}

void Crystal::fcn(double t, double *x){
    double c = c0;

    for(int i=0; i<nCrystals; i++){
        c += mu[i]*xStar[i]*xStar[i]*xStar[i] - mu[i]*x[i]*x[i]*x[i];
    }

    for(int i=0; i<nCrystals; i++){
        x[i] = k*(c - cStar*std::exp(gamma/x[i]));
    }

    return; 
}

int Crystal::test(double *x){
    for(int i=0; i<nCrystals; i++){
        if(x[i] < 0) return(-1);
        else return(0);    
    }
}

void Crystal::createParams(const char *fname, double Yf){
    std::ofstream params(fname, std::ios::trunc);
    params << std::setprecision(9);
    params << std::fixed;
    int neq = 1;
    params << neq << '\t' << this->cStar << '\t' << this->gamma 
            << '\t' << (this->c0)/(this->cStar) << '\t' << this->k 
            << '\t' << this->mu[1] << '\t' <<  Yf;
    params.close();
}
