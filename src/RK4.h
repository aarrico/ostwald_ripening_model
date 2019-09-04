template<class T> 
class RK4{
    public:
        RK4(const char *);
        ~RK4();
	    double* calc(T &,const char *);
        void createParams(const char *, double, double, int);
    private:
        int neq, steps;
        double t0, tEnd, *Yn; 
};

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>
#include<omp.h>

template<class T>
inline
RK4<T>::RK4(const char *fname){
    std::string line;

//make sure file opened correctly
    std::ifstream params (fname);	
    if(!params.is_open()){	
        throw "params.dat did not open";
    }

   while(std::getline(params,line)){
        if((line.at(0)=='#')) continue;
        else{
            std::stringstream ss(line);
            ss >> neq >> t0 >> tEnd >> steps; 

            Yn = new double[neq];   //allocate arrays based on neq

            for(int i=0;i<neq;i++){
                ss >> Yn[i];     //read in initial Y for each eqn
            }
        }
    }
    params.close();
}

template <class T>
inline
RK4<T>::~RK4(){
    if(this->neq > 1){
        delete[] Yn;
    }
}

template <class T>
inline
void RK4<T>::createParams(const char *fname, double Yf, double tEnd, int steps){
    std::ofstream params(fname,std::ios::trunc);
    params << std::setprecision(16);
    params << std::fixed;
    int neq = 1;
    params << neq << '\t' << this->tEnd << '\t' << tEnd 
            << '\t' << steps << '\t' << Yf;
    params.close();
}

template <class T>
inline
double* RK4<T>::calc(T &obj, const char *fname){
    std::ofstream fout(fname,std::ios::trunc);
    if(!fout.is_open()){	
        throw "output.dat did not open";
    }
    fout << std::setprecision(16);
    fout << std::fixed;

    double dt = (tEnd-t0)/steps, dt2 = dt/2.0;
    double tn = 0.0;
    double k1[neq],k2[neq],k3[neq],k4[neq];
    int numThreads;

    for(int i=0; tn <= tEnd; i++){
        //Set tn for current step
        tn = t0 + i*dt;
        //Report current n values
        fout << i << '\t' << tn << '\t';
        for(int j=0; j<neq; j++) fout << Yn[j] << '\t';
        fout << std::endl;
        
        numThreads = omp_get_num_procs();
        omp_set_num_threads(numThreads);

        #pragma omp parallel for
        for(int j=0; j<neq; j++){
            k1[j] = Yn[j];
        }

        //Runge-Kutta stages
        //Function accepts the k values as pointers
        //So they will be returned properly     
        obj.fcn(tn,k1);        

        #pragma omp parallel for
        for(int j=0; j<neq; j++) k2[j] = Yn[j] + dt2*k1[j]; 
        obj.fcn(tn+dt2,k2);
            
        #pragma omp parallel for
        for(int j=0; j<neq; j++) k3[j] = Yn[j] + dt2*k2[j];
        obj.fcn(tn+dt2,k3);

        #pragma omp parallel for
        for(int j=0; j<neq; j++) k4[j] = Yn[j] + dt*k3[j];
        obj.fcn(tn+dt,k4);

        #pragma omp parallel for
        for(int j=0; j<neq; j++){
            Yn[j] += (dt/6.0)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);            
        }
        if(obj.test(Yn) == -1){
            fout.close();
            return Yn;
        }
    }   
 
    this->tEnd = tn;
    fout.close();

    return Yn;
}

