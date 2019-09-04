#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include"Crystal.h"
#include"RK4.h"

int main(){

//create rk4 object with the Crystal type
//initialization will get parameters from
//specified data file
    char* fRK = "params1.dat";
    char* fOst = "ostwald.dat";
    double* Yf;

    try{
        //Create RK4 and Crystal objects, which will initialize data
        RK4<Crystal> rk4Objn(fRK);
        Crystal c1(fOst);

        //Run the RK4 calculation for both crystals
        Yf = rk4Objn.calc(c1,"out1.dat");
        
        //Set file names for new sets of parameters
        fRK = "params2.dat";
        fOst = "ostwald2.dat";

        //Create new files for parameters
        //for the surviving crystal
        c1.createParams(fOst, Yf[1]);
        rk4Objn.createParams(fRK, Yf[1], 3.0, 10000);

        //Create new objects based on these parameters
        Crystal c2(fOst);
        RK4<Crystal> rk4Obj1(fRK);

        Yf = rk4Obj1.calc(c2,"out2.dat");
        std::cout << "back from calc2.\n";
        delete[] Yf;
    }
    catch(std::string err){
        std::cout << "Exception: " << err << std::endl;
    }

	return(0);
}
