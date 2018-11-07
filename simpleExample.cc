#include <iostream>
#include <array>

#include "y15c_subroutine.h"

using namespace std;

extern"C" {
  void __yasso_MOD_mod5c(float *theta, float *time, float *climate, float *init, float *b, float *d, float *leach, float *xt, int *steadystate_pred);
  void __yasso07_MOD_mod5c(float *theta, float *time, float *climate, float *init, float *b, float *d, float *leach, float *xt);
}

int main() {
  { //Original Yasso07
    float theta[27] = {-0.49,-4.9,-0.24,-0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,0.091,-0.00021,-1.8,-0.15,-0.02,-0.92,-0.44,1.3,-0.26,-0.0013,0.0046};
    float time {1.}; //Time to run
    float climate[3] {10., 600., 12.}; //Temp annual average [C], precip annual summ [mm], amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    float init[5] {0.,0.,0.,0.,0.}; //Initial State
    float b[5] {0.5,0.1,0.1,0.2,0.}; //Infall
    float xt[5] {0.,0.,0.,0.,0.}; //Result
    float d {2.}; //size [cm]
    float leach {0.}; //Leaching

    cout << "Origin Yasso07" << endl;
    for(int year=0; year<10; ++year) {
      __yasso07_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt);
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << xt[i];}
      cout << endl;
      for(int i=0; i<5; ++i) {init[i] = xt[i];}
    }
    time = 9999.;
    __yasso07_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << xt[i];}
    cout << endl;

  }
  
  { //Original Yasso15
    float theta[35] = {0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26}; //Parameters
    float time {1.}; //Time to run
    float climate[3] {10., 600., 12.}; //Temp annual average [C], precip annual summ [mm], amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    float init[5] {0.,0.,0.,0.,0.}; //Initial State
    float b[5] {0.5,0.1,0.1,0.2,0.}; //Infall
    float xt[5] {0.,0.,0.,0.,0.}; //Result
    float d {2.}; //size [cm]
    float leach {0.}; //Leaching
    int steadystate_pred {0}; //set to true if ignore 'time' and compute solution in steady-state conditions (which sould give equal solution as if time is set large enough)
    
    cout << "Origin Yasso15" << endl;
    for(int year=0; year<10; ++year) {
      __yasso_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt, &steadystate_pred);
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << xt[i];}
      cout << endl;
      for(int i=0; i<5; ++i) {init[i] = xt[i];}
    }
    steadystate_pred = 1;
    __yasso_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt, &steadystate_pred);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << xt[i];}
    cout << endl;
  }

  { //Rewritten Yasso15
    array<double, 35> theta = {0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26}; //Parameters
    double time {1.}; //Time to run
    double avgT = 10.; //Temp annual average [C]
    double sumP = 600.; //Precip annual summ [mm]
    double ampT = 12.; //Amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    array<double, 5> init {0.,0.,0.,0.,0.}; //Initial State
    array<double, 5> infall {0.5,0.1,0.1,0.2,0.}; //Infall
    array<double, 5> result;
    double diam {2.}; //size [cm]
    double leach {0.}; //Leaching

    yasso::yasso15 yasso(theta);
    
    cout << "\nRewritten Yasso15" << endl;
    yasso.setClimSizeLeach(avgT, sumP, ampT, diam, leach);
    yasso.setTimespan(time);
    for(int year=0; year<10; ++year) {
      yasso.getNextTimestep(init, infall, init, 0); //Write result in init
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << init[i];}
      cout << endl;
    }
    yasso.getSpin(infall, result);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << result[i];}
    cout << endl;
  }
    
  return(0);
}
