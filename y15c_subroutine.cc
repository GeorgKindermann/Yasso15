#include <cmath>
#include "y15c_subroutine.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

  struct complex {double realnum, imagnum;};

  extern"C" {
    void dgchbv_(int *m, double *t, double *H, int *ldh, double *y, complex *wsp, int *iwsp, int *iflag);
  }
  
  //https://github.com/ngocson2vn/MPM/blob/master/547_546-parallel-GIMP/0main/fit.cpp
  //http://www.sci.utah.edu/~wallstedt/LU.htm
  //Philip Wallstedt
  //Decompose
  template<const int d>inline void Crout(const std::array<double,d*d>& a, std::array<double,d*d>& A){
    for(int k=0;k<d;++k){
      for(int i=k;i<d;++i){
	double sum=0.;
	for(int p=0;p<k;++p)sum+=A[i*d+p]*A[p*d+k];
	A[i*d+k]=a[i*d+k]-sum;
      }
      for(int j=k+1;j<d;++j){
	double sum=0.;
	for(int p=0;p<k;++p)sum+=A[k*d+p]*A[p*d+j];
	A[k*d+j]=(a[k*d+j]-sum)/A[k*d+k];
      }
    }
  }
  template<const int d>inline void solveCrout(const std::array<double,d*d>& A, const std::array<double,d>& b, std::array<double,d>& x){
    std::array<double,d> y;
    for(int i=0;i<d;++i){
      double sum=0.;
      for(int k=0;k<i;++k)sum+=A[i*d+k]*y[k];
      y[i]=(b[i]-sum)/A[i*d+i];
    }
    for(int i=d-1;i>=0;--i){
      double sum=0.;
      for(int k=i+1;k<d;++k)sum+=A[i*d+k]*x[k];
      x[i]=(y[i]-sum);
    }
  }

  //Matrixmultiplication
  template<int N, int M, int P> void MATMUL(const std::array<double,N*M> &A, const std::array<double,M*P> &B, std::array<double,N*P> &C) {
    if(&B == &C) {
      std::array<double,N*P> RES{0.};
      for (int n=0; n<N; ++n) {
	for (int p=0; p<P; ++p) {
	  for (int m=0; m<M; ++m) {
	    RES[n*P+p] += A[n*M+m] * B[m*P+p];
	  }
	}
      }
      for (int i=0; i<N*P; ++i) {C[i] = RES[i];}
    } else {
      C.fill(0.);
      for (int n=0; n<N; ++n) {
	for (int p=0; p<P; ++p) {
	  for (int m=0; m<M; ++m) {
	    C[n*P+p] += A[n*M+m] * B[m*P+p];
	  }
	}
      }
    }
  }

  //Functions for solving the diff. equation, adapted for the Yasso case
  void matrixnorm(const std::array<double, 5*5>& A, double& b) {
    //returns elementwise (i.e. Frobenius) norm of a square matrix
    const int n=5;
    b = 0.;
    for(int i=0; i<n*n; ++i) {b += A[i]*A[i];}
    b = sqrt(b);
  }
  void matrixexp(const std::array<double, 5*5>& A, std::array<double, 5*5>& B, const size_t& q) {
    //Approximated matrix exponential using Taylor series with scaling & squaring
    //Accurate enough for the Yasso case
    const int n=5;
    //int q = 11; // #terms in Taylor  gk: possibility to reduce to speed up
    //expect B is initialized:  for(int i=0; i<25; ++i) {B[i] = 0.;}
    for(int i=0; i<n; ++i) {B[i*6] = 1.;}
    double normiter = 2.; // Amount of scaling & squaring
    int j=2;
    double p;
    matrixnorm(A, p);
    while(p>normiter) {
      normiter *= 2.;
      ++j;
    }
    std::array<double, 5*5> C,D;
    for(int i=0; i<25; ++i) {C[i] = A[i]/normiter;} //scale
    for(int i=0; i<25; ++i) {B[i] += C[i];}
    for(int i=0; i<25; ++i) {D[i] = C[i];}
    for(size_t i=2; i<q; ++i) { //compute Taylor expansion
      MATMUL<5,5,5>(C,D,D);
      for(int j=0; j<25; ++j) {D[j] /= double(i);}
      for(int i=0; i<25; ++i) {B[i] += D[i];}
    }
    for(int i=1; i<j; ++i) { //square
      MATMUL<5,5,5>(B,B,B);
    }
  }

  //https://github.com/blackrim/lagrange/tree/master/src
  //Krylov Chebyshev matrix exp
  //https://www.maths.uq.edu.au/expokit/
  //Compute exp(A * t) * v using Chebyshev
  void dgchbv(std::array<double, 5*5>& A, const std::array<double, 5>& v, std::array<double, 5>& e) {
    int m{5}, ldh{5}, iflag, iwsp[5];
    double t{1.};
    complex wsp[2*5*(5+2)];
    for(int i=0; i<5; ++i) {e[i] = v[i];}
    dgchbv_(&m, &t, A.data(), &ldh, e.data(), wsp, iwsp, &iflag);
  }

}

namespace yasso {

  yasso15::yasso15() : A{0.} {
    setTheta({4.8971473e-01,4.9138734e+00,2.4197346e-01,9.4876416e-02,4.3628932e-01,2.4997402e-01,9.1512685e-01,9.9258227e-01,8.3853738e-02,1.1476783e-02,6.0831497e-04,4.7612821e-04,6.6037729e-02,7.7134168e-04,1.0401742e-01,6.4880756e-01  -1.5487177e-01  -1.9568024e-02  -9.1717130e-01  -4.0359430e-04  -1.6707272e-04,9.0598047e-02  -2.1440956e-04,4.8772465e-02  -7.9136021e-05,3.5185492e-02  -2.0899057e-04  -1.8089202e+00  -1.1725473e+00  -1.2535951e+01,4.5964720e-03,1.3025826e-03  -4.3892271e-01,1.2674668e+00,2.5691424e-01});
  }
  
  yasso15::yasso15(const std::array<double, 35> &atheta) : A{0.} {
    setTheta(atheta);
  }
  
  void yasso15::setTheta(const std::array<double, 35> &atheta) {
    for(int i=0; i<35; ++i) {theta[i] = atheta[i];}
    theta[31] = -std::fabs(theta[31]);
    theta[34] = -std::fabs(theta[34]);
    for(int i=0; i<4; ++i) {theta[i] = -std::fabs(theta[i]);}
    aIsSet = false;
  }

  void yasso15::setClimSizeLeach(const double& avgT, const double& sumP, const double& ampT, const double& diam, const double& leach) {
    double tem{0.};
    double temN{0.};
    double temH{0.};
    const double m3 = sumP/1000.;
    {
      const double pi=M_PI;
      //temperature annual cycle approximation
      const double m0 = (1./sqrt(2.)-1.)/pi;
      const double m1 = 1./sqrt(2.)/pi;
      const double m2 = (1.-1./sqrt(2.))/pi;
      const double clim = 4. * ampT;
      double te[4];
      te[0] = avgT + clim * m0;
      te[1] = avgT - clim * m1;
      te[2] = avgT + clim * m2;
      te[3] = avgT + clim * m1;
      
      //Average temperature dependence
      double te2[4];
      for(int i=0; i<4; ++i) {te2[i] = te[i] * te[i];}
      for(int i=0; i<4; ++i) {
	tem += exp(theta[21]*te[i]+theta[22]*te2[i]);
	temN += exp(theta[23]*te[i]+theta[24]*te2[i]);
	temH += exp(theta[25]*te[i]+theta[26]*te2[i]);
      }
      tem /= 4.; temN /= 4.; temH /= 4.;
      
      //Precipitation dependence
      tem *= 1.-exp(theta[27] * m3);
      temN *= 1.-exp(theta[28] * m3);
      temH *= 1.-exp(theta[29] * m3);
    }
    
    //! check rare case where no decomposition happens for some compartments 
    //! (basically, if no rain)
    //gk: Changed that there is also a solution for spinnup
    if(tem < tol) {
      noDecomposition = true;
      return;
    }
    noDecomposition = false;
    
    //Size class dependence -- no effect if d == 0.0
    if(diam > 0.) { //gk: orig calculates also for negative diam
      const double size_dep = std::min(1., pow(1. + theta[32]*diam + theta[33]*diam*diam, theta[34]));
      //Calculating matrix a (will work ok despite the sign of alphas)
      for(int i=0; i<3; ++i) {A[i*6] = theta[i]*tem*size_dep;}
      A[3*6] = theta[3]*temN*size_dep;
    } else {
      for(int i=0; i<3; ++i) {A[i*6] = theta[i]*tem;}
      A[3*6] = theta[3]*temN;
    }
    A[24] = theta[31]*temH; //no size effect in humus
    const std::array<double,4> dAbs{std::fabs(A[0]), std::fabs(A[6]), std::fabs(A[2*6]), std::fabs(A[3*6])};
    for(int i=0, idx=3; i<4; ++i) {
      for(int j=0; j<4; ++j) {
	if(i!=j) {A[i*5+j] = theta[++idx] * dAbs[j];}
      }
      //gk: default 0 init:  A[i*5+4] = 0.; //no mass flows from H -> AWEN
    }
    //mass flows AWEN -> H (size effect is present here)
    for(int i=0; i<4; ++i) {A[20+i] = theta[30] * dAbs[i];}
    //Leaching (no leaching for humus) 
    if(leach < 0.) { //gk: orig calculates also for positive leach
      const double aux = leach * m3;
      for(int i=0; i<4; ++i) {A[6*i] += aux;}
    }
    aIsSet=true;
    aIsDecomp=false;
  }

  void yasso15::getSpin(const std::array<double, 5>& infall, std::array<double, 5>& result) {
    if(noDecomposition) {
      for(int i=0; i<5; ++i) {
	if(infall[i] > 0.) {result[i] = std::numeric_limits<double>::infinity();
	} else {result[i] = 0.;}
      }
      return;
    }
    if(aIsSet) {
      if(!aIsDecomp) {
	Crout<5>(A, Adecomp);
	aIsDecomp = true;
      }
      solveCrout<5>(Adecomp, infall, result);
      for(int i=0; i<5; ++i) {result[i] *= -1.;}
    } else {
      for(int i=0; i<5; ++i) {
	result[i] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  void yasso15::getNextTimestep(const std::array<double, 5>& init, const std::array<double, 5>& infall, std::array<double, 5>& result, const double timespan, const int fun) {
    if(noDecomposition) {
      for(int i=0; i<5; ++i) {result[i] += infall[i];}
      return;
    }
    if(aIsSet) {
      //Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init
      //Solve DE in given time
      //Solve Matrix Differential equation Taylor x'(t) = Ax(t) + b
      MATMUL<5,5,1>(A,init,z1);
      for(int i=0; i<5; ++i) {z1[i] += infall[i];}
      for(int i=0; i<5*5; ++i) {At[i] = A[i] * timespan;}
      switch(fun) {
      case 1: //Using Expotit
	dgchbv(At, z1, z2);
	break;
      default: //Use default function from Yasso15
	mexpAt.fill(0.);
	matrixexp(At, mexpAt, taylorTerms);
	MATMUL<5,5,1>(mexpAt,z1,z2);
      }
      for(int i=0; i<5; ++i) {z2[i] -= infall[i];}
      if(!aIsDecomp) {
	Crout<5>(A, Adecomp);
	aIsDecomp = true;
      }
      solveCrout<5>(Adecomp, z2, result);
    } else {
      for(int i=0; i<5; ++i) {
	result[i] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

}
