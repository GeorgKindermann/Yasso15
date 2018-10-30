#ifndef YASSO_H
#define YASSO_H

//Rewritten Yasso15 core code based on Fortran code from Mikko Tuomi and Marko Järvenpää

#include <array>
#include <cmath>
#include <limits>

namespace yasso {

  class yasso15 {
  public:
    yasso15();
    yasso15(const std::array<double, 35> &theta);
    void setTheta(const std::array<double, 35> &theta);
    void setClimSizeLeach(const double& avgT, const double& sumP, const double& ampT, const double& diam, const double& leach);
    bool isThereDecomposition() {return(!noDecomposition);}
    void getSpin(const std::array<double, 5>& infall, std::array<double, 5>& result);
    void getNextTimestep(const std::array<double, 5>& init, const std::array<double, 5>& infall, std::array<double, 5>& result, const double timespan=1., const int fun=0);
    size_t setTaylorTerms(const size_t& n) {taylorTerms = n+1; return(taylorTerms-1);}
  private:
    std::array<double, 35> theta;
    std::array<double, 5*5> A;
    bool aIsSet=false;
    std::array<double, 5*5> Adecomp;
    bool aIsDecomp=false;
    double tol = 1.E-12;
    bool noDecomposition = true;
    std::array<double, 5> z1, z2;
    std::array<double, 5*5> At;
    std::array<double, 5*5> mexpAt;
    size_t taylorTerms = 11;
  };
  
}

#endif
