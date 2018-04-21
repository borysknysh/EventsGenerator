#ifndef GKINEMATICS_H
#define GKINEMATICS_H

#include "GInterfaceForMathFunctions.h"

namespace Gamapola
{
  class GKinematics: virtual public GInterfaceForMathFunctions{
  public:
    double GSjkMin(const double& sij, const double& mi, const double& mj, const double& mk, const double& M);
    double GSjkMax(const double& sij, const double& mi, const double& mj, const double& mk, const double& M);
  private:
    int fArg;
    std::complex<double> fE0PiPj;
    std::complex<double> fE0Pi;
    std::complex<double> fE0Pj;
  protected:
    GKinematics();
    virtual ~GKinematics();
    
    auto GMomentaA(double s, double sij, double massMes)->decltype(s);
    auto GPiPjVec4(double s, double sij, double sjk)->decltype(s);
    auto GEnergyA(double s, double sij, double massMes)->decltype(s);
    auto GSij(double s, double sij, double sik)->decltype(s);
    
    void GPiPjVec3();
    std::complex<double> GScPrd(const double* pi, const std::complex<double>* pj);
    double GPiPj();
    void GCalculateKinematics(double* p);
    double GExponentialM(double *p, int nRes, int nFnc);
    double GCosThetaij(const double& cosThstar, const double& sij, const double& mi, 
                       const double& mj, const double& mk, const double& M);
    void GKinematicsConstanstsCalculation(int nRes);
    
    std::vector<std::vector<std::vector<double> > > fConstParsForMomA;
  };
}
#endif