#ifndef GCURRENTS_H
#define GCURRENTS_H

#include "GInterfaceForMathFunctions.h"

namespace Gamapola
{
  class GCurrents: virtual public GInterfaceForMathFunctions{
  private:
    int fArg;
    static constexpr double fK = 1.;
    static constexpr double fRho = 1.;
    std::complex<double> cC[2][fMaxNRes][3][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro, Kappa}; {C1,C2}; {iPhS,iPhD}
  protected:
     GCurrents();
    virtual ~GCurrents();
    void GCalculateHelicityAmplitudes(double* x,double* p);
    void GCalculateCouplings(double*x, double* p);
    
    std::complex<double> GC1Kst(double* p, double s, double sij, double sjk, int nRes);
    std::complex<double> GC2Kst(double* p, double s, double sij, double sjk, int nRes);
    std::complex<double> GC1Rho(double* p, double s, int nRes);
    std::complex<double> GC2Rho(double* p, double s, int nRes);
    std::complex<double> GC1Kappa(double* p, double sjk, int nRes);
    std::complex<double> GC2Kappa(double* p, double sij, int nRes);
    
    std::complex<double> GC3Kst(double s, double sjk, int nRes);
    std::complex<double> GC3Rho(double s, int nRes);
    std::complex<double> GC4Kst(double s, double sij, int nRes);
    std::complex<double> GC5Kst(double s, double sjk, int nRes);
    std::complex<double> GC4Rho(double s, int nRes);
    
    std::complex<double> GC1(double* p, double s, double sij, double sjk, int nRes);
    std::complex<double> GC2(double* p, double s, double sij, double sjk, int nRes);
    std::complex<double> GC3(double s, double sjk, int nRes);
    std::complex<double> GC4(double s, int nRes);
    std::complex<double> GC5(double s, int nRes);
    
    double GJ2(double* p, double lambda, double s);
    std::complex<double>* GJVec(double* p, int charge, int nRes);
    std::complex<double>* GLVec(int charge, int nRes);
    std::complex<double>* GKVec(int charge, int nRes);
    
    double fLeft[fMaxNRes];
    double fRight[fMaxNRes];
    std::complex<double> fc1[2][fMaxNCharges];
    std::complex<double> fc2[2][fMaxNCharges];
    std::complex<double> fc3[2][fMaxNCharges];
    std::complex<double> fc4[2][fMaxNCharges];
    std::complex<double> fc5[2][fMaxNCharges];
    std::complex<double> fc4Rho;
    std::complex<double> fc4Kstr;
    std::complex<double> fc5Kstr;
    std::complex<double> fJVec3[fMaxNRes][fMaxNCharges][3];  
    double fDelta[fMaxNRes];
    std::complex<double> fDelta2[fMaxNRes];
    std::complex<double> fCCouplings[fMaxNRes][fMaxNDecays];
    double fFF[fMaxNRes][fMaxNDecays];
    double fK1[2];
    int fCharge[fMaxNRes];
    double fOmega;
    std::complex<double> cRL[2][fMaxNRes][3][2][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro, Kappa}; {C1,C2}; {R,L}; {iPhS,iPhD}
  };
}
#endif