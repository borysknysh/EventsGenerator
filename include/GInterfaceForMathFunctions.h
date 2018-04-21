#ifndef GINTERFACEFORMATHFUNCTIONS_H
#define GINTERFACEFORMATHFUNCTIONS_H

#include <iostream>
#include <cmath>
#include <TMath.h>
#include <complex>
#include <cstring>
#include "Math/WrappedMultiTF1.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "GPhysicsPDG.h"

namespace Gamapola
{
  class GInterfaceForMathFunctions{
//     https://stackoverflow.com/questions/7304511/partial-derivatives --- will be useful=)
//     https://dou.ua/lenta/articles/it-position-ds-ml/ --- machine learning
//     virtual void GCoeffsForIS() = 0;
  protected:
    GInterfaceForMathFunctions();
    virtual ~GInterfaceForMathFunctions(); 
    virtual void GCalculateQPMC(double* p,const int& iRes, const int& swap) = 0;
    virtual void GCalculateFormFactorsKstar(double* x,double* p, int nRes, int nSwap) = 0;
    virtual void GCalculateFormFactorsRho(double* x,double* p, int nRes) = 0;
    virtual void GCalculateHelicityAmplitudes(double* x,double* p) = 0;
    virtual void GCalculateCouplings(double*x,double* p) = 0;
    virtual std::complex<double> GAKstar(double* p, int nRes, int nSwap) = 0;
    virtual std::complex<double> GBKstar(double* x,double* p, int nRes, int nSwap) = 0;
    virtual std::complex<double> GARho(double* p, int nRes, int nSwap) = 0;
    virtual std::complex<double> GBRho(double* x,double* p, int nRes, int nSwap) = 0;
    
    virtual void GCouplingConstantsCalculation() = 0;
//     Breit-Wigner
    virtual std::complex<double> GBWKstr(double sij) = 0;
    virtual std::complex<double> GBWRho(double sij) = 0;
    virtual std::complex<double> GBWKappa(double sij) = 0;
    virtual std::complex<double> GBWK1(double s, int nRes) = 0;
    virtual inline double GSjkMin(const double& sij, const double& mi, const double& mj, const double& mk, const double& M) = 0;
    virtual inline double GSjkMax(const double& sij, const double& mi, const double& mj, const double& mk, const double& M) = 0;
    virtual double GMomentaA(double s, double sij, double massMes) = 0;
    virtual double GPiPjVec4(double s, double sij, double sjk) = 0;
    virtual void GPiPjVec3() = 0;
    virtual std::complex<double> GScPrd(const double* pi, const std::complex<double>* pj) = 0;
    virtual double GEnergyA(double s, double sij, double massMes) = 0;
    virtual double GSij(double s, double sij, double sik) = 0;
    virtual double GPiPj() = 0;
    virtual void GCalculateKinematics(double* p) = 0;
    virtual double GExponentialM(double *p, int nRes, int nFnc) = 0;   
    
    virtual double GProcessingComputationOfPDF(double* x, double* p) = 0; 
    
    static const int fMaxNRes = 7;
    static const int fMaxNDecays = 3;
    static const int fMaxNCharges = 3;
    static const int fMaxNWaves = 3;
    double fMasses[4];
    const std::complex<double>* fNormalizationIntegrals;
    int fNIntegrals;
//     std::vector<double (Gamapola::GMathFunctions::*) (double* x, double* p)> fFunctionPointers; 
//     For I0, I1:
    int fIndex[fMaxNRes];
    
    std::vector<std::vector<double> >  fMomentaConst;
    double fMomentaA[fMaxNDecays];
    double fResWidths[fMaxNRes];
    double fMK0star[fMaxNRes][fMaxNDecays];
    double fMRho0[fMaxNRes][fMaxNDecays];
//     Coupling GCouplingConstantsCalculation
    double fgKstarKpi;
    double fgRhoPiPi;
    double fgK1KappaPi[2];
    double fCGKstr[fMaxNRes];
    double fCGRho[fMaxNRes];
    
//     Model parameter
//     Some kinematical variables
    double fEnergy1;
    double fEnergy2;
    double fPiPjVec4;
    double fPiPjVec[3];
    double fSij;
    double fSik;
    double fMomi;
    double fMomj;
    double fMomij;
    double fPhii;
    double fPhij;
    double fMomiV[3];
    double fMomjV[3];
    std::complex<double> fK2KinC4combo[3];
    std::complex<double> fK2KinC5combo[3];
//     Polarization vectors
    std::complex<double> fA[fMaxNRes][fMaxNDecays][2];
    std::complex<double> fB[fMaxNRes][fMaxNDecays][2];
//     std::complex<double> fA2[fMaxNRes][fMaxNDecays][2];
//     std::complex<double> fB2[fMaxNRes][fMaxNDecays][2];
//     Some intermediate shit for acceleration of symbolic calculation
    double cM[2][2][2][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro}; {S,D}; {0,1}(nSwap)
    std::complex<double> cF[2][2][2][2][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro}; {f,h}; {0,1}(nSwap); {iPhS,iPhD}
  };
}
#endif