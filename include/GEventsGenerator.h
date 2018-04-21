#ifndef GEVENTSGENERATOR_H
#define GEVENTSGENERATOR_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <cmath>
#include <TF1.h>
#include <TMath.h>
#include "GInterfaceForMinimization.h"

namespace Gamapola
{
  class GEventsGenerator: public GInterfaceForMinimization{
  public:
    GEventsGenerator();
    ~GEventsGenerator();
    void GSetEventsNumber(const int& nEvents);
    void GGenerateEvents();
    void GSetDecayMode(const std::vector<std::string>& resName, const std::vector<int>& charges);
//     void GSetCharge(const int& charge);
    
    void GSetPhotonPolarization(const double& lambda);
    void GSetQPMC(const double& gammaQPCM, const double& f2, 
                const double&  thetaK1, const double& phiDK, const double& phiSRho, 
                const double& phiDRho, const std::complex<double>& gg);
    void GSetK1_1270(const int& charge);
    void GSetK1_1400(const int& charge, const std::complex<double>& contribution);
    void GSetKStr_1410(const int& charge, const std::complex<double>& contribution,
                       const std::complex<double>& fractionOfKstr,const std::complex<double>& fractionOfRho);
    void GSetKStr_1680(const int& charge, const std::complex<double>& contribution,
                       const std::complex<double>& fractionOfKstr,const std::complex<double>& fractionOfRho);
    void GSetK2_1430(const int& charge, const std::complex<double>& contribution,
                       const std::complex<double>& fractionOfKstr,const std::complex<double>& fractionOfRho);
    
    const std::vector<double>& GGetCosThetaPDF() const
    {
      return fCosThetaPDF;
    }
    const std::vector<double>& GGetPhiPDF() const
    {
      return fPhiPDF;
    }
    const std::vector<double>& GGetPhiPDF2() const
    {
      return fPhiPDF2;
    }
    const std::vector<double>& GGetPDF1() const
    {
      return fPDF1;
    }
    const std::vector<double>& GGetOmega() const
    {
      return fOmega;
    }
    const std::vector<double>& GGetOmega2() const
    {
      return fOmega2;
    }
    
    const double& GGetIntegralOfPDF() const
    {
      return fIntegralOfPDF;
    }
    const int& GGetEventsNumber() const
    {
      return fNEvents;
    }
    const std::vector<double>& GGetCosThetaPDFFlat() const
    {
      return fCosThetaPDFFlat;
    }
    const std::vector<double>& GGetPhiPDFFlat() const
    {
      return fPhiPDFFlat;
    }
    const std::vector<double>& GGetPDF1Flat() const
    {
      return fPDF1Flat;
    }
    const std::vector<double>& GGetOmegaFlat() const
    {
      return fOmegaFlat;
    }
    const std::vector<double>& GGetOmega2Flat() const
    {
      return fOmega2Flat;
    }
    const std::vector<double>& GGetSijFlat() const
    {
      return fSijFlat;
    }
    const std::vector<double>& GGetSjkFlat() const
    {
      return fSjkFlat;
    }
    const std::vector<double>& GGetSikFlat() const
    {
      return fSikFlat;
    }
    const std::vector<double>& GGetSij() const
    {
      return fSij;
    }
    const std::vector<double>& GGetSjk() const
    {
      return fSjk;
    }
    const std::vector<double>& GGetSik() const
    {
      return fSik;
    }
    const std::vector<double>& GGetMFlat() const
    {
      return fMFlat;
    }
    const std::vector<double>& GGetM() const
    {
      return fM;
    }
    const std::complex<double>* GGetNormalizationIntegrals()
    {
      return fNormalizationIntegrals;
    }
  private:
    int fArg;
    int fNEvents;
//     limits on generation of statistical parameters
    double fLowCosThetaLimit;
    double fUpCosThetaLimit;
    double fLowPhiLimit;
    double fUpPhiLimit;
    double fLowSijLimit;
    double fUpSijLimit;
    double fLowSjkLimit;
    double fUpSjkLimit;
    double fLowSikLimit;
    double fUpSikLimit;
    double fLowSLimit;
    double fUpSLimit;
//     generated kinematic variables
    std::vector<double> fCosThetaPDF;
    std::vector<double> fPhiPDF;
    std::vector<double> fPhiPDF2;
    std::vector<double> fPDF1;
    std::vector<double> fOmega;
    std::vector<double> fOmega2;
//     generated flat kinematic variables
    std::vector<double> fCosThetaPDFFlat;
    std::vector<double> fPhiPDFFlat;
    std::vector<double> fPDF1Flat;
    std::vector<double> fOmegaFlat;
    std::vector<double> fOmega2Flat;
    std::vector<double> fSijFlat;
    std::vector<double> fSjkFlat;
    std::vector<double> fSikFlat;
    std::vector<double> fSij;
    std::vector<double> fSjk;
    std::vector<double> fSik;
    std::vector<double> fMFlat;
    std::vector<double> fM;
    //     Coefficients for importance sampling
    double fK1coeffOfij;
    double fAtancoeffOfij;
    double fK1coeffOfik;
    double fAtancoeffOfik;
    double fK1coeffOfjk;
    double fAtancoeffOfjk;
    double fPars[100];
    std::vector<std::string> fResNames;
    std::vector<int> fCharge;
  };
}
#endif