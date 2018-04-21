#ifndef GFORMFACTORS_H
#define GFORMFACTORS_H

#include "GInterfaceForMathFunctions.h"

namespace Gamapola
{
  class GFormFactors: virtual public GInterfaceForMathFunctions{
  private:
    int fArg;
    double fsqrt05;
  protected:
    GFormFactors();
    virtual ~GFormFactors();
    std::complex<double> GAKstar(double* p, int nRes, int nSwap);
    std::complex<double> GBKstar(double* x,double* p, int nRes, int nSwap);
    std::complex<double> GARho(double* p, int nRes, int nSwap);
    std::complex<double> GBRho(double* x,double* p, int nRes, int nSwap);
    
    void GCalculateFormFactorsKstar(double* x,double* p, int nRes, int nSwap);
    void GCalculateFormFactorsRho(double* x,double* p, int nRes);
    
    double fEnergyA[fMaxNDecays];
  };
}
#endif