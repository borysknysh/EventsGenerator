#ifndef GBREITWIGNER_H
#define GBREITWIGNER_H

#include "GInterfaceForMathFunctions.h"
#include "GKinematics.h"
#include <memory>

namespace Gamapola
{
  class GBreitWigner: virtual public GInterfaceForMathFunctions{
  public:
    GBreitWigner();
    virtual ~GBreitWigner();
    std::complex<double> GBWKstr(double sij);
    std::complex<double> GBWRho(double sij);
    std::complex<double> GBWK1(double s, int nRes);
    std::complex<double> GBWKappa(double sij);
  private:
    int fArg;
  protected:
    double fResMasses[fMaxNRes];
  };
}
#endif
