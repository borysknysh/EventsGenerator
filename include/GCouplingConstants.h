#ifndef GCOUPLINGCONSTANTS_H
#define GCOUPLINGCONSTANTS_H

#include "GInterfaceForMathFunctions.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"

namespace Gamapola
{
  class GCouplingConstants: virtual public GInterfaceForMathFunctions{
  private:
    int fArg;
    double fs;
  protected:
    GCouplingConstants();
    virtual ~GCouplingConstants();
    void GCouplingConstantsCalculation();
    double GCouplingKappa(double sLow, double sUp);
    double GInnerIntegral(double sKappa);
    double GOuterIntegral(double s);
    void GHadronicFormFactor(double s, int nRes);
    
    double fResBrs[fMaxNRes][fMaxNDecays];
    double fHadronicFF[fMaxNRes][fMaxNDecays];
  };
}
#endif
