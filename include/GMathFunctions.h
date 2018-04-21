#ifndef GMATHFUNCTIONS_H
#define GMATHFUNCTIONS_H

#include "GInterfaceForMathFunctions.h"
#include "GKinematics.h"
#include "GBreitWigner.h"
#include "GCouplingConstants.h"
#include "GQuarkPairCreationModel.h"
#include "GFormFactors.h"
#include "GCurrents.h"

namespace Gamapola
{
  class GMathFunctions: public GKinematics, 
                        public GBreitWigner, 
                        public GCouplingConstants,
                        public GQuarkPairCreationModel,
                        public GFormFactors,
                        public GCurrents{
  public:
    GMathFunctions();
    virtual ~GMathFunctions();
    double GOmegaPDF1();
    double GCovarianceMatrix(double* x, double* p);
    double GLogLikelihoodPDF2(double* x, double* p);
    double GProcessingComputationOfPDF(double* x, double* p);
    double GMomentaConst(double* u);    
    void GCalculateResonance(double* x, double* p, const int& iRes);
//     https://stackoverflow.com/questions/7304511/partial-derivatives --- will be useful=)
//     https://dou.ua/lenta/articles/it-position-ds-ml/ --- machine learning
    void GSetMasses(const double& mi, const double& mj, const double& mk, const double& M);
    void GSetNumberOfNormalizationIntegrals(const int& nTerms);
    void GSetNormalizationIntegrals(const std::complex<double>* normIntegrals);
    void GSetResonances(const std::vector<std::string>& resName, const std::vector<int>& charges);
    const double* GGetMasses() const
    {
      return fMasses;
    }
    const std::complex<double>* GGetNormalizationIntegrals() const
    {
      return fNormalizationIntegrals;
    }
    const std::vector<double (Gamapola::GMathFunctions::*) (double* x, double* p)>& GGetFunctionPointers() const
    {
      return fFunctionPointers;
    }
    std::vector<std::complex<double> > GGetKinematicalCoefficients();
    
  private:
    std::vector<double (Gamapola::GMathFunctions::*) (double* x, double* p)> fFunctionPointers;    
    std::vector<std::vector<std::vector<double> > > fConstParsForI01; 
    std::vector<std::string> fResNames; 
    int fNRes;
    int fIndexControl[fMaxNRes];
  };
}
#endif