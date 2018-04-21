#ifndef GINTERFACEFORMINIMIZATION_H
#define GINTERFACEFORMINIMIZATION_H

#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "GMathFunctions.h"
#include "Math/WrappedMultiTF1.h"
#include "Minuit2/Minuit2Minimizer.h"

namespace Gamapola
{
  class GInterfaceForMinimization{
  public:
    GInterfaceForMinimization();
    virtual ~GInterfaceForMinimization();
    virtual void GSetFunction(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p), GMathFunctions** mf,
                      const char* memFunctionName, const char* memFunctionClassName);
    virtual void GSetFunctionParameters(const int& nPars, double* pars);
    virtual void GSetFunctionVariables(const int& nVars, double* initialVarValues, 
                               double* lowVarsLimits, double* upVarsLimits, 
                               double* varsSteps, const char* kindsOfVars[], const char* namesOfVars[]);
    virtual void GSetMinimizer(const char* nameOfMinimizer = "mc");
    virtual void GSetNRuns(const int& nRuns);
    void GWriteNormalizationIntegrals(const char* filename);
    void GReadNormalizationIntegrals(const char* filename);
    const std::vector<std::vector<double> >& GGetErrorMatrix() const
    {
      return fErrorMatrix;
    }
    
    const double& GGetMinFuncValue() const
    {
      return fMinFunctionValue;
    }
    const double& GGetMaxFuncValue() const
    {
      return fMaxFunctionValue;
    }
    double* GGetFitParameters()
    {
      return fFitParameters;
    }
    TF1* GGetFuncPointer() const
    {
      return fFuncPointer;
    }
    const double* GGetCovarianceMatrix()
    {
      return fCovarianceMatrix;
    }
    const int& GGetStatus() const
    {
      return fStatus;
    }
    const double& GGetEdm()
    {
      return fEdm;
    }
    void GPrintErrorMatrix();
  protected:
    void GMinimizeWithMinuit();
    void GMinimizeWithMonteCarlo();
    
    int fArg;
    double (Gamapola::GMathFunctions::*fMemFunction)(double* x, double* p);
    GMathFunctions** fMf;
    const char* fMemFunctionName;
    const char* fMemFunctionClassName;
    const char* fNameOfMinimizer;
    
    int fNParameters;
    double* fParsValues;
    
    int fNVariables;
    double* fInitialVarValues;
    double* fLowVarsLimits;
    double* fUpVarsLimits; 
    double* fVarsSteps;
    const char** fKindsOfVars;
    const char** fNamesOfVars;
    
    double fMinFunctionValue;
    double fMaxFunctionValue;
    double* fFitParameters;
    double* fCovarianceMatrix;
    std::vector<std::vector<double> > fErrorMatrix;
    std::vector<const char*> fparNames;
    TF1* fFuncPointer;
    int fStatus;
    //     some generated values
    double fIntegralOfPDF;
    std::complex<double>* fNormalizationIntegrals;
    int fNIntegrals;
    bool fReadWrite;
    const char* fFileWithIntegrals;
    int fNRuns;
    double fEdm;
  };
}
#endif