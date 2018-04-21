#include "GInterfaceForMinimization.h"
#include <iomanip>

namespace Gamapola{
  GInterfaceForMinimization::GInterfaceForMinimization():
  fArg(0),
  fMf(0),
  fMemFunctionName(0),
  fMemFunctionClassName(0),
  fNameOfMinimizer(0),
    
  fNParameters(0),
  fParsValues(0),
    
  fNVariables(0),
  fInitialVarValues(0),
  fLowVarsLimits(0),
  fUpVarsLimits(0),
  fVarsSteps(0),
  fKindsOfVars(0),
  fNamesOfVars(0),
    
  fMinFunctionValue(0),
  fMaxFunctionValue(0),
  fFitParameters(0),
  fCovarianceMatrix(0),
  fFuncPointer(0)
  {
    fNRuns = 1e6;
    std::cout << "GInterfaceForMinimization constructor calling. . ." << std::endl;
  }
  GInterfaceForMinimization::~GInterfaceForMinimization()
  {
    delete [] fFitParameters;
    delete fFuncPointer;
    delete [] fCovarianceMatrix;
    delete [] fParsValues;
    std::cout << "GInterfaceForMinimization destructor calling. . ." << std::endl;
  }
  void GInterfaceForMinimization::GSetFunction(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p), 
                                               GMathFunctions** mf,const char* memFunctionName, const char* memFunctionClassName)
  {
    fMf = mf;
    fMemFunction = memFunction;
    fMemFunctionName = memFunctionName;
    fMemFunctionClassName = memFunctionClassName;
    std::cout << fMemFunctionName << "  " << fMemFunctionClassName << std::endl;
  }
  void GInterfaceForMinimization::GSetFunctionParameters(const int& nPars, double* pars)
  {
    fNParameters = nPars;
    fParsValues = new double[fNParameters];
    for(auto i = 0; i < fNParameters; i++)
      fParsValues[i] = pars[i];
  }
  
  void GInterfaceForMinimization::GSetFunctionVariables(const int& nVars, double* initialVarValues, 
                               double* lowVarsLimits, double* upVarsLimits, 
                               double* varsSteps, const char* kindsOfVars[], const char* namesOfVars[])
  {
    fInitialVarValues = initialVarValues;
    fLowVarsLimits = lowVarsLimits;
    fUpVarsLimits = upVarsLimits;
    fVarsSteps = varsSteps;
    fKindsOfVars = kindsOfVars;
    fNamesOfVars = namesOfVars;
    fNVariables = nVars;
    fFitParameters = new double[fNVariables];
    fCovarianceMatrix = new double[fNVariables*fNVariables];
  }
  
  void GInterfaceForMinimization::GSetMinimizer(const char* nameOfMinimizer)
  {
    fNameOfMinimizer = nameOfMinimizer;
  }
  
  void GInterfaceForMinimization::GSetNRuns(const int& nRuns)
  {
    fNRuns = nRuns;
  }
  void GInterfaceForMinimization::GPrintErrorMatrix()
  {
    std::ofstream ofs;
    ofs.open("CovarianceMatrix.txt", std::ios::out);
    ofs << "____________________________________________________Error Matrix_________________________________________________" << std::endl;
    std::vector<double> errorVector;
    for(int i = 0; i < fNVariables; ++i)
    {
      for(int j = 0; j < fNVariables; ++j)
        if(*(fCovarianceMatrix + i*fNVariables + j) != 0)
          errorVector.push_back(*(fCovarianceMatrix + i*fNVariables + j));
//         std::cout << *(fCovarianceMatrix + i*fNVariables + j) << "   ";
//       std::cout << std::endl;
    }
    ofs << std::setw(20);
    for(size_t i = 0; i < fErrorMatrix.size(); i++)
      ofs << fparNames[i] << std::setw(15);
    ofs << std::endl;
    for(size_t i = 0; i < fErrorMatrix.size(); i++)
    {
      ofs << std::setw(10) << fparNames[i] << "   ";
      for(size_t j = 0; j < fErrorMatrix.size(); ++j)
      {
        fErrorMatrix[i][j] = errorVector[ i * fErrorMatrix.size() + j];
        if(i != j)
          fErrorMatrix[i][j] = fErrorMatrix[i][j] / 
          sqrt(errorVector[ i * fErrorMatrix.size() + i] * errorVector[ j * fErrorMatrix.size() + j]);
        else
          fErrorMatrix[i][j] = sqrt(errorVector[ i * fErrorMatrix.size() + j]);
        ofs << fErrorMatrix[i][j] << "   ";
      }
      ofs << std::endl;
    }
     ofs << "_________________________________________________________________________________________________________________" << std::endl;
     ofs.close();
  }
  void GInterfaceForMinimization::GMinimizeWithMinuit()
  {
    fFuncPointer = new TF1(fMemFunctionName, *fMf, fMemFunction, 0, 0, fNParameters, fMemFunctionClassName, fMemFunctionName);
    fFuncPointer->SetParameters(fParsValues);
    ROOT::Minuit2::Minuit2Minimizer globalMin(ROOT::Minuit2::kMigrad);
    ROOT::Math::WrappedMultiTF1 g1(*fFuncPointer, fNVariables);
    globalMin.SetFunction(g1);
    globalMin.SetMaxIterations(1000);
    globalMin.SetTolerance(0.5);
    unsigned int numberOfFreedVariables = 0;
    for(int i = 0; i < fNVariables; ++i)
    {
      std::cout << fNamesOfVars[i] << "   " << fKindsOfVars[i] << "  " << fInitialVarValues[i] << std::endl;
      if(std::strcmp(fKindsOfVars[i], "limited") == 0)
      {
        globalMin.SetLimitedVariable(i, fNamesOfVars[i], fInitialVarValues[i],fVarsSteps[i], fLowVarsLimits[i], fUpVarsLimits[i]);
        ++numberOfFreedVariables;
        fparNames.push_back(fNamesOfVars[i]);
      }
      else if(std::strcmp(fKindsOfVars[i], "fixed") == 0)
        globalMin.SetFixedVariable(i, fNamesOfVars[i], fInitialVarValues[i]);
      else
      {
        globalMin.SetVariable(i, fNamesOfVars[i], fInitialVarValues[i],fVarsSteps[i]);
        fparNames.push_back(fNamesOfVars[i]);
        numberOfFreedVariables++;
      }
    }
    fErrorMatrix.resize(numberOfFreedVariables);
    for(size_t i = 0; i < numberOfFreedVariables; ++i)
      fErrorMatrix[i].resize(numberOfFreedVariables, 0.);
    globalMin.Minimize();
    globalMin.GetCovMatrix(fCovarianceMatrix);
    std::cout << std::scientific;
    GPrintErrorMatrix();
    globalMin.PrintResults();
    for(int i = 0; i < fNVariables; ++i)
      fFitParameters[i] = globalMin.X()[i];
    fMinFunctionValue = globalMin.MinValue();
    fStatus = globalMin.Status();
    fEdm = globalMin.Edm();
  }
//   only applied for p.d.f > 0, for -2Log(L) won't work...
  void GInterfaceForMinimization::GMinimizeWithMonteCarlo() // non universal method!!!
  {
    fFuncPointer = new TF1(fMemFunctionName, *fMf, fMemFunction, 0, 0, fNParameters, fMemFunctionClassName, fMemFunctionName);
    double maxFunctionValue = -1e300;
    double minFunctionValue = 1e300;
    double* varsMin = new double[fNVariables];
    double* varsMax = new double[fNVariables];
    ROOT::Math::ParamFunctor functor = ROOT::Math::ParamFunctor(*fMf,fMemFunction);
    double sij,sjk, sik; // Kpi1, pi1pi2, K1pi2
    double lowsjk, upsjk;
    double mi = (*fMf)->GGetMasses()[0];
    double mj = (*fMf)->GGetMasses()[1];
    double mk = (*fMf)->GGetMasses()[2];
    double M = (*fMf)->GGetMasses()[3];
    double pdfVal = 0.; 
    for(int i = 0; i < fNRuns; ++i)
    {
      fInitialVarValues[4] = fLowVarsLimits[4] + (fUpVarsLimits[4]-fLowVarsLimits[4])*double(rand())/RAND_MAX;
      M = sqrt(fInitialVarValues[4]);
      fUpVarsLimits[2] = ( M - mk ) * ( M - mk );
      fUpVarsLimits[3] = ( M - mi ) * ( M - mi );
      
      sjk = fLowVarsLimits[3] + (fUpVarsLimits[3] - fLowVarsLimits[3])*double(rand())/RAND_MAX;
      sij = fLowVarsLimits[2] + (fUpVarsLimits[2] - fLowVarsLimits[2])*double(rand())/RAND_MAX;
      lowsjk = (*fMf)->GSjkMin(sij, mi, mj, mk, M);
      upsjk = (*fMf)->GSjkMax(sij, mi, mj, mk, M);
      if((sjk < lowsjk) || (sjk > upsjk))
        continue;
      sik = mi*mi+mj*mj+mk*mk+M*M-sjk-sij;
      fInitialVarValues[0] = fLowVarsLimits[0] + (fUpVarsLimits[0]-fLowVarsLimits[0])*double(rand())/RAND_MAX;
      fInitialVarValues[1] = fLowVarsLimits[1] + (fUpVarsLimits[1]-fLowVarsLimits[1])*double(rand())/RAND_MAX;
      fInitialVarValues[2] = sij;
      fInitialVarValues[3] = sjk;
      fInitialVarValues[4] = M*M;
      pdfVal = functor(fInitialVarValues, fParsValues);
      if(maxFunctionValue < pdfVal)
      {
        maxFunctionValue = pdfVal;
        std::cout << "MAX: " <<  maxFunctionValue << std::endl;
        for(int iVars = 0; iVars < fNVariables; iVars++)
          varsMax[iVars] = fInitialVarValues[iVars];
      }
      if(minFunctionValue > pdfVal)
      {
        minFunctionValue = pdfVal;
        std::cout << "MIN: " << minFunctionValue << "  " << i << std::endl;
        for(int iVars = 0; iVars < fNVariables; iVars++)
        {
          varsMin[iVars] = fInitialVarValues[iVars];
          fFitParameters[iVars] = fInitialVarValues[iVars];
        }
      }
    }
    fMinFunctionValue = minFunctionValue;
    fMaxFunctionValue = maxFunctionValue;
    delete varsMax;
    delete varsMin;
  }
  void GInterfaceForMinimization::GWriteNormalizationIntegrals(const char* filename)
  {
    fReadWrite = false;
    fFileWithIntegrals = filename;
  }
  void GInterfaceForMinimization::GReadNormalizationIntegrals(const char* filename)
  {
    fReadWrite = true;
    fFileWithIntegrals = filename;
  }
}