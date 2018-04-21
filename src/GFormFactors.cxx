#include "GFormFactors.h"

namespace Gamapola{
  GFormFactors::GFormFactors(): GInterfaceForMathFunctions()
  {
    std::cout << "GFormFactors constructor calling. . ." << std::endl;
  }
  GFormFactors::~GFormFactors()
  {
    std::cout << "GFormFactors destructor calling. . ." << std::endl;
  }
  
  void GFormFactors::GCalculateFormFactorsKstar(double* x,double* p, int nRes, int nSwap)
  {
    fA[nRes][0][nSwap] = GAKstar(p,nRes, nSwap);
    fB[nRes][0][nSwap] = GBKstar(x,p,nRes,nSwap);
  }
  void GFormFactors::GCalculateFormFactorsRho(double* x,double* p, int nRes)
  {
    fA[nRes][1][0] = GARho(p,nRes, 0);
    fB[nRes][1][0] = GBRho(x,p,nRes,0);
  }
  //   *****************************************f_V and h_V ********************************************
  std::complex<double> GFormFactors::GAKstar(double* p, int nRes, int nSwap)
  {
    std::complex<double> expPhaseDKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> aMatrixElement = fMK0star[nRes][0] + sqrt(0.5) * fMK0star[nRes][1] * expPhaseDKstar;
//     cF[0][nRes][0][0][nSwap] = -(cM[0][nRes][0][0][nSwap] + sqrt(0.5) * cM[0][nRes][0][1][nSwap] * expPhaseDKstar);
//     cF[1][nRes][0][0][nSwap] = -(cM[1][nRes][0][0][nSwap] + sqrt(0.5) * cM[1][nRes][0][1][nSwap] * expPhaseDKstar);
    cF[0][nRes][0][0][nSwap][0] = -cM[0][nRes][0][0][nSwap];
    cF[0][nRes][0][0][nSwap][1] = -sqrt(0.5) * cM[0][nRes][0][1][nSwap];
    cF[1][nRes][0][0][nSwap][0] = -cM[1][nRes][0][0][nSwap];
    cF[1][nRes][0][0][nSwap][1] = -sqrt(0.5) * cM[1][nRes][0][1][nSwap];
//     std::complex<double> aMatrixElement = (cF[0][nRes][0][0][nSwap][0] * sin(p[2]) + cF[1][nRes][0][0][nSwap][0] * cos(p[2]))
//     + (cF[0][nRes][0][0][nSwap][1] * sin(p[2]) + cF[1][nRes][0][0][nSwap][1] * cos(p[2])) * expPhaseDKstar;
//     std::cout << " In form-factors: " << (cF[0][nRes][0][0][nSwap][0] * sin(p[2]) + cF[1][nRes][0][0][nSwap][0] * cos(p[2]))
//     + (cF[0][nRes][0][0][nSwap][1] * sin(p[2]) + cF[1][nRes][0][0][nSwap][1] * cos(p[2])) * expPhaseDKstar << " " << aMatrixElement << std::endl;
    return -aMatrixElement;
  }
  std::complex<double> GFormFactors::GBKstar(double* x,double* p, int nRes, int nSwap)
  {
    double mij = sqrt(x[2]);
    double M = sqrt(x[4]);
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sBrace1 = 1 - mij / fEnergyA[0];
    double sBrace2 =  (1 + 2 * mij / fEnergyA[0]);
    double sTerm = fEnergyA[0] / ( M * fMomentaA[0] * fMomentaA[0] );
    std::complex<double> aMatrixElement = (fMK0star[nRes][0] * sBrace1 + 
    sqrt(0.5) * fMK0star[nRes][1] * sBrace2 * expPhaseSKstar) * sTerm;
    
    cF[0][nRes][0][1][nSwap][0] = -cM[0][nRes][0][0][nSwap] * sBrace1 * sTerm;
    cF[0][nRes][0][1][nSwap][1] = -sqrt(0.5) * cM[0][nRes][0][1][nSwap] * sBrace2 * sTerm;
    cF[1][nRes][0][1][nSwap][0] = -cM[1][nRes][0][0][nSwap] * sBrace1 * sTerm;
    cF[1][nRes][0][1][nSwap][1] = -sqrt(0.5) * cM[1][nRes][0][1][nSwap] * sBrace2 * sTerm;
    
//     -(cM[0][nRes][0][0][nSwap] * sBrace1[nSwap][0] + 
//     sqrt(0.5) * cM[0][nRes][0][1][nSwap] * sBrace2[nSwap][0] * expPhaseSKstar) * sTerm[nSwap][0];
    
//     cF[1][nRes][0][1][nSwap][] = -(cM[1][nRes][0][0][nSwap] * sBrace1[nSwap][0] + 
//     sqrt(0.5) * cM[1][nRes][0][1][nSwap] * sBrace2[nSwap][0] * expPhaseSKstar) * sTerm[nSwap][0];
//     std::complex<double> aMatrixElement = (cF[0][nRes][0][1][nSwap][0] * sin(p[2]) + cF[1][nRes][0][1][nSwap][0] * cos(p[2]))
//     + (cF[0][nRes][0][1][nSwap][1] * sin(p[2]) + cF[1][nRes][0][1][nSwap][1] * cos(p[2])) * expPhaseSKstar;
    return -aMatrixElement;
  }
  std::complex<double> GFormFactors::GARho(double* p, int nRes, int nSwap)
  {
//     std::complex<double> ff2(p[6], p[10]);
//     std::cout << ff2 << std::endl;
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> aMatrixElement = /*ff2**/(fMRho0[nRes][0] + sqrt(0.5) * fMRho0[nRes][1] * expPhaseDRho) * expPhaseSRho;
//     cF[0][nRes][1][0][nSwap] = -(cM[0][nRes][1][0][nSwap] + sqrt(0.5) * cM[0][nRes][1][1][nSwap] * expPhaseDRho) * expPhaseSRho;
//     cF[1][nRes][1][0][nSwap] = -(cM[1][nRes][1][0][nSwap] + sqrt(0.5) * cM[1][nRes][1][1][nSwap] * expPhaseDRho) * expPhaseSRho;
    
    cF[0][nRes][1][0][nSwap][0] = -cM[0][nRes][1][0][nSwap];
    cF[0][nRes][1][0][nSwap][1] = -sqrt(0.5) * cM[0][nRes][1][1][nSwap];
    cF[1][nRes][1][0][nSwap][0] = -cM[1][nRes][1][0][nSwap];
    cF[1][nRes][1][0][nSwap][1] = -sqrt(0.5) * cM[1][nRes][1][1][nSwap];
//     std::complex<double> aMatrixElement = (cF[0][nRes][1][0][nSwap][0] * sin(p[2]) + cF[1][nRes][1][0][nSwap][0] * cos(p[2])) * expPhaseSRho
//     + (cF[0][nRes][1][0][nSwap][1] * sin(p[2]) + cF[1][nRes][1][0][nSwap][1] * cos(p[2])) * expPhaseSRho * expPhaseDRho;
//     std::cout << -aMatrixElement << std::endl;
    return -aMatrixElement;
  }
  std::complex<double> GFormFactors::GBRho(double* x,double* p, int nRes, int nSwap)
  {
//     std::complex<double> ff2(p[6], p[10]);
//     std::cout << ff2 << std::endl;
    double mij = sqrt(fSij);
    double M = sqrt(x[4]);
    std::complex<double> expPhaseDKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sBrace1 = (1 - mij / fEnergyA[1]);
    double sBrace2 =  (1 + 2 * mij / fEnergyA[1]);
    double sTerm = fEnergyA[1] / ( M * fMomentaA[1] * fMomentaA[1] );
    std::complex<double> aMatrixElement = /*ff2**/(fMRho0[nRes][0] * sBrace1 + 
    sqrt(0.5) * fMRho0[nRes][1] * sBrace2 * expPhaseDRho) * expPhaseSRho * sTerm;
    
//     cF[0][nRes][1][1][nSwap] = -( cM[0][nRes][1][0][nSwap] * sBrace1[nSwap][1] + 
//     sqrt(0.5) * cM[0][nRes][1][1][nSwap] * sBrace2[nSwap][1] * expPhaseDRho) * expPhaseSRho * sTerm[nSwap][1];
//     
//     cF[1][nRes][1][1][nSwap] = -( cM[1][nRes][1][0][nSwap] * sBrace1[nSwap][1] + 
//     sqrt(0.5) * cM[1][nRes][1][1][nSwap] * sBrace2[nSwap][1] * expPhaseDRho) * expPhaseSRho * sTerm[nSwap][1];
    cF[0][nRes][1][1][nSwap][0] = -cM[0][nRes][1][0][nSwap] * sBrace1* sTerm;
    cF[0][nRes][1][1][nSwap][1] = -sqrt(0.5) * cM[0][nRes][1][1][nSwap] * sBrace2 * sTerm;
    cF[1][nRes][1][1][nSwap][0] = -cM[1][nRes][1][0][nSwap] * sBrace1 * sTerm;
    cF[1][nRes][1][1][nSwap][1] = -sqrt(0.5) * cM[1][nRes][1][1][nSwap] * sBrace2 * sTerm;
    
//     std::complex<double> aMatrixElement = (cF[0][nRes][1][1][nSwap][0] * sin(p[2]) + cF[1][nRes][1][1][nSwap][0] * cos(p[2])) * expPhaseSRho
//     + (cF[0][nRes][1][1][nSwap][1] * sin(p[2]) + cF[1][nRes][1][1][nSwap][1] * cos(p[2])) * expPhaseSRho * expPhaseDRho;
    
    return -aMatrixElement;
  }
}