#include "GQuarkPairCreationModel.h"

namespace Gamapola{
  GQuarkPairCreationModel::GQuarkPairCreationModel(): GInterfaceForMathFunctions()
  {
    std::cout << "GQuarkPairCreationModel constructor calling. . ." << std::endl;
    fFunctionPointersKstrNR[0][0] = &GQuarkPairCreationModel::GMSKstarNR1270;
    fFunctionPointersKstrNR[0][1] = &GQuarkPairCreationModel::GMDKstarNR1270;
    fFunctionPointersRhoNR[0][0] = &GQuarkPairCreationModel::GMSRhoNR1270;
    fFunctionPointersRhoNR[0][1] = &GQuarkPairCreationModel::GMDRhoNR1270;
    fFunctionPointersKstr[0][0] = &GQuarkPairCreationModel::GMSKstar1270;
    fFunctionPointersKstr[0][1] = &GQuarkPairCreationModel::GMDKstar1270;
    fFunctionPointersRho[0][0] = &GQuarkPairCreationModel::GMSRho1270;
    fFunctionPointersRho[0][1] = &GQuarkPairCreationModel::GMDRho1270;
    fFunctionPointersKstrNR[1][0] = &GQuarkPairCreationModel::GMSKstarNR1400;
    fFunctionPointersKstrNR[1][1] = &GQuarkPairCreationModel::GMDKstarNR1400;  
    fFunctionPointersRhoNR[1][0] = &GQuarkPairCreationModel::GMSRhoNR1400;
    fFunctionPointersRhoNR[1][1] = &GQuarkPairCreationModel::GMDRhoNR1400;      
    fFunctionPointersKstr[1][0] = &GQuarkPairCreationModel::GMSKstar1400;
    fFunctionPointersKstr[1][1] = &GQuarkPairCreationModel::GMDKstar1400;    
    fFunctionPointersRho[1][0] = &GQuarkPairCreationModel::GMSRho1400;
    fFunctionPointersRho[1][1] = &GQuarkPairCreationModel::GMDRho1400;
  }
  GQuarkPairCreationModel::~GQuarkPairCreationModel()
  {
    std::cout << "GQuarkPairCreationModel destructor calling. . ." << std::endl;
  }
  
  
  void GQuarkPairCreationModel::GCalculateQPMC(double* p, const int& iRes, const int& swap)
  {
    for(size_t i = 0; i < 2/*fMomentaConst[fIndex[ires]].size()-1*/; i++)
      { 
        fI0[0][i] = GI0(iRes, i);
        fI1[0][i] = GI1(iRes, i);
        fMS[0][i] = GMS(p,iRes,i,swap);
        fMD[0][i] = GMD(p,iRes,i,swap);
//         std::cout <<"S,D waves: " << fMS[0][i] << "   " << fMD[0][i] << std::endl;
      }
    for(int iWave = 0; iWave < 2; iWave++)
    {
//       Non-relativistic matrix elements
      fMK0starNR[iRes][iWave] = (this->*fFunctionPointersKstrNR[iRes][iWave])(swap);
      fMRho0NR[iRes][iWave] = (this->*fFunctionPointersRhoNR[iRes][iWave])(swap);   
//       Relativistic matrix elements
      fMK0star[iRes][iWave] = (this->*fFunctionPointersKstr[iRes][iWave])(swap);
      fMRho0[iRes][iWave] = (this->*fFunctionPointersRho[iRes][iWave])(swap);  // why zero?
//       std::cout << "MK0star: " << fMK0star[iRes][iWave] << " for wave " << iWave << " res " << iRes << std::endl;
//       std::cout << "MRho0: " << fMRho0[iRes][iWave] << " for wave " << iWave << " res " << iRes << std::endl;
    }
  }
  
  void GQuarkPairCreationModel::GQPCMConstanstsCalculation(std::vector<double>* u, int nRes)
  {
    double ra = u[0][0];
    double rb = u[0][1];
    double rc = u[0][2];
    double r2 = (ra*ra+rb*rb+rc*rc);
    double coeff = -4. / sqrt( sqrt( TMath::Pi() * TMath::Pi() * TMath::Pi() * TMath::Pi() * TMath::Pi() ) / 3. );
    fcoeff[nRes][0] = coeff * sqrt( ra * ra * ra * ra * ra ) * sqrt( rb * rc * rb * rc * rb * rc ) / sqrt( r2 * r2 * r2 * r2 * r2 );
    fbracePart[nRes][0] = ( r2 + ra * ra ) * ( rb * rb + rc * rc ) / ( 4 * r2 );
    fexpPart[nRes][0] = ra * ra * (rb * rb + rc * rc) / (8. * r2);
    
    fcoeff[nRes][1] = fcoeff[nRes][0];
    fbracePart[nRes][1] = fbracePart[nRes][0];
    fexpPart[nRes][1] = fexpPart[nRes][0];
    
    fcoeffAdd[nRes] = -sqrt(2.) / ( 3. * TMath::Pi() * sqrt( sqrt( TMath::Pi() ) ) ) * ra * ra * rb * rc * rb * sqrt( ra ) * sqrt( rb * rc  ) / 
    ( r2 * r2 * sqrt( r2 ) );
    
    fcoeff[nRes][2] = sqrt(2.) * fcoeffAdd[nRes] * (3 * ra * ra - 2 * ( rb * rb + rc * rc ) ) / r2;
    
    fbracePart[nRes][2] = ra * ra * ( rb * rb + rc * rc ) * ( 2 * ra * ra + rb * rb + rc * rc ) / 
    ( 4 * r2 * (3 * ra * ra - 2 * ( rb * rb + rc * rc ) ) );
    fexpPart[nRes][2] = fexpPart[nRes][0];
    fsqrt2 = sqrt(2.);
//     std::cout << r2 << "  " << fcoeffAdd[nRes] << "   " << fbracePart[nRes][2] << "   " << fexpPart[nRes][2] << std::endl;
  }
  double GQuarkPairCreationModel::GI0(int nRes, int nFnc)
  {    
    double mom2 = fMomentaA[nFnc] * fMomentaA[nFnc];
    double brace = 1 - mom2 * fbracePart[nRes][nFnc];
    double expo = exp( -mom2 * fexpPart[nRes][nFnc] );
//     std::cout << fMomentaA[nFnc] << std::endl;
    return fcoeff[nRes][nFnc] * brace * expo;
  }
  double GQuarkPairCreationModel::GI1(int nRes, int nFnc)
  {
    double mom2 = fMomentaA[nFnc] * fMomentaA[nFnc];
    double expo = exp( -mom2 * fexpPart[nRes][nFnc] );
    return -fcoeff[nRes][nFnc] * expo;
  }
  double GQuarkPairCreationModel::GMS(double* p, int nRes, int nFunc, int nSwap)
  {
    double gammaQPC = p[0];
    double f2 = p[1];
    double mom2 = fMomentaA[nFunc] * fMomentaA[nFunc];
    double sMS = sqrt( 1.5 ) * ( 2 * fI1[0][nFunc] - fI0[0][nFunc] ) / 18.;
    sExp[nSwap][nRes][nFunc] = fMomentaConst[nRes][nFunc] * fMomentaConst[nRes][nFunc] - mom2;
    double MS = gammaQPC * sMS * exp( sExp[nSwap][nRes][nFunc] * f2 );
//     std::cout << mom2 << "  " << std::endl;
    return MS;
  }
  double GQuarkPairCreationModel::GMD(double* p, int nRes, int nFunc, int nSwap)
  {
    double gammaQPC = p[0];
    double f2 = p[1];
    double sMD = sqrt( 1.5 ) * ( fI1[0][nFunc] + fI0[0][nFunc] ) / 18.;
    double MD = gammaQPC * sMD * exp( sExp[nSwap][nRes][nFunc] * f2 );
    return MD;
  }
//   Non-relativistic S and D matrix elements ************************************************************
  double GQuarkPairCreationModel::GMSKstarNR1270(int nSwap)
  {
    cM[0][0][0][0][nSwap] = fMS[0][0] * fsqrt2;
    cM[1][0][0][0][nSwap] = -fMS[0][0];
//     std::cout << fMS[0][0] * ( fsqrt2 * sin( thetaK1 ) - cos( thetaK1 ) ) << std::endl;
    return cM[0][0][0][0][nSwap] * fSinThK1 + cM[1][0][0][0][nSwap] * fCosThK1;
//     return fMS[0][0] * ( fsqrt2 * sin( p[2] ) - cos( p[2] ) )/* * fExponentialM[0][0]*/;
  }
  double GQuarkPairCreationModel::GMSKstarNR1400(int nSwap)
  {
    cM[0][1][0][0][nSwap] = fMS[0][0];
    cM[1][1][0][0][nSwap] = fMS[0][0] * fsqrt2;
//     std::cout << fMS[0][0]  << std::endl;
    return cM[0][1][0][0][nSwap] * fSinThK1 + cM[1][1][0][0][nSwap] * fCosThK1;
//     return fMS[0][0] * ( fsqrt2 * cos( p[2] ) + sin( p[2] ) )/* * fExponentialM[1][0]*/;
  }
  double GQuarkPairCreationModel::GMDKstarNR1270(int nSwap)
  {
    cM[0][0][0][1][nSwap] = -fMD[0][0];
    cM[1][0][0][1][nSwap] = -fsqrt2 * fMD[0][0];
    return cM[0][0][0][1][nSwap] * fSinThK1 + cM[1][0][0][1][nSwap] * fCosThK1;
//     return fMD[0][0] * ( -sin( p[2] ) - fsqrt2 * cos( p[2] ) )/* * fExponentialM[0][0]*/;
  }
  double GQuarkPairCreationModel::GMDKstarNR1400(int nSwap)
  {
    cM[0][1][0][1][nSwap] = fMD[0][0] * fsqrt2;
    cM[1][1][0][1][nSwap] = -fMD[0][0];
    return cM[0][1][0][1][nSwap] * fSinThK1 + cM[1][1][0][1][nSwap] * fCosThK1;
//     return fMD[0][0] * ( -cos( p[2] ) + fsqrt2 * sin( p[2] ) )/* * fExponentialM[1][0]*/;
  }
  double GQuarkPairCreationModel::GMSRhoNR1270(int nSwap)
  {
    cM[0][0][1][0][nSwap] = fMS[0][1] * fsqrt2;
    cM[1][0][1][0][nSwap] = fMS[0][1];
//     std::cout << fMS[0][1] << std::endl;
    return cM[0][0][1][0][nSwap] * fSinThK1 + cM[1][0][1][0][nSwap] * fCosThK1;
//     return fMS[0][1] * ( fsqrt2 * sin( p[2] ) + cos( p[2] ) )/* * fExponentialM[0][1]*/;
  }
  double GQuarkPairCreationModel::GMSRhoNR1400(int nSwap)
  {
    cM[0][1][1][0][nSwap] = -fMS[0][1];
    cM[1][1][1][0][nSwap] = fMS[0][1] * fsqrt2;
//     std::cout << fExponentialM[1][1] << std::endl;
    return cM[0][1][1][0][nSwap] * fSinThK1 + cM[1][1][1][0][nSwap] * fCosThK1;
//     return fMS[0][1] * ( fsqrt2 * cos( p[2] ) - sin( p[2] ) )/* * fExponentialM[1][1]*/;
  }
  double GQuarkPairCreationModel::GMDRhoNR1270(int nSwap)
  {
    cM[0][0][1][1][nSwap] = -fMD[0][1];
    cM[1][0][1][1][nSwap] = fMD[0][1] * fsqrt2;
    return cM[0][0][1][1][nSwap] * fSinThK1 + cM[1][0][1][1][nSwap] * fCosThK1;
//     return fMD[0][1] * ( -sin( p[2] ) + fsqrt2 * cos( p[2] ) ) /** fExponentialM[0][1]*/;
  }
  double GQuarkPairCreationModel::GMDRhoNR1400(int nSwap)
  {
    cM[0][1][1][1][nSwap] = -fMD[0][1] * fsqrt2;
    cM[1][1][1][1][nSwap] = -fMD[0][1];
    return cM[0][1][1][1][nSwap] * fSinThK1 + cM[1][1][1][1][nSwap] * fCosThK1;
//     return fMD[0][1] * ( -cos( p[2] ) - fsqrt2 * sin( p[2] ) )/* * fExponentialM[1][1]*/;
  }
//   ***************************Relativistic S and D matrix elements******************************************************************
double GQuarkPairCreationModel::GMSKstar1270(int nSwap)
  {
    cM[0][0][0][0][nSwap] *= fEnergyCoeff[0][0];
    cM[1][0][0][0][nSwap] *= fEnergyCoeff[0][0];
//     std::cout << fEnergyCoeff[0][0] * fMK0starNR[0][0] << std::endl;
    return fEnergyCoeff[0][0] * fMK0starNR[0][0];
  }
  double GQuarkPairCreationModel::GMSKstar1400(int nSwap)
  {
    cM[0][1][0][0][nSwap] *= fEnergyCoeff[0][0];
    cM[1][1][0][0][nSwap] *= fEnergyCoeff[0][0];
//     std::cout << sMS[nSwap][1][0] << std::endl;
    return fEnergyCoeff[0][0] * fMK0starNR[1][0];
  }
  double GQuarkPairCreationModel::GMDKstar1270(int nSwap)
  {
    cM[0][0][0][1][nSwap] *= fEnergyCoeff[0][0];
    cM[1][0][0][1][nSwap] *= fEnergyCoeff[0][0];
    return fEnergyCoeff[0][0] * fMK0starNR[0][1];
  }
  double GQuarkPairCreationModel::GMDKstar1400(int nSwap)
  {
    cM[0][1][0][1][nSwap] *= fEnergyCoeff[0][0];
    cM[1][1][0][1][nSwap] *= fEnergyCoeff[0][0];
    return fEnergyCoeff[0][0] * fMK0starNR[1][1];
  }
  double GQuarkPairCreationModel::GMSRho1270(int nSwap)
  {
    cM[0][0][1][0][nSwap] *= fEnergyCoeff[0][1];
    cM[1][0][1][0][nSwap] *= fEnergyCoeff[0][1];
    return fEnergyCoeff[0][1] * fMRho0NR[0][0];
  }
  double GQuarkPairCreationModel::GMSRho1400(int nSwap)
  {
    cM[0][1][1][0][nSwap] *= fEnergyCoeff[0][1];
    cM[1][1][1][0][nSwap] *= fEnergyCoeff[0][1];
    return fEnergyCoeff[0][1] * fMRho0NR[1][0];
  }
  double GQuarkPairCreationModel::GMDRho1270(int nSwap)
  {
    cM[0][0][1][1][nSwap] *= fEnergyCoeff[0][1];
    cM[1][0][1][1][nSwap] *= fEnergyCoeff[0][1];
    return fEnergyCoeff[0][1] * fMRho0NR[0][1];
  }
  double GQuarkPairCreationModel::GMDRho1400(int nSwap)
  {
    cM[0][1][1][1][nSwap] *= fEnergyCoeff[0][1];
    cM[1][1][1][1][nSwap] *= fEnergyCoeff[0][1];
    return fEnergyCoeff[0][1] * fMRho0NR[1][1];
  }
}