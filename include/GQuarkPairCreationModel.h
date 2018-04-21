#ifndef GQUARKPAIRCREATIONMODEL_H
#define GQUARKPAIRCREATIONMODEL_H

#include "GInterfaceForMathFunctions.h"

namespace Gamapola
{
  class GQuarkPairCreationModel: virtual public GInterfaceForMathFunctions{
  private:
    int fArg;
    double fsqrt2;
    double fsqrt15;
    double sExp[2][fMaxNRes][fMaxNDecays];
    double fcoeff[fMaxNRes][fMaxNDecays];
    double fbracePart[fMaxNRes][fMaxNDecays];
    double fexpPart[fMaxNRes][fMaxNDecays];
    double fcoeffAdd[2];
    double fbracePartAdd[2];
    double fexpPartAdd[2];
  protected:
    GQuarkPairCreationModel();
    virtual ~GQuarkPairCreationModel();
    
    void GCalculateQPMC(double* p,const int& iRes, const int& swap);
    
    double GI0(int nRes, int nFnc);
    double GI1(int nRes, int nFnc);
    double GMS(double* p, int nRes, int nFnc, int nSwap);
    double GMD(double* p, int nRes, int nFnc, int nSwap);
    //     Non-relativistic S and D waves
    double GMSKstarNR1270(int nSwap);
    double GMDKstarNR1270(int nSwap);
    double GMSRhoNR1270(int nSwap);
    double GMDRhoNR1270(int nSwap);
    double GMSKstarNR1400(int nSwap);
    double GMDKstarNR1400( int nSwap);
    double GMSRhoNR1400( int nSwap);
    double GMDRhoNR1400(int nSwap);
//     Relativistic S and D waves
    double GMSKstar1270( int nSwap);
    double GMDKstar1270( int nSwap);
    double GMSRho1270(int nSwap);
    double GMDRho1270(int nSwap);
    double GMSKstar1400(int nSwap);
    double GMDKstar1400(int nSwap);
    double GMSRho1400(int nSwap);
    double GMDRho1400(int nSwap);
    void GQPCMConstanstsCalculation(std::vector<double>* u, int nRes);
    
    double (Gamapola::GQuarkPairCreationModel::*fFunctionPointersKstrNR[fMaxNRes][fMaxNWaves])(int nSwap); 
    double (Gamapola::GQuarkPairCreationModel::*fFunctionPointersKstr[fMaxNRes][fMaxNWaves]) (int nSwap);
    double (Gamapola::GQuarkPairCreationModel::*fFunctionPointersRhoNR[fMaxNRes][fMaxNWaves]) (int nSwap);
    double (Gamapola::GQuarkPairCreationModel::*fFunctionPointersRho[fMaxNRes][fMaxNWaves]) (int nSwap);
    double fI0[fMaxNRes][fMaxNDecays];
    double fI1[fMaxNRes][fMaxNDecays];
    double fMS[fMaxNRes][fMaxNDecays];
    double fMD[fMaxNRes][fMaxNDecays];
    double fEnergyCoeff[fMaxNRes][fMaxNDecays];
    double fMK0starNR[fMaxNRes][fMaxNDecays];
    double fMRho0NR[fMaxNRes][fMaxNDecays];
    double fSinThK1;
    double fCosThK1;
  };
}
#endif