#include "GCurrents.h"

namespace Gamapola{
  GCurrents::GCurrents(): GInterfaceForMathFunctions()
  {
    std::cout << "GCurrents constructor calling. . ." << std::endl;
    fFF[2][0] = kFF1Kst_1410;
    fFF[2][1] = kFF2Kst_1410;
    fFF[3][0] = kFF1Kst_1680;
    fFF[3][1] = kFF2Kst_1680;
    fFF[4][0] = kFF1K2_1430;
    fFF[4][1] = kFF2K2_1430;
        
    fLeft[0] = 1;
    fLeft[1] = (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
    fLeft[2] = (kMBmeson * kMBmeson - kMKst_1410 * kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
    fLeft[3] = (kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
    fLeft[4] = 1;
        
    fRight[0] =  1;
    fRight[1] =  (kMBmeson * kMBmeson - kMK1_1400 *  kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
    fRight[2] = -(kMBmeson * kMBmeson - kMKst_1410 *  kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
    fRight[3] = -(kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
    fRight[4] = -1;
  }
  GCurrents::~GCurrents()
  {
    std::cout << "GCurrents destructor calling. . ." << std::endl;
  }
  void GCurrents::GCalculateHelicityAmplitudes(double* x,double* p)
  {
    fc1[0][fCharge[0]] = GC1(p,x[4], x[2], fSik, 0) /** fMomi*/;
    fc2[0][fCharge[0]] = GC2(p,x[4], x[2], fSik, 0) /** fMomj*/;
    fc3[0][fCharge[2]] = GC3(x[4], fSik, 2);
    fc1[1][fCharge[1]] = GC1(p,x[4], x[2], fSik, 1) /** fMomi*/;
    fc2[1][fCharge[1]] = GC2(p,x[4], x[2], fSik, 1) /** fMomj*/;
    fc3[1][fCharge[3]] = GC3(x[4], fSik, 3);
    fc4Rho = GC4Rho(x[4], 4);
    fc4Kstr = GC4Kst(x[4], x[2], 4);
    fc5Kstr = GC5Kst(x[4], fSik, 4);
    fc4[0][fCharge[4]] = GC4(x[4], 4);
    fc5[0][fCharge[4]] = GC5(x[4], 4);
  }
  void GCurrents::GCalculateCouplings(double*x,double* p)
  {
    fDelta2[0] = std::complex<double>(p[29],0);
    fDelta2[1] = std::complex<double>(p[8],p[9]);
    fDelta2[2] = std::complex<double>(p[11],p[12]);
    fDelta2[3] = std::complex<double>(p[13],p[14]);
    fDelta2[4] = std::complex<double>(p[15],p[16]);
    fCCouplings[2][0] = std::complex<double>(p[17],p[18]);
    fCCouplings[2][1] = std::complex<double>(p[19],p[20]);
    fCCouplings[3][0] = std::complex<double>(p[21],p[22]);
    fCCouplings[3][1] = std::complex<double>(p[23],p[24]);
    fCCouplings[4][0] = std::complex<double>(p[25],p[26]);
    fCCouplings[4][1] = std::complex<double>(p[27],p[28]);
    fLeft[4] = (kMBmeson * sqrt(2 * x[4])) / (kMBmeson * kMBmeson - x[4]) ;
    fRight[4] = -fLeft[4]; 
  }
  std::complex<double> GCurrents::GC1Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sSmallBrace = (1 + ( kMPion * kMPion - kMKaon * kMKaon ) / sik );
    double sBigBrace = ( 2 * fPiPjVec4 - sSmallBrace * ( sqrt(s) * fEnergy1 - kMPion * kMPion ) );
    std::complex<double> c1Kstr = fCGKstr[nRes] * fgKstarKpi * 
    ( (fA[nRes][0][1] * sSmallBrace + fB[nRes][0][1] * sBigBrace) * GBWKstr(sik) - 2. * fA[nRes][0][0] * GBWKstr(sij) * fDelta[nRes] );
//     std::cout << c1Kstr << std::endl;
    for(int iphSKstr = 0; iphSKstr < 2; iphSKstr++)
    {
      cC[0][nRes][0][0][iphSKstr] = fCGKstr[nRes] * fgKstarKpi * 
      ( (cF[0][nRes][0][0][1][iphSKstr] * sSmallBrace + cF[0][nRes][0][1][1][iphSKstr] * sBigBrace) * GBWKstr(sik)
      - 2. * cF[0][nRes][0][0][0][iphSKstr] * GBWKstr(sij) * fDelta[nRes] );
    
      cC[1][nRes][0][0][iphSKstr] = fCGKstr[nRes] * fgKstarKpi * 
      ( (cF[1][nRes][0][0][1][iphSKstr] * sSmallBrace + cF[1][nRes][0][1][1][iphSKstr] * sBigBrace) * GBWKstr(sik)
      - 2. * cF[1][nRes][0][0][0][iphSKstr] * GBWKstr(sij) * fDelta[nRes] );
    }
//     std::cout << "C1 of K*" << c1Kstr << "  " << nRes << std::endl;
//     std::complex<double> c1Kstr = (cC[0][nRes][0][0][0] * sin(p[2]) + cC[1][nRes][0][0][0] * cos(p[2])) + 
//     (cC[0][nRes][0][0][1] * sin(p[2]) + cC[1][nRes][0][0][1] * cos(p[2])) * expPhaseSKstar;
    return c1Kstr;
  }
  std::complex<double> GCurrents::GC2Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sSmallBrace = (1 + ( kMPion * kMPion - kMKaon * kMKaon ) / sij );
    double sBigBrace = ( 2 * fPiPjVec4 - sSmallBrace * ( sqrt(s) * fEnergy2 - kMPion * kMPion ) );
    std::complex<double> c2Kstr = fCGKstr[nRes] * fgKstarKpi * 
    ( (fA[nRes][0][0] * sSmallBrace + fB[nRes][0][0] * sBigBrace) * GBWKstr(sij) * fDelta[nRes]
    - 2. * fA[nRes][0][1] * GBWKstr(sik) );
//     std::cout << c2Kstr << std::endl;
    for(int iphSKstr = 0; iphSKstr < 2; iphSKstr++)
    {
      cC[0][nRes][0][1][iphSKstr] = fCGKstr[nRes] * fgKstarKpi * 
      ( (cF[0][nRes][0][0][0][iphSKstr] * sSmallBrace + cF[0][nRes][0][1][0][iphSKstr] * sBigBrace) * GBWKstr(sij) * fDelta[nRes]
      - 2. * cF[0][nRes][0][0][1][iphSKstr] * GBWKstr(sik) );
      
      cC[1][nRes][0][1][iphSKstr] = fCGKstr[nRes] * fgKstarKpi * 
      ( (cF[1][nRes][0][0][0][iphSKstr] * sSmallBrace + cF[1][nRes][0][1][0][iphSKstr] * sBigBrace) * GBWKstr(sij) * fDelta[nRes]
      - 2. * cF[1][nRes][0][0][1][iphSKstr] * GBWKstr(sik) );
    }
//     std::cout << "C2 of K*" << c2Kstr << "  " << nRes << std::endl;
//     std::complex<double> c2Kstr = (cC[0][nRes][0][1][0] * sin(p[2]) + cC[1][nRes][0][1][0] * cos(p[2])) + 
//     (cC[0][nRes][0][1][1] * sin(p[2]) + cC[1][nRes][0][1][1] * cos(p[2])) * expPhaseSKstar;
    return c2Kstr;
  }
  std::complex<double> GCurrents::GC1Rho(double* p, double s, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sVerySmallBrace = sqrt( s ) * ( fEnergy1 - fEnergy2 );
    std::complex<double> c1Rho = fCGRho[nRes] * fgRhoPiPi * ( fA[nRes][1][0] - fB[nRes][1][0] * sVerySmallBrace) * GBWRho(fSij);
//     std::cout << c1Rho << std::endl;
    for(int iphRho = 0; iphRho < 2; iphRho++)
    {
      cC[0][nRes][1][0][iphRho] = fCGRho[nRes] * fgRhoPiPi * ( cF[0][nRes][1][0][0][iphRho] - 
      cF[0][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
      
      cC[1][nRes][1][0][iphRho] = fCGRho[nRes] * fgRhoPiPi * ( cF[1][nRes][1][0][0][iphRho] - 
      cF[1][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
    }
//     std::cout << "C1 of Rho" << c1Rho << "  " << nRes << std::endl;
//     std::complex<double> c1Rho = (cC[0][nRes][1][0][0] * sin(p[2]) + cC[1][nRes][1][0][0] * cos(p[2])) * expPhaseDRho + 
//     (cC[0][nRes][1][0][1] * sin(p[2]) + cC[1][nRes][1][0][1] * cos(p[2])) * expPhaseDRho * expPhaseSRho;
    return c1Rho;
  }
  std::complex<double> GCurrents::GC2Rho(double* p, double s, int nRes)
  {
   std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sVerySmallBrace = sqrt( s ) * ( fEnergy1 - fEnergy2 );
    std::complex<double> c2Rho = fCGRho[nRes] * fgRhoPiPi * ( fA[nRes][1][0] + fB[nRes][1][0] * sVerySmallBrace) * GBWRho(fSij);
//     std::cout << c2Rho << std::endl;
    for(int iphRho = 0; iphRho < 2; iphRho++)
    {
      cC[0][nRes][1][1][iphRho] = fCGRho[nRes] * fgRhoPiPi * ( cF[0][nRes][1][0][0][iphRho] + 
      cF[0][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
      
      cC[1][nRes][1][1][iphRho] = fCGRho[nRes] * fgRhoPiPi * ( cF[1][nRes][1][0][0][iphRho] + 
      cF[1][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
    }
//     std::complex<double> c2Rho = (cC[0][nRes][1][1][0] * sin(p[2]) + cC[1][nRes][1][1][0] * cos(p[2])) * expPhaseDRho + 
//     (cC[0][nRes][1][1][1] * sin(p[2]) + cC[1][nRes][1][1][1] * cos(p[2])) * expPhaseDRho * expPhaseSRho;
//     std::cout << "C2 of Rho" << c2Rho << "  " << nRes << std::endl;
    return c2Rho;
  }
  std::complex<double> GCurrents::GC1Kappa(double* p, double sik, int nRes)
  {
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> c1Kappa = -fCGKstr[nRes] * 2 * fgK1KappaPi[nRes] * GBWKappa(sik);
    for(int iph = 0; iph < 2; iph++)
    {
      cC[0][nRes][2][0][iph] = c1Kappa;
      cC[1][nRes][2][0][iph] = c1Kappa;
    }
//      std::cout << " C1Kappa " << c1Kappa << "  " << nRes << "  " << fDelta[nRes] << "   " << fCGKstr[nRes] << std::endl;
    return ff2 * c1Kappa;
  }
  std::complex<double> GCurrents::GC2Kappa(double* p, double sij, int nRes)
  {
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> c2Kappa = -fCGKstr[nRes] * 2 * fgK1KappaPi[nRes] * GBWKappa(sij) * fDelta[nRes];
    for(int iph = 0; iph < 2; iph++)
    {
      cC[0][nRes][2][1][iph] = c2Kappa;
      cC[1][nRes][2][1][iph] = c2Kappa;
    }
//     std::cout << " C2Kappa " << c2Kappa << "  " << nRes << "  " << fDelta[nRes] << "   " << fCGKstr[nRes] << std::endl;
    return ff2 * c2Kappa;
  }
  std::complex<double> GCurrents::GC3Kst(double s, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0,-2.) * sqrt(s) * fFF[nRes][0];
    std::complex<double> c3Kstr = fCGKstr[nRes] * fgKstarKpi * GBWKstr(sik) * coeff;
    cC[0][nRes][0][0][0] = c3Kstr;
//     std::cout << "C3 gg3Kstr " << fFF[nRes][0] << "   "<< fCCouplings[nRes][0] << std::endl;
//     std::cout << std::complex<double>(0,-2.) * sqrt(s) * fFF[nRes][0] * fgKstarKpi * GBWKstr(sij) * fCGKstr[nRes] << std::endl;
//     std::cout << fCCouplings[nRes][0] * c3Kstr * fPiPjVec[0] * fEpsilonL[0] << std::endl;
//     std::cout << fCGKstr[nRes] << "  " << fgKstarKpi << "   " << GBWKstr(sij) << "   " << fFF[nRes][0] << std::endl;
    return fCCouplings[nRes][0] * c3Kstr;
  }
  std::complex<double> GCurrents::GC3Rho(double s, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0,2.) * sqrt(s) * fFF[nRes][1];
    std::complex<double> c3Rho = fCGRho[nRes] * fgRhoPiPi * GBWRho(fSij) * coeff;
    cC[0][nRes][1][0][0] = c3Rho;
//     std::cout << "C3 gg3Rho " << fFF[nRes][1] << "  " << fCCouplings[nRes][1] << std::endl;
//     std::cout << fCCouplings[nRes][1] * c3Rho * fPiPjVec[0] * fEpsilonL[0] << std::endl;
    return fCCouplings[nRes][1] * c3Rho;
  }
  std::complex<double> GCurrents::GC4Kst(double s, double sij, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0, sqrt(8.)) * fgKstarKpi * sqrt(s) * fFF[nRes][0];
    std::complex<double> c4Kstr = fCGKstr[nRes] * GBWKstr(sij) * coeff * fDelta[nRes];
//     std::cout << "C4 gg5Kstr " << fFF[nRes][0] << "  " << fCCouplings[nRes][0] << std::endl;
    cC[0][nRes][0][0][0] = c4Kstr;
    
    return fCCouplings[nRes][0] * c4Kstr;
  }
  std::complex<double> GCurrents::GC5Kst(double s, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0, sqrt(8.)) * fgKstarKpi * sqrt(s) * fFF[nRes][0];
    std::complex<double> c5Kstr = fCGKstr[nRes] * GBWKstr(sik) * coeff;
//     std::cout << "C5 gg5Kstr " << fFF[nRes][0] << "  " << fCCouplings[nRes][0] << std::endl;
    cC[0][nRes][0][1][0] = c5Kstr;
    
    return fCCouplings[nRes][0] * c5Kstr;
  }
  std::complex<double> GCurrents::GC4Rho(double s, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0, sqrt(8.)) * fgRhoPiPi * sqrt(s) * fFF[nRes][1];
    std::complex<double> c4Rho = fCGRho[nRes] * GBWRho(fSij) * coeff;
    cC[0][nRes][1][0][0] = c4Rho;
    cC[0][nRes][1][1][0] = c4Rho;
//     std::cout << "C4 gg5Rho " << fFF[nRes][1] << "  " << fCCouplings[nRes][1] << std::endl;
    return fCCouplings[nRes][1] * c4Rho;
  }
  
  std::complex<double> GCurrents::GC1(double* p, double s, double sij, double sik, int nRes)
  {
   std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> currentC1 = ( GC1Kst(p,s, sij, sik, nRes) + GC1Rho(p,s, nRes) + 
    GC1Kappa(p,sik, nRes) ) * GBWK1(s, nRes);
//     std::cout <<"C1: " << GC1Kst(p,s, sij, sik, charge, nRes) * GBWK1(s, nRes) << "  " <<
//     GC1Rho(p,s, nRes) * GBWK1(s, nRes) << "  " << 
//      GC1Kappa(p,sik,nRes) * GBWK1(s, nRes) << std::endl;
         
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 3; j++)
        for(int k = 0; k < 2; k++)
          cC[i][nRes][j][0][k] *= GBWK1(s, nRes)/sqrt(s);
    
//     std::cout << 
//     ((cC[0][nRes][0][0][0] * sin(p[2]) + cC[1][nRes][0][0][0] * cos(p[2])) + 
//     (cC[0][nRes][0][0][1] * sin(p[2]) + cC[1][nRes][0][0][1] * cos(p[2])) * expPhaseSKstar + 
//     (cC[0][nRes][1][0][0] * sin(p[2]) + cC[1][nRes][1][0][0] * cos(p[2])) * expPhaseDRho + 
//     (cC[0][nRes][1][0][1] * sin(p[2]) + cC[1][nRes][1][0][1] * cos(p[2])) * expPhaseDRho * expPhaseSRho + 
//     ff2 * cC[0][nRes][2][0][0])/* * GBWK1(s, nRes)*/ << "   " << currentC1/sqrt(s) << std::endl;
    return currentC1;
  }
  std::complex<double> GCurrents::GC2(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> currentC2 = ( GC2Kst(p,s, sij, sik, nRes) + GC2Rho(p,s, nRes) +
    GC2Kappa(p,sij, nRes) ) * GBWK1(s, nRes);
//     std::cout << "C2 :" << GC2Kst(p,s, sij, sik, charge, nRes) * GBWK1(s, nRes) << "  " <<
//     GC2Rho(p,s, nRes) * GBWK1(s, nRes) << "  " << 
//      GC2Kappa(p,sij, nRes) * GBWK1(s, nRes) << std::endl;
     
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 3; j++)
        for(int k = 0; k < 2; k++)
          cC[i][nRes][j][1][k] *= GBWK1(s, nRes)/sqrt(s);
     
//     std::cout << ((cC[0][nRes][0][1][0] * sin(p[2]) + cC[1][nRes][0][1][0] * cos(p[2])) + 
//     (cC[0][nRes][0][1][1] * sin(p[2]) + cC[1][nRes][0][1][1] * cos(p[2])) * expPhaseSKstar + 
//     (cC[0][nRes][1][1][0] * sin(p[2]) + cC[1][nRes][1][1][0] * cos(p[2])) * expPhaseDRho + 
//     (cC[0][nRes][1][1][1] * sin(p[2]) + cC[1][nRes][1][1][1] * cos(p[2])) * expPhaseDRho * expPhaseSRho + 
//     ff2 * cC[0][nRes][2][1][0]) << "  " << currentC2/sqrt(s) << std::endl;
    return currentC2;
  }
  std::complex<double> GCurrents::GC3(double s, double sik, int nRes)
  {
    std::complex<double> currentC3 = (GC3Kst(s, sik, nRes) + GC3Rho(s,nRes)) * GBWK1(s, nRes);
    
      for(int j = 0; j < 2; j++)
        cC[0][nRes][j][0][0] *= GBWK1(s, nRes)/sqrt(s); // remove sqrt when checking
      
//      std::cout << (fRho * cC[0][nRes][1][0][0] + fK * cC[0][nRes][0][0][0]) /** sqrt(s)*/ <<  "   " << currentC3 / sqrt(s) << "  " << nRes << std::endl;
     return currentC3;
  }
  
  std::complex<double> GCurrents::GC4(double s, int nRes)
  {
    std::complex<double> currentC4 = (fc4Kstr + fc4Rho) * GBWK1(s, nRes);
//     std::cout << " 1_4: "<< fc4Rho * GBWK1(s, nRes) / sqrt(s) << std::endl;
      for(int j = 0; j < 2; j++)
        cC[0][nRes][j][0][0] *= GBWK1(s, nRes)/sqrt(s);// remove sqrt when checking
//     std::cout << " 2_4: "<< cC[0][nRes][1][0][0] << std::endl;     
//     std::cout << cC[0][nRes][0][0][0] + cC[0][nRes][1][0][0] << "  " << currentC4 / sqrt(s) << std::endl;
    return currentC4;
  }
  std::complex<double> GCurrents::GC5(double s, int nRes)
  {
    std::complex<double> currentC5 = (fc5Kstr + fc4Rho) * GBWK1(s, nRes);
//     std::cout << " 1_5: "<< fc4Rho * GBWK1(s, nRes) / sqrt(s) << std::endl;
      for(int j = 0; j < 2; j++)
        cC[0][nRes][j][1][0] *= GBWK1(s, nRes)/sqrt(s);// remove sqrt when checking
//     std::cout << " 2_5: "<< cC[0][nRes][1][1][0] << std::endl;   
//     std::cout << cC[0][nRes][0][1][0] + cC[0][nRes][1][1][0] << "  " << currentC5/sqrt(s) << std::endl;
    return currentC5;
  }
    
  std::complex<double>* GCurrents::GJVec(double* p, int charge, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
//     std::cout << fDelta2[nRes] << "  " << nRes << std::endl;
    for(int iTrig = 0; iTrig < 2; iTrig++)
      for(int jC = 0; jC < 3; jC++)
        for(int iPh = 0; iPh < 2; iPh++)
        {
          
          cRL[iTrig][nRes][jC][0][0][iPh] = (fEpsilonR[0] * fMomiV[0] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonR[1] * fMomiV[1] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonR[2] * fMomiV[2] * cC[iTrig][nRes][jC][0][iPh]) * fRight[nRes];
        
          cRL[iTrig][nRes][jC][1][0][iPh] = (-fEpsilonR[0] * fMomjV[0] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonR[1] * fMomjV[1] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonR[2] * fMomjV[2] * cC[iTrig][nRes][jC][1][iPh]) * fRight[nRes];
          
          cRL[iTrig][nRes][jC][0][1][iPh] = (fEpsilonL[0] * fMomiV[0] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonL[1] * fMomiV[1] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonL[2] * fMomiV[2] * cC[iTrig][nRes][jC][0][iPh]) * fLeft[nRes];
          
          cRL[iTrig][nRes][jC][1][1][iPh] = (-fEpsilonL[0] * fMomjV[0] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonL[1] * fMomjV[1] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonL[2] * fMomjV[2] * cC[iTrig][nRes][jC][1][iPh]) * fLeft[nRes];
        }
        
//         for(size_t iCurr = 0; iCurr < 2; iCurr++)
//           for(size_t iRL = 0; iRL < 2; iRL++)
//             for(size_t iMod = 0; iMod < 2; iMod++)
//               for(size_t iRes = 0; iRes < 2; iRes++)
//                 for(size_t iPh = 0; iPh < 2; iPh++)
//                   std::cout << cRL[0][iRes][2][iCurr][iRL][iPh]  << std::endl;
//         std::complex<double> rl[2]; 
//         for(size_t i3 = 0; i3 < 2; i3++) // {RL}
//         {// {sin, cos} {C1,2}, {RL}, {K*ro},  {K1}, {phase S,D}
//           rl[i3] = 
//                ( ((cRL[0][0][0][0][i3][0] * sin(p[2]) + cRL[1][0][0][0][i3][0] * cos(p[2])) + expPhaseSKstar * 
//                 (cRL[0][0][0][0][i3][1] * sin(p[2]) + cRL[1][0][0][0][i3][1] * cos(p[2]))) +  
//                 ((cRL[0][0][1][0][i3][0] * sin(p[2]) + cRL[1][0][1][0][i3][0] * cos(p[2])) * expPhaseSRho + expPhaseSRho * expPhaseDRho *
//                 (cRL[0][0][1][0][i3][1] * sin(p[2]) + cRL[1][0][1][0][i3][1] * cos(p[2])) ) + ff2 * cRL[0][0][2][0][i3][0] +
//                 ((cRL[0][0][0][1][i3][0] * sin(p[2]) + cRL[1][0][0][1][i3][0] * cos(p[2])) + expPhaseSKstar * 
//                 (cRL[0][0][0][1][i3][1] * sin(p[2]) + cRL[1][0][0][1][i3][1] * cos(p[2]))) +  
//                 ((cRL[0][0][1][1][i3][0] * sin(p[2]) + cRL[1][0][1][1][i3][0] * cos(p[2])) * expPhaseSRho + expPhaseSRho * expPhaseDRho *
//                 (cRL[0][0][1][1][i3][1] * sin(p[2]) + cRL[1][0][1][1][i3][1] * cos(p[2])) ) + ff2 * cRL[0][0][2][1][i3][0] )
//                + std::complex<double>(p[8],p[10]) * 
//                ( ((cRL[0][1][0][0][i3][0] * sin(p[2]) + cRL[1][1][0][0][i3][0] * cos(p[2])) + expPhaseSKstar * 
//                 (cRL[0][1][0][0][i3][1] * sin(p[2]) + cRL[1][1][0][0][i3][1] * cos(p[2]))) +  
//                 ((cRL[0][1][1][0][i3][0] * sin(p[2]) + cRL[1][1][1][0][i3][0] * cos(p[2])) * expPhaseSRho + expPhaseSRho * expPhaseDRho *
//                 (cRL[0][1][1][0][i3][1] * sin(p[2]) + cRL[1][1][1][0][i3][1] * cos(p[2])) ) + ff2 * cRL[0][1][2][0][i3][0] +
//                 ((cRL[0][1][0][1][i3][0] * sin(p[2]) + cRL[1][1][0][1][i3][0] * cos(p[2])) + expPhaseSKstar * 
//                 (cRL[0][1][0][1][i3][1] * sin(p[2]) + cRL[1][1][0][1][i3][1] * cos(p[2]))) +  
//                 ((cRL[0][1][1][1][i3][0] * sin(p[2]) + cRL[1][1][1][1][i3][0] * cos(p[2])) * expPhaseSRho + expPhaseSRho * expPhaseDRho *
//                 (cRL[0][1][1][1][i3][1] * sin(p[2]) + cRL[1][1][1][1][i3][1] * cos(p[2])) ) + ff2 * cRL[0][1][2][1][i3][0] );
//         }
//         std::cout << rl[0] << "   " << rl[1] << std::endl;
//         std::cout << (rl[0] * std::conj(rl[0]) + rl[1] * std::conj(rl[1]) + 
//    p[7] * (rl[0] * std::conj(rl[0]) - rl[1] * std::conj(rl[1]))) << std::endl;
//     std::cout << fDelta2[nRes] * (((cRL[0][nRes][0][0][1] * sin(p[2]) + cRL[1][nRes][0][0][1] * cos(p[2]))
//                          + p[6] * (cRL[0][nRes][1][0][1] * sin(p[2]) + cRL[1][nRes][1][0][1] * cos(p[2])))
//                          +        ((cRL[0][nRes][0][1][1] * sin(p[2]) + cRL[1][nRes][0][1][1] * cos(p[2]))
//                          + p[6] * (cRL[0][nRes][1][1][1] * sin(p[2]) + cRL[1][nRes][1][1][1] * cos(p[2])))) << std::endl;  
    for(size_t iAxis = 0; iAxis < 3; iAxis++)
      fJVec3[nRes][charge][iAxis] = fDelta2[nRes] * (fc1[nRes][charge] * fMomiV[iAxis] - fc2[nRes][charge] * fMomjV[iAxis]); 
    
//     std::cout << (fJVec3[nRes][charge][0] * fEpsilonL[0] + fJVec3[nRes][charge][1] * fEpsilonL[1] + fJVec3[nRes][charge][2] * fEpsilonL[2]) << "  " <<
//                  (fJVec3[nRes][charge][0] * fEpsilonR[0] + fJVec3[nRes][charge][1] * fEpsilonR[1] + fJVec3[nRes][charge][2] * fEpsilonR[2]) << "  " << nRes << std::endl;
    return fJVec3[nRes][charge];
  }
  std::complex<double>* GCurrents::GLVec(int charge, int nRes)
  {
   for(size_t iAxis = 0; iAxis < 3; iAxis++)
      fJVec3[nRes][charge][iAxis] = fDelta2[nRes] * fc3[nRes-2][charge] * fPiPjVec[iAxis];
   for(int jC = 0; jC < 2; jC++)
   {
          
    cRL[0][nRes][jC][0][0][0] = cC[0][nRes][jC][0][0] * (fEpsilonR[0] * fPiPjVec[0] + fEpsilonR[1] * fPiPjVec[1]
    + fEpsilonR[2] * fPiPjVec[2] ) * fRight[nRes];
        
    cRL[0][nRes][jC][0][1][0] = cC[0][nRes][jC][0][0] * (fEpsilonL[0] * fPiPjVec[0] + fEpsilonL[1] * fPiPjVec[1]
    + fEpsilonL[2] * fPiPjVec[2]) * fLeft[nRes];
//     std::cout << cRL[0][nRes][jC][0][0][0] << "   "  << 0 << "  " << 0 << "  " << jC << "  " << nRes << std::endl;
//     std::cout << cRL[0][nRes][jC][0][1][0] << "   "  << 0 << "  " << 1 << "  " << jC << "  " << nRes << std::endl;
   }
//    std::cout << (cRL[0][nRes][0][0][1][0] + cRL[0][nRes][1][0][1][0]) * fDelta2[nRes] <<
//    "   " << 
//    fJVec3[nRes][charge][0] * fEpsilonL[0] + fJVec3[nRes][charge][1] * fEpsilonL[1] + fJVec3[nRes][charge][2] * fEpsilonL[2]  << std::endl;

//    for(size_t iCurr = 0; iCurr < 2; iCurr++) // C1,2
//         for(size_t iRL = 0; iRL < 2; iRL++) // {RL}
//           std::cout << fDelta2[nRes] * (cRL[0][nRes][0][0][iRL][0] + cRL[0][nRes][1][0][iRL][0]) << "  " << iRL << "  " << nRes << std::endl;
    return fJVec3[nRes][charge];
  }
  std::complex<double>* GCurrents::GKVec(int charge, int nRes)
  {
    for(size_t iAxis = 0; iAxis < 3; iAxis++)
      fJVec3[nRes][charge][iAxis] = fDelta2[nRes] * (fc4[nRes-4][charge] * fK2KinC4combo[iAxis] +
                                                     fc5[nRes-4][charge] * fK2KinC5combo[iAxis] );
    for(int jC = 0; jC < 2; jC++)
        {
          cRL[0][nRes][jC][0][0][0] = cC[0][nRes][jC][0][0] * 
          (fEpsilonR[0] * fK2KinC4combo[0] + fEpsilonR[1] * fK2KinC4combo[1] + fEpsilonR[2] * fK2KinC4combo[2]) * fRight[nRes];
        
          cRL[0][nRes][jC][1][0][0] = cC[0][nRes][jC][1][0] * 
          (fEpsilonR[0] * fK2KinC5combo[0] + fEpsilonR[1] * fK2KinC5combo[1] + fEpsilonR[2] * fK2KinC5combo[2]) * fRight[nRes];
          
          cRL[0][nRes][jC][0][1][0] = cC[0][nRes][jC][0][0] * 
          (fEpsilonL[0] * fK2KinC4combo[0] + fEpsilonL[1] * fK2KinC4combo[1] + fEpsilonL[2] * fK2KinC4combo[2]) * fLeft[nRes];
        
          cRL[0][nRes][jC][1][1][0] = cC[0][nRes][jC][1][0] * 
          (fEpsilonL[0] * fK2KinC5combo[0] + fEpsilonL[1] * fK2KinC5combo[1] + fEpsilonL[2] * fK2KinC5combo[2]) * fLeft[nRes];
//           std::cout << cRL[0][nRes][jC][0][0][0] << "   "  << 0 << "  " << 0 << "  " << jC << "  " << nRes << std::endl;
//           std::cout << cRL[0][nRes][jC][1][0][0] << "   "  << 1 << "  " << 0 << "  " << jC << "  " << nRes << std::endl;
//           std::cout << cRL[0][nRes][jC][0][1][0] << "   "  << 0 << "  " << 1 << "  " << jC << "  " << nRes << std::endl;
//           std::cout << cRL[0][nRes][jC][1][1][0] << "   "  << 1 << "  " << 1 << "  " << jC << "  " << nRes << std::endl;
        }
//        for(size_t iRL = 0; iRL < 2; iRL++) // {RL}
//           std::cout << fDelta2[nRes] * (cRL[0][nRes][0][0][iRL][0] + cRL[0][nRes][1][0][iRL][0] + 
//             cRL[0][nRes][0][1][iRL][0] + cRL[0][nRes][1][1][iRL][0]) << "  " << iRL << "  " << nRes << std::endl;
            
//    std::cout << fDelta2[nRes] * (cRL[0][nRes][0][0][0][0] + cRL[0][nRes][1][0][0][0] + cRL[0][nRes][0][1][0][0] + cRL[0][nRes][1][1][0][0]) << 
//    "  " << fJVec3[nRes][charge][0] * fEpsilonR[0] + fJVec3[nRes][charge][1] * fEpsilonR[1] + fJVec3[nRes][charge][2] * fEpsilonR[2] << std::endl;
    return fJVec3[nRes][charge];
  }
  double GCurrents::GJ2(double* p, double lambda, double s)
  { 
    std::complex<double> matrixElementLX(0.,0.);
    std::complex<double> matrixElementLY(0.,0.);
    std::complex<double> matrixElementLZ(0.,0.);
    
    std::complex<double> matrixElementRX(0.,0.);
    std::complex<double> matrixElementRY(0.,0.);
    std::complex<double> matrixElementRZ(0.,0.);
    GJVec(p,fCharge[0], fIndex[0]);
    GJVec(p,fCharge[1], fIndex[1]);
    GLVec(fCharge[2], fIndex[2]);
    GLVec(fCharge[3], fIndex[3]);
    GKVec(fCharge[4], fIndex[4]);
    for(auto&& nRes = 0; nRes < 5; nRes++)
    {
      matrixElementLX += fLeft[nRes] * fJVec3[fIndex[nRes]][fCharge[nRes]][0];
      matrixElementLY += fLeft[nRes] * fJVec3[fIndex[nRes]][fCharge[nRes]][1];
      matrixElementLZ += fLeft[nRes] * fJVec3[fIndex[nRes]][fCharge[nRes]][2];
      
      matrixElementRX += fRight[nRes] * fJVec3[fIndex[nRes]][fCharge[nRes]][0];
      matrixElementRY += fRight[nRes] * fJVec3[fIndex[nRes]][fCharge[nRes]][1];
      matrixElementRZ += fRight[nRes] * fJVec3[fIndex[nRes]][fCharge[nRes]][2];
//       std::cout << fCCouplings[nRes][0] << "  " << fCCouplings[nRes][1] << "  " << fLeft[nRes] << "  " << fRight[nRes] << "  "<< nRes << std::endl;
//       std::cout << fJVec3[nRes][fCharge[nRes]][0] << "   " << fJVec3[nRes][fCharge[nRes]][1] <<
//       "   " << fJVec3[nRes][fCharge[nRes]][2] << "  " << nRes << std::endl;
    }
    std::complex<double> scprL = matrixElementLX * fEpsilonL[0] + matrixElementLY * fEpsilonL[1] + matrixElementLZ * fEpsilonL[2];
    std::complex<double> scprR = matrixElementRX * fEpsilonR[0] + matrixElementRY * fEpsilonR[1] + matrixElementRZ * fEpsilonR[2];
//     std::cout << scprR << "  " << scprL << std::endl;
    double pdf = -1e-6*( ( 1 - lambda ) * std::real( scprL * std::conj( scprL ) ) + 
                         ( 1 + lambda ) * std::real( scprR * std::conj( scprR ) ) ) / s; // could be...
    fOmega = ( -std::real( scprL * std::conj( scprL ) ) + std::real( scprR * std::conj( scprR ) ) ) / 
             ( std::real( scprR * std::conj( scprR ) ) + std::real( scprL * std::conj( scprL ) ) );
             
//     std::cout << pdf << std::endl;
    return pdf;
//     std::cout << matrixElementX << "  " << matrixElementY << "  " << matrixElementZ << std::endl;
//       return -1e-5*std::real(matrixElementX * std::conj(matrixElementX) + 
//              matrixElementY * std::conj(matrixElementY) + 
//              matrixElementZ * std::conj(matrixElementZ)) / s;
  }
}