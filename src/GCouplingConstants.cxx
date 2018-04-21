#include "GCouplingConstants.h"

namespace Gamapola{
  GCouplingConstants::GCouplingConstants(): GInterfaceForMathFunctions()
  {
    std::cout << "GCouplingConstants constructor calling. . ." << std::endl;
  }
  GCouplingConstants::~GCouplingConstants()
  {
    std::cout << "GCouplingConstants destructor calling. . ." << std::endl;
  }
  void GCouplingConstants::GCouplingConstantsCalculation()
  {
    double momKst = GMomentaA( kMK0star_892 * kMK0star_892, kMPion * kMPion, kMKaon);
    double momRho = GMomentaA( kMRho0_775 * kMRho0_775, kMPion * kMPion, kMPion);
    
    fgKstarKpi = sqrt( 3. * 2. * TMath::Pi() * kMK0star_892 * kMK0star_892 * kGammaK0star_892 / ( momKst * momKst * momKst ) );
    fgRhoPiPi = -sqrt( 3. * 2. * TMath::Pi() * kMRho0_775 * kMRho0_775 * kGammaRho0_775 / ( momRho * momRho * momRho ) );
//     double sLow = ( 2 * kMPion + kMKaon ) * ( 2 * kMPion + kMKaon );
//     double sUp = 10.;
    fgK1KappaPi[0] = 43.0666;
//     sqrt(0.28 * kGammaK1_1270 / GCouplingKappa(sLow, sUp));
    fgK1KappaPi[1] = 0;
  }
  double GCouplingConstants::GCouplingKappa(double sLow, double sUp)
  {
    ROOT::Math::Functor1D wf(this,&Gamapola::GCouplingConstants::GOuterIntegral);
    ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR); 
    ig.SetFunction(wf);
    double gKappaKpi = ig.Integral(sLow,sUp);
    double coeff = 1./ ( 8 * TMath::Pi() ) * 1./3. * 4 * kMKappa_800 * kGammaKappa_800 / TMath::Pi() * 
    kMK1_1270 * kGammaK1_1270 / TMath::Pi();
    std::cout << coeff << std::endl;
    return gKappaKpi * coeff;
  }
  double GCouplingConstants::GInnerIntegral(double sKappa)
  {
    double kappaBW2 = 1./ ( ( sKappa - kMKappa_800 * kMKappa_800 ) * ( sKappa - kMKappa_800 * kMKappa_800 )
    + ( kMKappa_800 * kGammaKappa_800 ) * ( kMKappa_800 * kGammaKappa_800 ) );
    double k1270 = 1./ ( ( fs - kMK1_1270 * kMK1_1270 ) * ( fs - kMK1_1270 * kMK1_1270 )
    + ( kMK1_1270 * kGammaK1_1270 ) * ( kMK1_1270 * kGammaK1_1270 ) );
    double mom3 = GMomentaA(fs, sKappa, kMPion) * GMomentaA(fs, sKappa, kMPion) * GMomentaA(fs, sKappa, kMPion);
    return kappaBW2 * k1270 * mom3 / fs;
  }
  double GCouplingConstants::GOuterIntegral(double s)
  {
    ROOT::Math::Functor1D wf(this,&Gamapola::GCouplingConstants::GInnerIntegral);
    ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR); 
    ig.SetFunction(wf);
    fs = s;
    double sKappaLow = ( kMKaon + kMPion ) * ( kMKaon + kMPion );
    double sKappaUp = (sqrt(fs) - kMPion) * (sqrt(fs) - kMPion);
    double val = ig.Integral( sKappaLow , sKappaUp );
    return val;
  }
  
  void GCouplingConstants::GHadronicFormFactor(double s, int nRes)
  {
    double pDecays[2];
    pDecays[0] = sqrt( ( s - (kMK0star_892 + kMPion) * (kMK0star_892 + kMPion) ) * 
                         ( s - (kMK0star_892 - kMPion) * (kMK0star_892 - kMPion) ) ) / ( 2 * sqrt(s) );
    pDecays[1] = sqrt( ( s - (kMRho0_775 + kMKaon) * (kMRho0_775 + kMKaon) ) * 
                        ( s - (kMRho0_775 - kMKaon) * (kMRho0_775 - kMKaon) ) ) / ( 2 * sqrt(s) );    
    for(auto&& iRes = 2; iRes < nRes; iRes++)
      for(size_t iDec = 0; iDec < 2; iDec++)
      {
        fHadronicFF[iRes][iDec] = sqrt( 8 * TMath::Pi() * fResWidths[iRes] * fResBrs[iRes][iDec]
        / ( pDecays[iDec] * pDecays[iDec] * pDecays[iDec] * fCGKstr[iRes] * fCGKstr[iRes] ) );
//         std::cout << fHadronicFF[iRes][iDec] << "  " << iRes << "  " << iDec << "  " << s << std::endl;
      }
//     fHadronicFF[3][0] = sqrt( 8 * TMath::Pi() * fResWidths[3] * fResBrs[3][0]
//     / ( pDecays[0] * pDecays[0] * pDecays[0] * fCGKstr[3] * fCGKstr[3] ) );
//     
//     fHadronicFF[4][0] = sqrt( 8 * TMath::Pi() * fResWidths[4] * fResBrs[4][0]
//     / ( pDecays[0] * pDecays[0] * pDecays[0] * fCGKstr[4] * fCGKstr[4] ) );    
  }
}