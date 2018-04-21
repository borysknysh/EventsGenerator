#include "GMathFunctions.h"

namespace Gamapola{
  GMathFunctions::GMathFunctions():
  GKinematics(), 
  GBreitWigner(), 
  GCouplingConstants(), 
  GQuarkPairCreationModel(),
  GFormFactors(),
  GCurrents()
  {
    std::cout << "GMathFunctions constructor calling. . ." << std::endl;
  }
  GMathFunctions::~GMathFunctions()
  {
    std::cout << "GMathFunctions destructor calling. . ." << std::endl;
  }
  void GMathFunctions::GSetMasses(const double& mi, const double& mj, const double& mk, const double& M)
  {
    fMasses[0] = mi;
    fMasses[1] = mj;
    fMasses[2] = mk;
    fMasses[3] = M;
  }
  void GMathFunctions::GSetNormalizationIntegrals(const std::complex<double>* normIntegrals)
  {
    fNormalizationIntegrals = normIntegrals;
  }
  void GMathFunctions::GSetNumberOfNormalizationIntegrals(const int& nTerms)
  {
    fNIntegrals = nTerms;
  }

  void GMathFunctions::GSetResonances(const std::vector<std::string>& resNames, const std::vector<int>& charges)
  {
    fConstParsForI01.resize(10);
    fConstParsForMomA.resize(10);        
    fMomentaConst.resize(fConstParsForI01.size());
    std::vector<double> pars;
    std::vector<std::vector<double> > parsForRes;
    fResNames = resNames;
    fNRes = resNames.size();
    for(int i = 0; i < fNRes; i++)
    {
      if(resNames[i] == "K1_1270" ||
         resNames[i] == "K1_1400" || 
         resNames[i] == "K*_1410" ||
         resNames[i] == "K*_1680" ||
         resNames[i] == "K2_1430")
      {
        double mK1;
        double gammaK1;
        int iRes = 0;
        
        if(resNames[i] == "K1_1270")
        {
          iRes = 0;
          fIndexControl[iRes] = 1;
        }
        else if(resNames[i] == "K1_1400")
        {
          iRes = 1;
          fIndexControl[iRes] = 1;
        }
        else if(resNames[i] == "K*_1410")
        {
          iRes = 2;
          fIndexControl[iRes] = 1;
        }
        else if(resNames[i] == "K*_1680")
        {
          iRes = 3;
          fIndexControl[iRes] = 1;
        }
        else if(resNames[i] == "K2_1430")
        {
          iRes = 4;
          fIndexControl[iRes] = 1;
        }
        fCharge[iRes] = charges[i];
        fCGKstr[iRes] = (fCharge[iRes]==0)?sqrt(2.)/3.:-2./3.;
        fCGRho[iRes] = (fCharge[iRes]==0)?1./sqrt(3.):-1./sqrt(6.);
        fDelta[iRes] = (fCharge[iRes]==0)?1.:0.;
        
        fDelta2[iRes] = std::complex<double>(1.,0.);
//         fIndex.push_back(iRes);
        fIndex[iRes] = iRes;
        GCouplingConstantsCalculation();
        switch(iRes)
        {
          case 0:
            mK1 = kMK1_1270;
            gammaK1 = kGammaK1_1270;
            fResBrs[iRes][0] = kBr1K1_1270;
            fResBrs[iRes][1] = kBr2K1_1270;
            break;
          case 1:
            mK1 = kMK1_1400;
            gammaK1 = kGammaK1_1400;
            fResBrs[iRes][0] = kBr1K1_1400;
            fResBrs[iRes][1] = kBr2K1_1400;
            break;
          case 2:
            mK1 = kMKst_1410;
            gammaK1 = kGammaKst_1410;
            fResBrs[iRes][0] = kBr1Kst_1410;
            fResBrs[iRes][1] = kBr2Kst_1410;
            break;
          case 3:
            mK1 = kMKst_1680;
            gammaK1 = kGammaKst_1680;
            fResBrs[iRes][0] = kBr1Kst_1680;
            fResBrs[iRes][1] = kBr2Kst_1680;
            break;
          case 4:
            mK1 = kMK2_1430;
            gammaK1 = kGammaK2_1430;
            fResBrs[iRes][0] = kBr1K2_1430;
            fResBrs[iRes][1] = kBr2K2_1430;
            break;
        }
          pars.push_back(kRK1); pars.push_back(kRKst); pars.push_back(kRpi); pars.push_back(kMKaon); pars.push_back(kMPion);
          parsForRes.push_back(pars);
          pars.clear();
          pars.push_back(kRK1); pars.push_back(kRrho); pars.push_back(kRK); pars.push_back(kMPion); pars.push_back(kMPion);
          parsForRes.push_back(pars);
          pars.clear();   
          pars.push_back(kRK1); pars.push_back(kRKst); pars.push_back(kRpi); pars.push_back(kMKaon); pars.push_back(kMPion);
          parsForRes.push_back(pars);
          fConstParsForI01[iRes] = parsForRes;
          pars.clear();
          parsForRes.clear();
          
          pars.push_back(kMPion); pars.push_back(kMK0star_892); pars.push_back(mK1); parsForRes.push_back(pars);
          pars.clear();
          pars.push_back(kMKaon); pars.push_back(kMRho0_775); pars.push_back(mK1); parsForRes.push_back(pars);
          pars.clear();  
          pars.push_back(kMPion); pars.push_back(kMK0star_1430); pars.push_back(mK1);
          parsForRes.push_back(pars);
          
          fConstParsForMomA[iRes] = parsForRes;
          pars.clear();
          parsForRes.clear();

          fResMasses[iRes] = mK1;
          fResWidths[iRes] = gammaK1;
          fMomentaConst[iRes].resize(fConstParsForI01[iRes].size());
      }
    }
    for(auto&& iRes = 0; iRes < fNRes; iRes++)
    {
      if(fIndexControl[iRes] > 0)
      {
        GQPCMConstanstsCalculation(&fConstParsForI01[fIndex[iRes]][0],fIndex[iRes]);
        GKinematicsConstanstsCalculation(iRes);
      }
    }
  }
  double GMathFunctions::GOmegaPDF1()
  {
    return fOmega;
  }
  void GMathFunctions::GCalculateResonance(double* x, double* p, const int& iRes)
  {
    fMomentaA[0] = GMomentaA(x[4], x[2], fConstParsForMomA[iRes][0][0]);
      fEnergyA[0] = GEnergyA(x[4], x[2], fConstParsForMomA[iRes][0][0]);
      fMomentaA[1] = GMomentaA(x[4],fSij, fConstParsForMomA[iRes][1][0]);
      fEnergyA[1] = GEnergyA(x[4], fSij, fConstParsForMomA[iRes][1][0]);
      fEnergyCoeff[0][0] = 8 * sqrt ( TMath::Pi() * TMath::Pi() * TMath::Pi() * 
        GEnergyA(x[4], x[2], fConstParsForMomA[iRes][0][0]) * 
        GEnergyA(x[4], fConstParsForMomA[iRes][0][0] * fConstParsForMomA[iRes][0][0], sqrt(x[2])) * sqrt( x[4] ) );
      
      fEnergyCoeff[0][1] = 8 * sqrt ( TMath::Pi() * TMath::Pi() * TMath::Pi() * 
        GEnergyA(x[4], fSij, fConstParsForMomA[iRes][1][0]) * 
        GEnergyA(x[4], fConstParsForMomA[iRes][1][0] * fConstParsForMomA[iRes][1][0], sqrt(fSij)) * sqrt( x[4] ) );
      
      GCalculateQPMC(p,iRes,0);
      GCalculateFormFactorsKstar(x,p,iRes,0);
      GCalculateFormFactorsRho(x,p,iRes); 
      
      auto sij = x[2];
      fMomentaA[0] = GMomentaA(x[4], fSik, fConstParsForMomA[iRes][0][0]);
      fEnergyA[0] = GEnergyA(x[4], fSik, fConstParsForMomA[iRes][0][0]);
      fEnergyCoeff[0][0] = 8 * sqrt ( TMath::Pi() * TMath::Pi() * TMath::Pi() * 
      GEnergyA(x[4], fSik, fConstParsForMomA[iRes][0][0]) * 
      GEnergyA(x[4], fConstParsForMomA[iRes][0][0] * fConstParsForMomA[iRes][0][0], sqrt(fSik)) * sqrt( x[4] ) );
      
      GCalculateQPMC(p,iRes,1);
      x[2] = fSik;
      GCalculateFormFactorsKstar(x,p,iRes,1);
      x[2] = sij;
  }
  double GMathFunctions::GProcessingComputationOfPDF(double* x, double* p)
  {
    double pdf;
    GCalculateCouplings(x,p);
//     GQPCM
    fSinThK1 = sin(p[2]);
    fCosThK1 = cos(p[2]);
    fSij = x[3];
    fSik = x[4] - x[2] - x[3] + kMKaon * kMKaon + kMPion * kMPion + kMPion * kMPion;
    for(auto iRes = 0; iRes < 2; iRes++)
      if(fIndexControl[iRes] > 0)
        GCalculateResonance(x,p,iRes);
    GCalculateKinematics(x);
    GCalculateHelicityAmplitudes(x,p);
    pdf = GJ2(p,p[7],x[4]);
    pdf = (pdf<0)?pdf:0;
        return pdf;
  }
  std::vector<std::complex<double> > GMathFunctions::GGetKinematicalCoefficients() 
  {
    std::vector<std::complex<double> > kinCoeffs;
    kinCoeffs.resize(240, 0);
    int iii = 0;
//     std::cout << "**********************************************" << std::endl;
    for(size_t iTrig = 0; iTrig < 2; iTrig++) // {sin, cos}
      for(size_t iCurr = 0; iCurr < 2; iCurr++) // C1,2
        for(size_t iRL = 0; iRL < 2; iRL++) // {RL}
          for(size_t iMod = 0; iMod < 3; iMod++) // {K*ro}
            for(size_t iRes = 0; iRes < 5; iRes++) // {K1}
              for(size_t iPh = 0; iPh < 2; iPh++) // {phase S,D}
              {
                if(isnan(std::abs(cRL[iTrig][iRes][iMod][iCurr][iRL][iPh])))
                {
                  std::fill(kinCoeffs.begin(), kinCoeffs.end(), 0);
                  break;
                }
                kinCoeffs[iii] = cRL[iTrig][iRes][iMod][iCurr][iRL][iPh]; // cRL[iTrig][nRes][jC][0][0][iPh]
//                 std::cout << kinCoeffs[iii] << "  " << iTrig << "  " << iRes << "  " << iMod << "  " << iCurr
//                 <<"  " << iRL << "  " << iPh << std::endl;
                iii++;
              }     
    return kinCoeffs;
  }
}