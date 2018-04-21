#include "GEventsGenerator.h"
#include "TTree.h"
#include "TFile.h"

namespace Gamapola{
  GEventsGenerator::GEventsGenerator():GInterfaceForMinimization(),
  fArg(0),
  fNEvents(0),
  fLowCosThetaLimit(0),
  fUpCosThetaLimit(0),
  fLowPhiLimit(0),
  fUpPhiLimit(0),
//     generated kinematic variables
  fCosThetaPDF(0),
  fPhiPDF(0),
  fPDF1(0),
  fOmega(0),
  fOmega2(0),
//     generated flat kinematic variables
  fCosThetaPDFFlat(0),
  fPhiPDFFlat(0),
  fPDF1Flat(0),
  fOmegaFlat(0),
  fOmega2Flat(0)
  {
    std::cout << "GEventsGenerator constructor calling. . ." << std::endl;
  }
  GEventsGenerator::~GEventsGenerator()
  {
//     fNormalizationIntegrals = NULL;
    delete fNormalizationIntegrals;
    
    fCosThetaPDFFlat.clear();
    fPhiPDFFlat.clear();
    fPDF1Flat.clear();
    fOmegaFlat.clear();
    fOmega2Flat.clear();
    fCosThetaPDF.clear();
    fPhiPDF.clear();
    fPDF1.clear();
    fOmega.clear();
    fOmega2.clear();
    std::cout << "GEventsGenerator destructor calling. . ." << std::endl;
  }
  void GEventsGenerator::GSetPhotonPolarization(const double& lambda)
  {
    fPars[7] = lambda;
  }
  void GEventsGenerator::GSetQPMC(const double& gammaQPCM, const double& f2, 
                const double&  thetaK1, const double& phiDK, const double& phiSRho, 
                const double& phiDRho, const std::complex<double>& gg)
  {
    fPars[0] = gammaQPCM;
    fPars[1] = f2;
    fPars[2] = thetaK1;
    fPars[3] = phiDK;
    fPars[4] = phiSRho;
    fPars[5] = phiDRho;
    fPars[6] = gg.real();
    fPars[10] = gg.imag();
  }
  void GEventsGenerator::GSetK1_1270(const int& charge)
  {
    fResNames.push_back("K1_1270");
    fCharge.push_back(charge);
    fPars[29] = 1;
  }
  void GEventsGenerator::GSetK1_1400(const int& charge, const std::complex<double>& ff)
  {
    fResNames.push_back("K1_1400");
    fCharge.push_back(charge);
    fPars[8] = ff.real();
    fPars[9] = ff.imag();
  }
  void GEventsGenerator::GSetKStr_1410(const int& charge, const std::complex<double>& ff,
                       const std::complex<double>& gg3Kstr,const std::complex<double>& gg3Rho)
  {
    fResNames.push_back("K*_1410");
    fCharge.push_back(charge);
    fPars[11] = ff.real();
    fPars[12] = ff.imag();
    fPars[17] = gg3Kstr.real();
    fPars[18] = gg3Kstr.imag();
    fPars[19] = gg3Rho.real();
    fPars[20] = gg3Rho.imag();
  }
  void GEventsGenerator::GSetKStr_1680(const int& charge, const std::complex<double>& ff,
                       const std::complex<double>& gg3Kstr,const std::complex<double>& gg3Rho)
  {
    fResNames.push_back("K*_1680");
    fCharge.push_back(charge);
    fPars[13] = ff.real();
    fPars[14] = ff.imag();
    fPars[21] = gg3Kstr.real();
    fPars[22] = gg3Kstr.imag();
    fPars[23] = gg3Rho.real();
    fPars[24] = gg3Rho.imag();
  }
  void GEventsGenerator::GSetK2_1430(const int& charge, const std::complex<double>& ff,
                       const std::complex<double>& gg3Kstr,const std::complex<double>& gg3Rho)
  {
    fResNames.push_back("K2_1430");
    fCharge.push_back(charge);
    fPars[15] = ff.real();
    fPars[16] = ff.imag();
    fPars[25] = gg3Kstr.real();
    fPars[26] = gg3Kstr.imag();
    fPars[27] = gg3Rho.real();
    fPars[28] = gg3Rho.imag();
  }
  void GEventsGenerator::GSetEventsNumber(const int& nEvents)
  {
    fNEvents = nEvents;
  }
  void GEventsGenerator::GSetDecayMode(const std::vector<std::string>& decMode, const std::vector<int>& charges)
  {
    (*fMf)->GSetResonances(decMode, charges);
  }
  
  void GEventsGenerator::GGenerateEvents()
  {
    GSetFunctionParameters(30, fPars);
    (*fMf)->GSetResonances(fResNames, fCharge);
    srand (time(NULL));
    int nEvents = 0, nAllEvents = 0;
    double ui;
    if(!fNameOfMinimizer)
    {
      std::cout << std::endl;
      std::cout << "Error!!! Please specify minimizer (\"minuit\" or \"mc\") in main.C file!" << std::endl;
      std::cout << std::endl;
      exit(1);
    }
    else if(std::strcmp(fNameOfMinimizer, "minuit") == 0)
      GMinimizeWithMinuit();
    else if(std::strcmp(fNameOfMinimizer, "mc") == 0)
      GMinimizeWithMonteCarlo();
    else
    {
      std::cout << std::endl;
      std::cout << "Error!!!Minimizer is unknown =( Please specify minimizer (\"minuit\" or \"mc\") in main.C file!" << std::endl;
      std::cout << std::endl;
      exit(1);
    }
    std::cout << -fMinFunctionValue << std::endl;
    double inVars[7];
    double maxPDFValue = -fMinFunctionValue;
    double sij,sjk, sik; // Kpi1, pi1pi2, Kpi2 
    double lowsjk, upsjk;
    
    double mi = (*fMf)->GGetMasses()[0];
    double mj = (*fMf)->GGetMasses()[1];
    double mk = (*fMf)->GGetMasses()[2];
    double M = (*fMf)->GGetMasses()[3];
    double s;
    double pdf;   
    clock_t tStart = clock();
    while(nEvents < fNEvents)
    {
      s = fLowVarsLimits[4] + (fUpVarsLimits[4]-fLowVarsLimits[4])*double(rand())/RAND_MAX; //
      M = sqrt(s);
      fUpVarsLimits[2] = ( M - mk ) * ( M - mk );
      fUpVarsLimits[3] = ( M - mi ) * ( M - mi ); //
      sjk = fLowVarsLimits[3] + (fUpVarsLimits[3] - fLowVarsLimits[3])*double(rand())/RAND_MAX; //
      sij = fLowVarsLimits[2] + (fUpVarsLimits[2] - fLowVarsLimits[2])*double(rand())/RAND_MAX; //
      lowsjk = (*fMf)->GSjkMin(sij, mi, mj, mk, M);
      upsjk = (*fMf)->GSjkMax(sij, mi, mj, mk, M); 
      if((sjk < lowsjk) || (sjk > upsjk))
        continue;
      sik = mi*mi+mj*mj+mk*mk+M*M-sjk-sij;
      inVars[0] = fLowVarsLimits[0] + (fUpVarsLimits[0]-fLowVarsLimits[0])*double(rand())/RAND_MAX; //
      inVars[1] = fLowVarsLimits[1] + (fUpVarsLimits[1]-fLowVarsLimits[1])*double(rand())/RAND_MAX; //
      inVars[2] = sij;
      inVars[3] = sjk; 
      inVars[4] = M*M;
      ui = double(rand())/RAND_MAX;
      nAllEvents++;
      pdf = -fFuncPointer->EvalPar(inVars,fParsValues);
      if(nAllEvents < fNEvents)
      {
        fCosThetaPDFFlat.push_back(inVars[0]);
        fPhiPDFFlat.push_back(inVars[1]);
        fPDF1Flat.push_back(pdf);
        fOmegaFlat.push_back((*fMf)->GOmegaPDF1());
        fOmega2Flat.push_back((*fMf)->GOmegaPDF1()*(*fMf)->GOmegaPDF1());
        fSijFlat.push_back(sij);
        fSjkFlat.push_back(sjk);
        fSikFlat.push_back(sik);
        fMFlat.push_back(M);
      }
      if(ui*maxPDFValue < pdf)
      {
        nEvents++;
        if(nEvents%10000 == 0 )
                  std::cout << nEvents << std::endl;
        fCosThetaPDF.push_back(inVars[0]);
        fPhiPDF.push_back(inVars[1]);
        fPhiPDF2.push_back(inVars[6]);
        fPDF1.push_back(pdf);
//         std::cout << pdf << std::endl;
        fOmega.push_back((*fMf)->GOmegaPDF1());
        fOmega2.push_back((*fMf)->GOmegaPDF1()*(*fMf)->GOmegaPDF1());
        fSij.push_back(sij);
        fSjk.push_back(sjk);
        fSik.push_back(sik);
        fM.push_back(M);
      }
    }
    std::cout << nAllEvents << std::endl;
    if(!fReadWrite)
    {
      std::ofstream ost;
      ost.open(fFileWithIntegrals);
      for(int iInts = 0; iInts < fNIntegrals; ++iInts)
        ost << std::real(fNormalizationIntegrals[iInts]) << "  " << std::imag(fNormalizationIntegrals[iInts]) << std::endl;
      ost.close();
    }
    (*fMf)->GSetNormalizationIntegrals(fNormalizationIntegrals);
    
    std::cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s " << std::endl;
  }
}