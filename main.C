#include "GEventsGenerator.h"
#include "GMathFunctions.h"
#include "GWriteToTree.h"
#include "GInterfaceForMinimization.h"
#include "GPhysicsPDG.h"
#include "TGraph2D.h"
#include "TTree.h"
#include <ctime>
#include <map>
#include <iterator>

using namespace Gamapola;

int main()
{
  clock_t tStart = clock();
  srand (time(NULL));
  double mK = 0.498, mPi1 = 0.140, mPi2 = 0.140, Ecm= 1.272;
  // set lower and upper limits on kinematic variables
  double s12Min = (mK + mPi1) * (mK + mPi1);
  double s23Min = (mPi1 + mPi2) * (mPi1 + mPi2);
  double s13Min = (mK + mPi2) * (mK + mPi2);
  double s12Max = (Ecm- mPi2) * (Ecm- mPi2);
  double s23Max = (Ecm- mK) * (Ecm- mK);
  double s13Max = (Ecm- mPi1) * (Ecm- mPi1);
  double sMin = 1.;
  double sMax = (Ecm+0.6)*(Ecm+0.6);
  // there could be 3 types of variables: "limited", "fixed", "unspecified" = unlimited
  const char* kindOfVars0[5] = {"limited", "limited", "limited","limited", "limited"};
  const char* namesOfVars0[5] = {"theta","phi1","s12", "s23", "s"};
  double initialValues[5] = {0,0, s12Min , s23Min, sMin};
  double lowBoundValues[5] = {-1,0, s12Min, s23Min, sMin};
  double upBoundValues[5] = {1,6.28,s12Max, s23Max, sMax};
  double stepValues[5] = {0.01,0.01, 0.05, 0.05, 0.05};
  std::shared_ptr<GEventsGenerator> eg(new GEventsGenerator());
  std::shared_ptr<GMathFunctions> mf(new GMathFunctions());
  std::shared_ptr<GWriteToTree> wt(new GWriteToTree());
  
  eg->GSetPhotonPolarization(0.8);
  
  // set hadronic parameters
  eg->GSetQPMC(4., // gammaQPCM
               0.5, // f2
               TMath::Pi()/3., // thetaK1
               TMath::Pi()/2, // phiDK*
               TMath::Pi()/2., // phiSRho
               TMath::Pi()/2, // phiDRho
               std::complex<double>(1.0,0.0)); // kappa-contribution
  // set K1_1270: charge = 1
  eg->GSetK1_1270(1);
  
  // set K1_1400
  eg->GSetK1_1400(1, // charge
                  std::complex<double>(0.47,0.0)); // coupling
  
  // set K*_1410
  eg->GSetKStr_1410(1, // charge
                    std::complex<double>(0.784,0.0), // coupling
                       std::complex<double>(1.0,0.0), // K*-contribution
                       std::complex<double>(1.0,0.0)); // rho-contribution
  
  // set K*_1680:
  eg->GSetKStr_1680(1, // charge
                    std::complex<double>(1.23,0.0), // coupling 
                       std::complex<double>(1.0,0.0), // K*-contribution
                       std::complex<double>(1.0,0.0)); // rho-contribution
 
  // set K2_1430:
  eg->GSetK2_1430(1, // charge
                  std::complex<double>(0.444,0.0), // coupling 
                       std::complex<double>(1.0,0.0), // K*-contribution
                       std::complex<double>(1.0,0.0)); // rho-contribution

  int nEvents = 50000;
  mf->GSetMasses(mK, mPi1, mPi2, Ecm);
  eg->GSetMinimizer(); // "mc" by default
  eg->GSetFunctionVariables(5, initialValues, lowBoundValues, upBoundValues, stepValues, kindOfVars0, namesOfVars0);
  eg->GSetEventsNumber(nEvents);
  eg->GSetFunction(&Gamapola::GMathFunctions::GProcessingComputationOfPDF, &mf, "GProcessingComputationOfPDF", "GMathFunctions");
  eg->GGenerateEvents();
  
  wt->GCreateFile("Events.root");
  wt->GMakeDirectory("non-flatEvents");
  wt->GChangeDirectory("non-flatEvents");
 
  wt->GMakeDirectory("flatEvents");
  wt->GChangeDirectory("flatEvents");
  wt->GAddHistogram("costhetaFlat", eg->GGetCosThetaPDFFlat(), -1.,1., 100);
  wt->GAddHistogram("phiFlat", eg->GGetPhiPDFFlat(), 0,6.28, 500);
  wt->GAddHistogram("omegaFlat", eg->GGetOmegaFlat(), -1.,1., 500);
  wt->GAddHistogram("omega2Flat", eg->GGetOmega2Flat(), 0,1, 500);
  wt->GAddHistogram("PDF1Flat", eg->GGetPDF1Flat(), 0., 50, 500);
  wt->GAddHistogram2D("s12s23", eg->GGetSijFlat(), 0,(Ecm+0.5)*(Ecm+0.5), 100, eg->GGetSjkFlat(), 0,(Ecm+0.5)*(Ecm+0.5), 100);
  wt->GAddHistogram2D("s12s13", eg->GGetSijFlat(), 0,(Ecm+0.5)*(Ecm+0.5), 100, eg->GGetSikFlat(), 0,(Ecm+0.5)*(Ecm+0.5), 100);
  wt->GAddHistogram2D("s23s13", eg->GGetSjkFlat(), 0,(Ecm+0.5)*(Ecm+0.5), 100, eg->GGetSikFlat(), 0,(Ecm+0.5)*(Ecm+0.5), 100);
  wt->GAddHistogram("s12", eg->GGetSijFlat(), s12Min, s12Max, 500);
  wt->GAddHistogram("s23", eg->GGetSjkFlat(), s23Min, s23Max, 500);
  wt->GAddHistogram("s13", eg->GGetSikFlat(), s13Min, s13Max, 500);
  wt->GAddHistogram("MFlat", eg->GGetMFlat(), (mK+mPi1+mPi2), (Ecm+1), 500);
  wt->GWriteHistograms();
  wt->GChangeDirectory("..");
  std::cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s " << std::endl;
  
  return 0;
}
