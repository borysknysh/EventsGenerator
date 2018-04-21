#ifndef GPHYSICSPDG_H
#define GPHYSICSPDG_H

#include <iostream>
using cd = const double;
namespace Gamapola
{
  cd ka = 10000;
  cd kMBmeson = 5.279;
  cd kMKaon = 0.498;
  cd kMPion = 0.140;
  
  cd kMK1_1270 = 1.272;
  cd kGammaK1_1270 = 0.090;
  cd kBr1K1_1270 = 0.16;
  cd kBr2K1_1270 = 0.42;
  
  cd kMK1_1400 = 1.403;
  cd kGammaK1_1400 = 0.174;
  cd kBr1K1_1400 = 0.94;
  cd kBr2K1_1400 = 0.03;
  
  cd kMKst_1410 = 1.414;
  cd kGammaKst_1410 = 0.232;
  cd kFF1Kst_1410 = 13.34;
  cd kFF2Kst_1410 = 9.2;
  cd kBr1Kst_1410 = 0.93;
  cd kBr2Kst_1410 = 0.07;
  
  cd kMKst_1680 = 1.717;
  cd kGammaKst_1680 = 0.32;
  cd kFF1Kst_1680 = 4.8;
  cd kFF2Kst_1680 = 9.04;
  cd kBr1Kst_1680 = 0.299;
  cd kBr2Kst_1680 = 0.314;
  
  cd kMK2_1600 = 1.605;
  cd kGammaK2_1600 = 0.115;
  cd kMK2_1770 = 1.773;
  cd kGammaK2_1770 = 0.186;
  
  cd kMK2_1430 = 1.426;
  cd kGammaK2_1430 = 0.099;
  cd kFF1K2_1430 = 3.05;
  cd kFF2K2_1430 = 4.47;
  cd kBr1K2_1430 = 0.247;
  cd kBr2K2_1430 = 0.087;
  
//   Kaons resonances
  cd kMK0star_892 = 0.892;
  cd kGammaK0star_892 = 0.05;
  cd kMK0star_1430 = 1.43;
  cd kGammaK0star_1430 = 0.27;
  
//   Kappa resonance
  cd kMKappa_800 = 0.8;
  cd kGammaKappa_800 = 0.5;
  
//   Rho resonances
  cd kMRho0_775 = 0.775;
  cd kGammaRho0_775 = 0.149;
  cd kMRhoPlus = 0.775;
  cd kGammaRhoPlus = 0.150;
  cd kMRhoMinus = 0.775;
  cd kGammaRhoMinus = 0.150;
  
//   
  cd kQPC = 4;
  cd kRpi = 1/0.4;
  cd kRK = 1/0.4;
  cd kRKst = 1/0.4;
  cd kRrho = 1/0.4;
  cd kRomega = 1/0.4;
  cd kRK1 = 1/0.4;
  cd kRb1 = 1/0.4;
  
  const std::complex<double> fEpsilonL[3] = 
  {
    std::complex<double>(1., 0.),  
    std::complex<double>(0., -1.), 
    std::complex<double>(0., 0.)
  };
  const std::complex<double> fEpsilonR[3] = 
  {
    std::complex<double>(-1., 0.),  
    std::complex<double>(0., -1.), 
    std::complex<double>(0., 0.)
    
  };
  const std::complex<double> fEpsilon0[3] = 
  {
    std::complex<double>(0., 0.),  
    std::complex<double>(0., 0.), 
    std::complex<double>(1., 0.)
  };
}
#endif