#include "GBreitWigner.h"

namespace Gamapola{
  GBreitWigner::GBreitWigner():GInterfaceForMathFunctions()
  {
    std::cout << "GBreitWigner constructor calling. . ." << std::endl;
  }
  GBreitWigner::~GBreitWigner()
  {
    std::cout << "GBreitWigner destructor calling. . ." << std::endl;
  }
  std::complex<double> GBreitWigner::GBWKstr(double sij)
  {
    std::complex<double> denom(sij - kMK0star_892 * kMK0star_892, kGammaK0star_892 * kMK0star_892);
    return 1./denom;
  }
  std::complex<double> GBreitWigner::GBWRho(double sij)
  {
    double rel = GMomentaA(sij, kMPion * kMPion, kMPion) / GMomentaA(kMRho0_775 * kMRho0_775,  kMPion * kMPion, kMPion);
    std::complex<double> denom(sij - kMRho0_775 * kMRho0_775, kGammaRho0_775 * kMRho0_775 * rel * rel * rel);
    return  1./denom;
  }
  std::complex<double> GBreitWigner::GBWK1(double s, int nRes)
  {
//     std::cout << fResMasses[nRes] << "  " << fResWidths[nRes] << "   " << nRes << std::endl;
    std::complex<double> denom(s - fResMasses[nRes] * fResMasses[nRes], fResWidths[nRes] * fResMasses[nRes]);
    return 1./denom;
  }
  std::complex<double> GBreitWigner::GBWKappa(double sij)
  {
    std::complex<double> denom(sij - kMKappa_800 * kMKappa_800, kGammaKappa_800 * kMKappa_800);
    return 1./denom;
  }
}
