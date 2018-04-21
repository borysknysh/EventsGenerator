#include "GKinematics.h"

namespace Gamapola{
  GKinematics::GKinematics():GInterfaceForMathFunctions()
  {
    std::cout << "GKinematics constructor calling. . ." << std::endl;
  }
  GKinematics::~GKinematics()
  {
    std::cout << "GKinematics destructor calling. . ." << std::endl;
  }
  void GKinematics::GKinematicsConstanstsCalculation(int nRes)
  {
    for(size_t i = 0; i < fMomentaConst[fIndex[nRes]].size(); i++)
      { 
//         std::cout << i << "   " << nRes << "  " << fIndex[nRes] << "  " << std::endl;
        fMomentaConst[fIndex[nRes]][i] = GMomentaA(fConstParsForMomA[fIndex[nRes]][i][2] * fConstParsForMomA[fIndex[nRes]][i][2],
        fConstParsForMomA[fIndex[nRes]][i][1]*fConstParsForMomA[fIndex[nRes]][i][1],fConstParsForMomA[fIndex[nRes]][i][0]);
        fMomentaConst[fIndex[nRes]][i] = (isnan(fMomentaConst[fIndex[nRes]][i]))?0:fMomentaConst[fIndex[nRes]][i];
      }
//       std::cout << fMomentaConst[fIndex[0]][1] << std::endl;
  }
  double GKinematics::GSjkMin(const double& sij, const double& mi, const double& mj, const double& mk, const double& M)
  {
    auto Ejstar = ( sij - mi * mi + mj * mj ) / ( 2 * sqrt ( sij ) );
    auto Ekstar = ( M * M - sij - mk * mk ) / ( 2 * sqrt ( sij ) );
    auto sjk = ( Ejstar + Ekstar ) * ( Ejstar + Ekstar ) - 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) + sqrt( Ekstar * Ekstar - mk * mk ) ) * 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) + sqrt( Ekstar * Ekstar - mk * mk ) );
    return  sjk ;
  }
  double GKinematics::GSjkMax(const double& sij, const double& mi, const double& mj, const double& mk, const double& M)
  {
    auto Ejstar = ( sij - mi * mi + mj * mj ) / ( 2 * sqrt ( sij ) );
    auto Ekstar = ( M * M - sij - mk * mk ) / ( 2 * sqrt ( sij ) );
    auto sjk = ( Ejstar + Ekstar ) * ( Ejstar + Ekstar ) - 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) - sqrt( Ekstar * Ekstar - mk * mk ) ) * 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) - sqrt( Ekstar * Ekstar - mk * mk ) );
    return  sjk ;
  } 
  auto GKinematics::GMomentaA(double s, double sij, double mb)->decltype(s)
  { 
    return sqrt ( ( s - ( sqrt( sij ) + mb) * (sqrt( sij ) + mb) ) * 
    (s - (sqrt( sij ) - mb) * (sqrt( sij ) - mb) ) / ( 4. * s ) );
  }
  auto GKinematics::GEnergyA(double s, double sij, double mb)->decltype(s)
  {
//     std::cout << sqrt(s) << "  " << mb << "  " << sqrt(sij) << std::endl;
    return ( s - mb * mb + sij ) / ( 2 * sqrt( s ) );
  } 
  auto GKinematics::GPiPjVec4(double s, double sij, double sjk)->decltype(s)
  {
//     std::cout << s << " " << sij << "  " << sjk << "  " <<  kMKaon << std::endl;
    return ( s - sij - sjk + kMKaon * kMKaon ) / 2. ;
  }
  void GKinematics::GPiPjVec3()
  {
    fPiPjVec[0] = fMomiV[1] * fMomjV[2] - fMomiV[2] * fMomjV[1];
    fPiPjVec[1] = fMomiV[2] * fMomjV[0] - fMomiV[0] * fMomjV[2];
    fPiPjVec[2] = fMomiV[0] * fMomjV[1] - fMomiV[1] * fMomjV[0];
//     std::cout << fPiPjVec[0] << "  " << fMomiV[1] * fMomjV[2] << "  " << fMomiV[2] * fMomjV[1] << std::endl;
  }
  auto GKinematics::GSij(double s, double sij, double sjk)->decltype(s)
  {
    return ( s - sij - sjk + kMKaon * kMKaon + kMPion * kMPion + kMPion * kMPion) ;
  }
  std::complex<double> GKinematics::GScPrd(const double* pi, const std::complex<double>* pj)
  {
    return pi[0] * pj[0] + pi[1] * pj[1] + pi[2] * pj[2]; 
  }
  double GKinematics::GPiPj()
  {
//     std::cout << fEnergy1 << "  " << fEnergy2 << "  " << (kMPion * kMPion + kMKaon * kMKaon - fSij) / 2. << std::endl;
    return fEnergy1 * fEnergy2 + (kMPion * kMPion + kMPion * kMPion - fSij) / 2.;
  }
  double GKinematics::GExponentialM(double* p, int nRes, int nFnc)
  {
    auto f2 = p[1];
//     std::cout << TMath::Exp( f2 * fMomentaConst[nRes][nFnc] * fMomentaConst[nRes][nFnc] ) << "  _____" <<std::endl;
    return TMath::Exp( f2 * fMomentaConst[nRes][nFnc] * fMomentaConst[nRes][nFnc] );
  }
  double GKinematics::GCosThetaij(const double& cosThstar, const double& sij, 
                                  const double& mi, const double& mj, 
                                  const double& mk, const double& M)
  {
    auto Eistar = ( sij + mi * mi - mj * mj ) / ( 2 * sqrt ( sij ) );
    auto Ejstar = ( sij + mj * mj - mi * mi ) / ( 2 * sqrt ( sij ) );
    
    auto momistar = sqrt( Eistar * Eistar - mi * mi );
    
    auto momiXstar = momistar * sqrt( 1 - cosThstar * cosThstar );
    auto momiZstar = momistar * cosThstar;
    
    auto momjXstar = -momiXstar;
    auto momjZstar = -momiZstar;
    
    auto Eij = M - ( M * M + mk * mk - sij ) / ( 2 * M );
    auto momij = sqrt( ( M * M - ( sqrt ( sij ) + mk ) * ( sqrt ( sij ) + mk ) ) *
    ( M * M - ( sqrt ( sij ) - mk ) * ( sqrt ( sij ) - mk ) ) ) / ( 2 * M );
    
    auto gammaij = Eij / sqrt ( sij );
    auto vij = momij / Eij;
    
    auto Ei = ( Eistar + vij * momiZstar ) * gammaij;
    auto Ej = ( Ejstar + vij * momjZstar ) * gammaij;
    
    auto momiZ = gammaij * ( vij * Eistar + momiZstar );
    auto momjZ = gammaij * ( vij * Ejstar + momjZstar );
    
    auto momi = sqrt( momiXstar * momiXstar + momiZ * momiZ );
    auto momj = sqrt( momjXstar * momjXstar + momjZ * momjZ );
    
    auto cosThetaij = ( mi * mi + mj * mj + 2 * Ei * Ej - sij ) / ( 2 * momi * momj );
    
    return fabs(cosThetaij);
  }
  void GKinematics::GCalculateKinematics(double* x)
  {
    auto sinTh = sin(acos(x[0]));
    fEnergy1 = GEnergyA(x[4], kMPion*kMPion, sqrt(fSik));
    fEnergy2 = GEnergyA(x[4], kMPion*kMPion, sqrt(x[2]));
    fPiPjVec4 = GPiPjVec4(x[4], x[2], fSik);
    fMomij = GPiPj();
    fMomi = sqrt(fEnergy1 * fEnergy1 - kMPion * kMPion);
    fMomj = sqrt(fEnergy2 * fEnergy2 - kMPion * kMPion);
        
    auto delta = acos(fMomij / (fMomi * fMomj));
    fPhii = (2 * x[1] - delta) / 2;
    fPhij = (2 * x[1] + delta) / 2;
//     std::cout <<"Delta: " << sin(delta) * sinTh * fMomi * fMomj << std::endl;
//     std::cout << fMomi << "  " << fMomj << "   " << fMomij << std::endl;
    fMomiV[0] = fMomi * x[0] * cos(fPhii);
    fMomiV[1] = fMomi * sin(fPhii);
    fMomiV[2] = -fMomi * sinTh * cos(fPhii);
    fMomjV[0] = fMomj * x[0] * cos(fPhij);
    fMomjV[1] = fMomj * sin(fPhij);
    fMomjV[2] = -fMomj * sinTh * cos(fPhij);
//     std::cout << fMomiV[1] << "  " << fMomiV[2] << std::endl;
//     std::cout << fMomjV[1] << "  " << fMomjV[2] << std::endl;
    GPiPjVec3();
    fE0PiPj = GScPrd(fPiPjVec, fEpsilon0);
    fE0Pi = GScPrd(fMomiV, fEpsilon0);
    fE0Pj = GScPrd(fMomjV, fEpsilon0);
    for(size_t iAxis = 0; iAxis < 3; iAxis++)
    {
      fK2KinC4combo[iAxis] = ( fE0Pi * fPiPjVec[iAxis] + fE0PiPj * fMomiV[iAxis] );
      fK2KinC5combo[iAxis] = ( fE0Pj * fPiPjVec[iAxis] + fE0PiPj * fMomjV[iAxis] );
    }
  }
}