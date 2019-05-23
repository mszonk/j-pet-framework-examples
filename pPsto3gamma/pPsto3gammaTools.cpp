/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file pPsto3gammaTools.cpp
 */

#include "pPsto3gammaTools.h"
#include <TMath.h>
#include <vector>

using namespace std;
double pPsto3gammaTools::calculateAnnihilationTime(const JPetHit& firstHit, const TVector3& Vtx)
{
  TVector3 hit1=calculateHitPositionVector(firstHit,Vtx);
  return (firstHit.getTime()/1000. -(hit1.Mag())/kLightVelocity_cm_ns);
}

TVector3 pPsto3gammaTools::calculateHitPositionVector(const JPetHit& hit1, const TVector3& Vtx)
{
  TVector3 Rvec(hit1.getPosX()-Vtx.X(),hit1.getPosY()-Vtx.Y(),hit1.getPosZ()-Vtx.Z());
  return Rvec;
}

double pPsto3gammaTools::calculateDotProduct(const TVector3& hit1, const TVector3& hit2)
{
  double DotProduct = hit1.X()*hit2.X() + hit1.Y()*hit2.Y() + hit1.Z()*hit2.Z();
  return DotProduct;
}

/**
* Method for determining type of event - back to back 2 gamma
*/
int pPsto3gammaTools::checkFor2Gamma(const JPetEvent& event, JPetStatistics& stats,
				     bool saveHistos, double b2bThetaDiff,double AnnihTimeDiffcut2g,
				     double RxyCut, double VtxzCut, std::vector<int>& TwogammaCandidatesId,
				     std::vector<double>& AnnihTimes)
{
  int Npairs=0;
  if (event.getHits().size() < 2) {
    return 0;
  }
 
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      JPetHit firstHit, secondHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        firstHit = event.getHits().at(i);
        secondHit = event.getHits().at(j);
      } else {
        firstHit = event.getHits().at(j);
        secondHit = event.getHits().at(i);
      }
        TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
	double AnnihTime1 = calculateAnnihilationTime(firstHit,annhilationPoint);
	double AnnihTime2 = calculateAnnihilationTime(secondHit,annhilationPoint);
	//
	//time diference between emmisions of the two photons
	double AnnihTimeDiff= AnnihTime1 - AnnihTime2;
	if(fabs(AnnihTimeDiff)<=AnnihTimeDiffcut2g){
	  double VtxRxy = sqrt(annhilationPoint.X()*annhilationPoint.X() + annhilationPoint.Y()*annhilationPoint.Y());
	  double Vtxz = annhilationPoint.Z();
	  if(VtxRxy<=RxyCut && Vtxz<=VtxzCut){ 
	    // Checking for back to back
	    TVector3 hit1=calculateHitPositionVector(firstHit,annhilationPoint);
	    TVector3 hit2=calculateHitPositionVector(secondHit,annhilationPoint);
	    double DProduct = calculateDotProduct(hit1,hit2);
	    double thetaDiff = DProduct/(hit1.Mag()*hit2.Mag());
	    thetaDiff =acos(thetaDiff);
	    double minTheta = 180.0 - b2bThetaDiff;
	    double maxTheta = 180.0 + b2bThetaDiff;
	    if (thetaDiff >= minTheta && thetaDiff <= maxTheta) { 
	      Npairs++;
	      AnnihTimes.push_back(AnnihTime1);
	      AnnihTimes.push_back(AnnihTime2);
	      TwogammaCandidatesId.push_back(i);
	      TwogammaCandidatesId.push_back(j);
	      if (saveHistos) {
		double distance = calculateDistance(secondHit, firstHit);
	   // TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
		stats.getHisto1D("2Gamma_Zpos")->Fill(firstHit.getPosZ());
		stats.getHisto1D("2Gamma_Zpos")->Fill(secondHit.getPosZ());
		stats.getHisto1D("2Gamma_TimeDiff")->Fill(AnnihTimeDiff);
		stats.getHisto1D("2Gamma_Dist")->Fill(distance);
		stats.getHisto1D("2Gamma_AnnihTimeDiff")->Fill(AnnihTimeDiff);
		stats.getHisto1D("Annih_TOF")->Fill(calculateTOF(firstHit, secondHit));
		stats.getHisto2D("AnnihPoint_XY")->Fill(annhilationPoint.X(), annhilationPoint.Y());
		stats.getHisto2D("AnnihPoint_XZ")->Fill(annhilationPoint.X(), annhilationPoint.Z());
		stats.getHisto2D("AnnihPoint_YZ")->Fill(annhilationPoint.Y(), annhilationPoint.Z());
	    }
	   }
	  }
	}
    }
  }
  return Npairs;
}

/**
* Method for determining type of event - 3Gamma
*/
bool pPsto3gammaTools::checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos)
{
  if (event.getHits().size() < 3) return false;
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      for (uint k = j + 1; k < event.getHits().size(); k++) {
        JPetHit firstHit = event.getHits().at(i);
        JPetHit secondHit = event.getHits().at(j);
        JPetHit thirdHit = event.getHits().at(k);

        vector<double> thetaAngles;
        thetaAngles.push_back(firstHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(secondHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(thirdHit.getBarrelSlot().getTheta());
        sort(thetaAngles.begin(), thetaAngles.end());

        vector<double> relativeAngles;
        relativeAngles.push_back(thetaAngles.at(1) - thetaAngles.at(0));
        relativeAngles.push_back(thetaAngles.at(2) - thetaAngles.at(1));
        relativeAngles.push_back(360.0 - thetaAngles.at(2) + thetaAngles.at(0));
        sort(relativeAngles.begin(), relativeAngles.end());
        double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
        double transformedY = relativeAngles.at(1) - relativeAngles.at(0);

        if (saveHistos) {
          stats.getHisto2D("3Gamma_Angles")->Fill(transformedX, transformedY);
        }
      }
    }
  }
  return true;
}

/**
* Method for determining type of event - prompt
*/
int  pPsto3gammaTools::checkForPrompt(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double deexTOTCutMin, double deexTOTCutMax, std::vector<int>& PromptCandidate)
{
  int nPrompts=0;
  for (unsigned i = 0; i < event.getHits().size(); i++) {
    double tot = calculateTOT(event.getHits().at(i));
    if (tot > deexTOTCutMin && tot < deexTOTCutMax) {
      PromptCandidate.push_back(nPrompts);
      nPrompts++;
      if (saveHistos) {
        stats.getHisto1D("Deex_TOT_cut")->Fill(tot);
      }
    }
  }
  return nPrompts;
}

/**
* Method for determining type of event - scatter
*/
bool pPsto3gammaTools::checkForScatter(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos, double scatterTOFTimeDiff
)
{
  if (event.getHits().size() < 2) {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      JPetHit primaryHit, scatterHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        primaryHit = event.getHits().at(i);
        scatterHit = event.getHits().at(j);
      } else {
        primaryHit = event.getHits().at(j);
        scatterHit = event.getHits().at(i);
      }

      double scattAngle = calculateScatteringAngle(primaryHit, scatterHit);
      double scattTOF = calculateScatteringTime(primaryHit, scatterHit);
      double timeDiff = scatterHit.getTime() - primaryHit.getTime();

      if (saveHistos) {
        stats.getHisto1D("ScatterTOF_TimeDiff")->Fill(fabs(scattTOF - timeDiff));
      }

      if (fabs(scattTOF - timeDiff) < scatterTOFTimeDiff) {
        if (saveHistos) {
          stats.getHisto2D("ScatterAngle_PrimaryTOT")->Fill(scattAngle, calculateTOT(primaryHit));
          stats.getHisto2D("ScatterAngle_ScatterTOT")->Fill(scattAngle, calculateTOT(scatterHit));
        }
        return true;
      }
    }
  }
  return false;
}

/**
* Calculation of the total TOT of the hit - Time over Threshold:
* the sum of the TOTs on all of the thresholds (1-4) and on the both sides (A,B)
*/
double pPsto3gammaTools::calculateTOT(const JPetHit& hit)
{
  double tot = 0.0;

  auto sigALead = hit.getSignalA().getRecoSignal().getRawSignal()
                  .getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrNum);
  auto sigBLead = hit.getSignalB().getRecoSignal().getRawSignal()
                  .getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrNum);
  auto sigATrail = hit.getSignalA().getRecoSignal().getRawSignal()
                   .getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrNum);
  auto sigBTrail = hit.getSignalB().getRecoSignal().getRawSignal()
                   .getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrNum);

  if (sigALead.size() > 0 && sigATrail.size() > 0) {
    for (unsigned i = 0; i < sigALead.size() && i < sigATrail.size(); i++) {
      tot += (sigATrail.at(i).getValue() - sigALead.at(i).getValue());
    }
  }
  if (sigBLead.size() > 0 && sigBTrail.size() > 0) {
    for (unsigned i = 0; i < sigBLead.size() && i < sigBTrail.size(); i++) {
      tot += (sigBTrail.at(i).getValue() - sigBLead.at(i).getValue());
    }
  }
  return tot;
}

/**
* Calculation of distance between two hits
*/
double pPsto3gammaTools::calculateDistance(const JPetHit& hit1, const JPetHit& hit2)
{
  return (hit1.getPos() - hit2.getPos()).Mag();
}

/**
* Calculation of time that light needs to travel the distance between primary gamma
* and scattered gamma. Return value in picoseconds.
*/
double pPsto3gammaTools::calculateScatteringTime(const JPetHit& hit1, const JPetHit& hit2)
{
  return 1000. * calculateDistance(hit1, hit2) / kLightVelocity_cm_ns;
}

/**
* Calculation of scatter angle between primary hit and scattered hit.
* This function assumes that source of first gamma was in (0,0,0).
* Angle is calculated from scalar product, return value in degrees.
*/
double pPsto3gammaTools::calculateScatteringAngle(const JPetHit& hit1, const JPetHit& hit2)
{
  return TMath::RadToDeg() * hit1.getPos().Angle(hit2.getPos() - hit1.getPos());
}

/**
* Calculation point in 3D, where annihilation occured
*/
TVector3 pPsto3gammaTools::calculateAnnihilationPoint(const JPetHit& hitA, const JPetHit& hitB)
{
  double tof = pPsto3gammaTools::calculateTOF(hitA, hitB);
  return calculateAnnihilationPoint(hitA.getPos(), hitB.getPos(), tof);
}

TVector3 pPsto3gammaTools::calculateAnnihilationPoint(const TVector3& hitA, const TVector3& hitB, double tof)
{
  TVector3 middleOfLOR = 0.5 * (hitA + hitB);
  TVector3 versorOnLOR = (hitB - hitA).Unit()  ;

  double shift = 0.5 * tof  * kLightVelocity_cm_ns / 1000.0;
  TVector3 annihilationPoint(middleOfLOR.X() + shift * versorOnLOR.X(),
                             middleOfLOR.Y() + shift * versorOnLOR.Y(),
                             middleOfLOR.Z() + shift * versorOnLOR.Z());
  return annihilationPoint;
}

double pPsto3gammaTools::calculateTOFByConvention(const JPetHit& hitA, const JPetHit& hitB)
{
  if (hitA.getBarrelSlot().getTheta() < hitB.getBarrelSlot().getTheta()) {
    return calculateTOF(hitA, hitB);
  } else {
    return calculateTOF(hitB, hitA);
  }
}

double pPsto3gammaTools::calculateTOF(const JPetHit& hitA, const JPetHit& hitB)
{
  return pPsto3gammaTools::calculateTOF(hitA.getTime(), hitB.getTime());
}

double pPsto3gammaTools::calculateTOF(double time1, double time2)
{
  return (time1 - time2);
}

/**
* Calculating distance from the center of the decay plane
*/
double pPsto3gammaTools::calculatePlaneCenterDistance(
  const JPetHit& firstHit, const JPetHit& secondHit, const JPetHit& thirdHit)
{
  TVector3 crossProd = (secondHit.getPos() - firstHit.getPos()).Cross(thirdHit.getPos() - secondHit.getPos());
  double distCoef = -crossProd.X() * secondHit.getPosX() - crossProd.Y() * secondHit.getPosY() - crossProd.Z() * secondHit.getPosZ();
  if (crossProd.Mag() != 0) {
    return fabs(distCoef) / crossProd.Mag();
  } else {
    ERROR("One of the hit has zero position vector - unable to calculate distance from the center of the surface");
    return -1.;
  }
}

/**
* Method for determining type of event for streaming - 2 gamma
* @todo: the selection criteria b2b distance from center needs to be checked
* and implemented again
*/
bool pPsto3gammaTools::stream2Gamma(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double b2bSlotThetaDiff, double b2bTimeDiff
)
{
  if (event.getHits().size() < 2) {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      JPetHit firstHit, secondHit;
      if (event.getHits().at(i).getTime() < event.getHits().at(j).getTime()) {
        firstHit = event.getHits().at(i);
        secondHit = event.getHits().at(j);
      } else {
        firstHit = event.getHits().at(j);
        secondHit = event.getHits().at(i);
      }
      // Checking for back to back
      double timeDiff = fabs(firstHit.getTime() - secondHit.getTime());
      double deltaLor = (secondHit.getTime() - firstHit.getTime()) * kLightVelocity_cm_ns / 2000.;
      double theta1 = min(firstHit.getBarrelSlot().getTheta(), secondHit.getBarrelSlot().getTheta());
      double theta2 = max(firstHit.getBarrelSlot().getTheta(), secondHit.getBarrelSlot().getTheta());
      double thetaDiff = min(theta2 - theta1, 360.0 - theta2 + theta1);
      if (saveHistos) {
        stats.getHisto1D("2Gamma_TimeDiff")->Fill(timeDiff / 1000.0);
        stats.getHisto1D("2Gamma_DLOR")->Fill(deltaLor);
        stats.getHisto1D("2Gamma_ThetaDiff")->Fill(thetaDiff);
      }
      if (fabs(thetaDiff - 180.0) < b2bSlotThetaDiff && timeDiff < b2bTimeDiff) {
        if (saveHistos) {
          TVector3 annhilationPoint = calculateAnnihilationPoint(firstHit, secondHit);
          stats.getHisto1D("2Annih_TimeDiff")->Fill(timeDiff / 1000.0);
          stats.getHisto1D("2Annih_DLOR")->Fill(deltaLor);
          stats.getHisto1D("2Annih_ThetaDiff")->Fill(thetaDiff);
          stats.getHisto2D("2Annih_XY")->Fill(annhilationPoint.X(), annhilationPoint.Y());
          stats.getHisto1D("2Annih_Z")->Fill(annhilationPoint.Z());
        }
        return true;
      }
    }
  }
  return false;
}

/**
* Method for determining type of event for streaming - 3 gamma annihilation
*/
bool pPsto3gammaTools::stream3Gamma(
  const JPetEvent& event, JPetStatistics& stats, bool saveHistos,
  double d3SlotThetaMin, double d3TimeDiff, double d3PlaneCenterDist
)
{
  if (event.getHits().size() < 3) {
    return false;
  }
  for (uint i = 0; i < event.getHits().size(); i++) {
    for (uint j = i + 1; j < event.getHits().size(); j++) {
      for (uint k = j + 1; k < event.getHits().size(); k++) {
        JPetHit firstHit = event.getHits().at(i);
        JPetHit secondHit = event.getHits().at(j);
        JPetHit thirdHit = event.getHits().at(k);

        vector<double> thetaAngles;
        thetaAngles.push_back(firstHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(secondHit.getBarrelSlot().getTheta());
        thetaAngles.push_back(thirdHit.getBarrelSlot().getTheta());
        sort(thetaAngles.begin(), thetaAngles.end());

        vector<double> relativeAngles;
        relativeAngles.push_back(thetaAngles.at(1) - thetaAngles.at(0));
        relativeAngles.push_back(thetaAngles.at(2) - thetaAngles.at(1));
        relativeAngles.push_back(360.0 - thetaAngles.at(2) + thetaAngles.at(0));
        sort(relativeAngles.begin(), relativeAngles.end());

        double transformedX = relativeAngles.at(1) + relativeAngles.at(0);
        double transformedY = relativeAngles.at(1) - relativeAngles.at(0);
        double timeDiff = fabs(thirdHit.getTime() - firstHit.getTime());
        double planeCenterDist = calculatePlaneCenterDistance(firstHit, secondHit, thirdHit);
        if (saveHistos) {
          stats.getHisto1D("3GammaTimeDiff")->Fill(timeDiff);
          stats.getHisto2D("3GammaThetas")->Fill(transformedX, transformedY);
          stats.getHisto1D("3GammaPlaneDist")->Fill(planeCenterDist);
        }
        if (transformedX > d3SlotThetaMin && timeDiff < d3TimeDiff && planeCenterDist < d3PlaneCenterDist) {
          if (saveHistos) {
            stats.getHisto1D("3AnnihPlaneDist")->Fill(planeCenterDist);
            stats.getHisto1D("3AnnihTimeDiff")->Fill(timeDiff);
          }
          return true;
        }
      }
    }
  }
  return false;
}
