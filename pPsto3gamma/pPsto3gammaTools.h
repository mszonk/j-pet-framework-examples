/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file pPsto3gammaTools.h
 */

#ifndef pPsto3gammaTools_H
#define pPsto3gammaTools_H

#include <JPetStatistics/JPetStatistics.h>
#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>

static const double kLightVelocity_cm_ns = 29.9792458;
static const double kUndefinedValue = 999.0;

/**
 * @brief Tools for Event Categorization
 *
 * Lots of tools in constatnt developement.
*/
class pPsto3gammaTools
{
public:
  static double calculateAnnihilationTime(const JPetHit& firstHit, const TVector3& Vtx);
  static TVector3 calculateHitPositionVector(const JPetHit& hit1, const TVector3& Vtx);
  static double calculateDotProduct(const TVector3& hit1, const TVector3& hit2);
  static int checkFor2Gamma(const JPetEvent& event, JPetStatistics& stats,
			    bool saveHistos, double b2bThetaDiff,double AnnihTimeDiffcut, double RxyCut, double VtxzCut, std::vector<int>& TwogammaCandidatesId,
			    std::vector<double>& AnnihTimes);
  static bool checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos);
  static int checkForPrompt(const JPetEvent& event, JPetStatistics& stats,
			    bool saveHistos, double deexTOTCutMin, double deexTOTCutMax,std::vector<int>& PromptCandidate);
  static bool checkForScatter(const JPetEvent& event, JPetStatistics& stats,
                              bool saveHistos, double scatterTOFTimeDiff);
  static double calculateTOT(const JPetHit& hit);
  static double calculateDistance(const JPetHit& hit1, const JPetHit& hit2);
  static double calculateScatteringTime(const JPetHit& hit1, const JPetHit& hit2);
  static double calculateScatteringAngle(const JPetHit& hit1, const JPetHit& hit2);
  /// Tof is calculated as  time1 -time2.
  static double calculateTOF(const JPetHit& hitA, const JPetHit& hitB);
  static double calculateTOF(double time1, double time2);
  /// Tof calculated with the ordered hits with respect to scintillator number.
  /// The first one will be hit with smaller theta angle.
  /// See also: http://koza.if.uj.edu.pl/petwiki/index.php/Coordinate_system_in_Big_Barrel
  static double calculateTOFByConvention(const JPetHit& hitA, const JPetHit& hitB);
  static TVector3 calculateAnnihilationPoint(const JPetHit& hitA, const JPetHit& hitB);
  static TVector3 calculateAnnihilationPoint(const TVector3& hitA, const TVector3& hitB, double tof);
  static double calculatePlaneCenterDistance(const JPetHit& firstHit,
      const JPetHit& secondHit, const JPetHit& thirdHit);
  static bool stream2Gamma(const JPetEvent& event, JPetStatistics& stats,
                           bool saveHistos, double b2bSlotThetaDiff, double b2bTimeDiff);
  static bool stream3Gamma(const JPetEvent& event, JPetStatistics& stats,
                           bool saveHistos, double d3SlotThetaMin, double d3TimeDiff, double d3DistanceFromCenter);
};

#endif /* !pPsto3gammaTools_H */
