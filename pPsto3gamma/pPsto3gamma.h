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
 *  @file pPsto3gamma.h
 */

#ifndef pPsto3gamma_H
#define pPsto3gamma_H

#include <JPetUserTask/JPetUserTask.h>
//#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include "pPsto3gammaTools.h" 
#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>
#include <vector>
#include <map>

class JPetWriter;

#ifdef __CINT__
#	define override
#endif

/**
 * @brief Task to search for the pPs->3gamma events
 *
 */
class pPsto3gamma: public JPetUserTask{
public:
	pPsto3gamma(const char * name);
	virtual ~pPsto3gamma();
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;

protected:
	const std::string kBack2BackSlotThetaDiffParamKey = "Back2Back_Categorizer_SlotThetaDiff_float";
	const std::string kScatterTOFTimeDiffParamKey = "Scatter_Categorizer_TOF_TimeDiff_float";
	const std::string kDeexTOTCutMinParamKey = "Deex_Categorizer_TOT_Cut_Min_float";
	const std::string kDeexTOTCutMaxParamKey = "Deex_Categorizer_TOT_Cut_Max_float";
	const std::string kSaveControlHistosParamKey = "Save_Control_Histograms_bool";
	const std::string kAnnihTimeDiff2gcutParamKey = "AnnihTimeDiffcut2g_float";
        const std::string kVxyCut2gParamKey = "AnnihilationPointXYCut_float";
	const std::string kVzCut2gParamKey = "AnnihilationPointZCut_float";
	void saveEvents(const std::vector<JPetEvent>& event);
	double fScatterTOFTimeDiff = 2000.0;
	double fB2BSlotThetaDiff = 3.0;
	double fDeexTOTCutMin = 30000.0;
	double fDeexTOTCutMax = 50000.0;
	bool fSaveControlHistos = true;
	void initialiseHistograms();
	double fAnnihTimeDiffcut2g=666.; //ns
	double fVxyCut2g=666; //cm
	double fVzCut2g=666; //cm 
//
};
#endif /* !pPsto3gamma_H */
