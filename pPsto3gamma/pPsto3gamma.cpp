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
 *  @file pPsto3gamma.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "pPsto3gamma.h"
#include "pPsto3gammaTools.h"
#include <iostream>

using namespace jpet_options_tools;
using namespace std;

pPsto3gamma::pPsto3gamma(const char* name): JPetUserTask(name) {}

pPsto3gamma::~pPsto3gamma() {}

bool pPsto3gamma::init()
{
  INFO("pPs->3g analysis started.");
  // Parameter for back to back categorization
  if (isOptionSet(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey)){
    fB2BSlotThetaDiff = getOptionAsFloat(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kBack2BackSlotThetaDiffParamKey.c_str(), fB2BSlotThetaDiff
    ));
  }
  // Parameter for scattering determination
  if (isOptionSet(fParams.getOptions(), kScatterTOFTimeDiffParamKey)) {
    fScatterTOFTimeDiff = getOptionAsFloat(fParams.getOptions(), kScatterTOFTimeDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kScatterTOFTimeDiffParamKey.c_str(), fScatterTOFTimeDiff
    ));
  }
  // Parameters for deexcitation TOT cut
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMinParamKey)) {
    fDeexTOTCutMin = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMinParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMinParamKey.c_str(), fDeexTOTCutMin
    ));
  }
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMaxParamKey)) {
    fDeexTOTCutMax = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMaxParamKey.c_str(), fDeexTOTCutMax
    ));
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  //2 gamma annihilation time cut
  if (isOptionSet(fParams.getOptions(),kAnnihTimeDiff2gcutParamKey)) {
   fAnnihTimeDiffcut2g = getOptionAsFloat(fParams.getOptions(), kAnnihTimeDiff2gcutParamKey);
  } else {
    WARNING(Form(
		 "No value of the %s parameter provided by the user. Using default value of %lf.",
		 kAnnihTimeDiff2gcutParamKey.c_str(),fAnnihTimeDiffcut2g
		 ));
  }
    //2 gamma vetex cuts
    if (isOptionSet(fParams.getOptions(),kVxyCut2gParamKey)) {
      fVxyCut2g = getOptionAsFloat(fParams.getOptions(), kVxyCut2gParamKey);
    } else{
      WARNING(Form(
		   "No value of the %s parameter provided by the user. Using default value of %lf.",
		   kVxyCut2gParamKey.c_str(),fVxyCut2g
		   ));
    }
    if (isOptionSet(fParams.getOptions(),kVzCut2gParamKey)) {
      fVzCut2g = getOptionAsFloat(fParams.getOptions(), kVzCut2gParamKey);
    } else{
      WARNING(Form(
                   "No value of the %s parameter provided by the user. Using default value of %lf.",
                   kVzCut2gParamKey.c_str(),fVzCut2g
                   ));
    }
  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  if(fSaveControlHistos) initialiseHistograms();
  return true;
}

bool pPsto3gamma::exec()
{
  std::vector <int>  PromptCandidatesId;
  std::vector <int> TwoGammaCandidatesId;//first two is the first pair, the  second two the second pair etc.
  std::vector<double> AnnihTimesCandidates;
  int Nhit =0;
  int N2GammaPairs=0;
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));
      Nhit = event.getHits().size();
      getStatistics().getHisto1D("NhitPerEvent")->Fill(Nhit);
      //Looking for the de-excitation photon
      int isPrompt = pPsto3gammaTools::checkForPrompt(
	  event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax,PromptCandidatesId );
      getStatistics().getHisto1D("NpromptAll")->Fill(isPrompt);
      if(isPrompt){
	//The normalization sample
	if(Nhit==3){
	  double LifeTime2g = -666.;
          int flags2g[3] = {-1,-1,-1};
          double AnnihTime = -666.;
	 	  N2GammaPairs = pPsto3gammaTools::checkFor2Gamma(                                                                                                                                         
	  						  event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff,fAnnihTimeDiffcut2g,
	  						  fVxyCut2g,fVzCut2g,TwoGammaCandidatesId,AnnihTimesCandidates);
		  std::cout << "\n N2gammaPairs= "<<N2GammaPairs << std::endl;
	  for (int i = 0; i < N2GammaPairs; i++) {
       	    JPetHit firstHit, secondHit,PromptHit;
	    firstHit = event.getHits().at(TwoGammaCandidatesId[2*i]);
	    secondHit = event.getHits().at(TwoGammaCandidatesId[2*i+1]);
	    TVector3 Vtx = pPsto3gammaTools::calculateAnnihilationPoint(firstHit, secondHit);
	    double DeexcTime=100000000.;
	    for (int j = 0; j < isPrompt; j++) {
	      if(!(PromptCandidatesId[j]==TwoGammaCandidatesId[2*i] 
		   || PromptCandidatesId[j]==TwoGammaCandidatesId[2*i+1])){
	        PromptHit = event.getHits().at(TwoGammaCandidatesId[j]);
	        double DeexcTimeTmp = pPsto3gammaTools::calculateAnnihilationTime(PromptHit,Vtx);
		AnnihTime = (AnnihTimesCandidates[2*i] + AnnihTimesCandidates[2*i+1])/2.;
		if(DeexcTimeTmp < DeexcTime && DeexcTimeTmp < AnnihTime){
		  DeexcTime = DeexcTimeTmp;		  
		  LifeTime2g = AnnihTime -DeexcTime;
		  flags2g[0] = PromptCandidatesId[j];
		  flags2g[1] = TwoGammaCandidatesId[2*i];
		  flags2g[2] = TwoGammaCandidatesId[2*i+1];
		}	    
	      }
	    }
	  }
	  //
	  //	  std::cout <<flags2g <<std::endl; 
	  getStatistics().getHisto1D("2glifetime")->Fill(LifeTime2g);
	  getStatistics().getHisto1D("Nprompt2g")->Fill(isPrompt);
	}
	if(Nhit==4){
	  getStatistics().getHisto1D("Nprompt3g")->Fill(isPrompt);
	}
      }
      // Check types of current event

      /*      bool is2Gamma = EventCategorizerTools::checkFor2Gamma(
        event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff
      );
      bool is3Gamma = EventCategorizerTools::checkFor3Gamma(
        event, getStatistics(), fSaveControlHistos
      );
      bool isPrompt = EventCategorizerTools::checkForPrompt(
        event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax
      );
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff
      );
      */
      JPetEvent newEvent = event;
      /*
      if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
      if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);
      */
      if(fSaveControlHistos){
        for(auto hit : event.getHits()){
          getStatistics().getHisto2D("All_XYpos")->Fill(hit.getPosX(), hit.getPosY());
        }
      }
      events.push_back(newEvent);
    }
    saveEvents(events);
  } else { return false; }
  return true;
}

bool pPsto3gamma::terminate()
{
  INFO("Event categorization completed.");
  return true;
}

void pPsto3gamma::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void pPsto3gamma::initialiseHistograms(){

  // General histograms
  getStatistics().createHistogram(
    new TH2F("All_XYpos", "Hit position XY", 242, -60.5, 60.5, 121, -60.5, 60.5));
  getStatistics().getHisto2D("All_XYpos")->GetXaxis()->SetTitle("Hit X position [cm]");
  getStatistics().getHisto2D("All_XYpos")->GetYaxis()->SetTitle("Hit Y position [cm]");

  // Histograms for 2Gamama category
  getStatistics().createHistogram(
    new TH1F("2Gamma_Zpos", "B2B hits Z position", 200, -50.0, 50.0));
  getStatistics().getHisto1D("2Gamma_Zpos")->GetXaxis()->SetTitle("Z axis position [cm]");
  getStatistics().getHisto1D("2Gamma_Zpos")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(
    new TH1F("2Gamma_TimeDiff", "B2B hits time difference", 100, -10000.0, 10000.0));
  getStatistics().getHisto1D("2Gamma_TimeDiff")->GetXaxis()->SetTitle("Time Difference [ps]");
  getStatistics().getHisto1D("2Gamma_TimeDiff")->GetYaxis()->SetTitle("Number of Hit Pairs");

  getStatistics().createHistogram(
    new TH1F("2Gamma_Dist", "B2B hits distance", 150, 0.0, 150.0));
  getStatistics().getHisto1D("2Gamma_Dist")->GetXaxis()->SetTitle("Distance [cm]");
  getStatistics().getHisto1D("2Gamma_Dist")->GetYaxis()->SetTitle("Number of Hit Pairs");

  getStatistics().createHistogram(
    new TH1F("Annih_TOF", "Annihilation pairs Time of Flight", 200, -3000.0, 3000.0));
  getStatistics().getHisto1D("Annih_TOF")->GetXaxis()->SetTitle("Time of Flight [ps]");
  getStatistics().getHisto1D("Annih_TOF")->GetYaxis()->SetTitle("Number of Annihilation Pairs");

  getStatistics().createHistogram(
     new TH2F("AnnihPoint_XY", "XY position of annihilation point", 121, -60.5, 60.5, 121, -60.5, 60.5));
  getStatistics().getHisto2D("AnnihPoint_XY")->GetXaxis()->SetTitle("X position [cm]");
  getStatistics().getHisto2D("AnnihPoint_XY")->GetYaxis()->SetTitle("Y position [cm]");

  getStatistics().createHistogram(
    new TH2F("AnnihPoint_XZ", "XZ position of annihilation point", 121, -60.5, 60.5, 121, -60.5, 60.5));
  getStatistics().getHisto2D("AnnihPoint_XZ")->GetXaxis()->SetTitle("X position [cm]");
  getStatistics().getHisto2D("AnnihPoint_XZ")->GetYaxis()->SetTitle("Z position [cm]");

  getStatistics().createHistogram(
    new TH2F("AnnihPoint_YZ", "YZ position of annihilation point", 121, -60.5, 60.5, 121, -60.5, 60.5));
  getStatistics().getHisto2D("AnnihPoint_YZ")->GetXaxis()->SetTitle("Y position [cm]");
  getStatistics().getHisto2D("AnnihPoint_YZ")->GetYaxis()->SetTitle("Z position [cm]");

  // Histograms for 3Gamama category
  getStatistics().createHistogram(
    new TH2F("3Gamma_Angles", "Relative angles - transformed", 251, -0.5, 250.5, 201, -0.5, 200.5));
  getStatistics().getHisto2D("3Gamma_Angles")->GetXaxis()->SetTitle("Relative angle 1-2");
  getStatistics().getHisto2D("3Gamma_Angles")->GetYaxis()->SetTitle("Relative angle 2-3");

  // Histograms for scattering category
  getStatistics().createHistogram(
    new TH1F("ScatterTOF_TimeDiff", "Difference of Scatter TOF and hits time difference",
      200, 0.0, 3.0*fScatterTOFTimeDiff));
  getStatistics().getHisto1D("ScatterTOF_TimeDiff")->GetXaxis()->SetTitle("Scat_TOF & time diff [ps]");
  getStatistics().getHisto1D("ScatterTOF_TimeDiff")->GetYaxis()->SetTitle("Number of Hit Pairs");

  getStatistics().createHistogram(
     new TH2F("ScatterAngle_PrimaryTOT", "Angle of scattering vs. TOT of primary hits",
      181, -0.5, 180.5, 200, 0.0, 40000.0));
  getStatistics().getHisto2D("ScatterAngle_PrimaryTOT")->GetXaxis()->SetTitle("Scattering Angle");
  getStatistics().getHisto2D("ScatterAngle_PrimaryTOT")->GetYaxis()->SetTitle("TOT of primary hit [ps]");

  getStatistics().createHistogram(
     new TH2F("ScatterAngle_ScatterTOT", "Angle of scattering vs. TOT of scattered hits",
      181, -0.5, 180.5, 200, 0.0, 40000.0));
  getStatistics().getHisto2D("ScatterAngle_ScatterTOT")->GetXaxis()->SetTitle("Scattering Angle");
  getStatistics().getHisto2D("ScatterAngle_ScatterTOT")->GetYaxis()->SetTitle("TOT of scattered hit [ps]");

  // Histograms for deexcitation
  getStatistics().createHistogram(
    new TH1F("Deex_TOT_cut", "TOT of all hits with deex cut (30,50) ns",
      200, 25000.0, 55000.0));
  getStatistics().getHisto1D("Deex_TOT_cut")->GetXaxis()->SetTitle("TOT [ps]");
  getStatistics().getHisto1D("Deex_TOT_cut")->GetYaxis()->SetTitle("Number of Hits");
  //
  //New histograms for 2 gamma
  //
  getStatistics().createHistogram(
				  new TH1F("2Gamma_AnnihTimeDiff", "Diference of the emmision time of 2 gamma hits ",
					   200, -100.0, 100.0));
 getStatistics().getHisto1D("2Gamma_AnnihTimeDiff")->GetXaxis()->SetTitle("Delta t [ns]");
 getStatistics().getHisto1D("2Gamma_AnnihTimeDiff")->GetYaxis()->SetTitle("Counts");
 //
 getStatistics().createHistogram(
				 new TH1F("NpromptAll", "Number of prompts per event",
					  11, -1.0, 10.0));
 getStatistics().getHisto1D("NpromptAll")->GetXaxis()->SetTitle("Number of prompt photons");
 getStatistics().getHisto1D("NpromptAll")->GetYaxis()->SetTitle("Counts");
 //
 getStatistics().createHistogram(
                                 new TH1F("Nprompt2g", "Number of prompts per event for 2g candidate events",
                                          11, -1.0, 10.0));
 getStatistics().getHisto1D("Nprompt2g")->GetXaxis()->SetTitle("Number of prompt photons");
 getStatistics().getHisto1D("Nprompt2g")->GetYaxis()->SetTitle("Counts");
 //                                                                                                                                                                                                          
 getStatistics().createHistogram(
                                 new TH1F("Nprompt3g", "Number of prompts per event for 3g candidate events",
                                          11, -1.0, 10.0));
 getStatistics().getHisto1D("Nprompt3g")->GetXaxis()->SetTitle("Number of prompt photons");
 getStatistics().getHisto1D("Nprompt3g")->GetYaxis()->SetTitle("Counts");
 //
 getStatistics().createHistogram(
                                 new TH1F("2glifetime", "Lifetime for the 2g events",
                                          8000, -200.0, 600.0));
 getStatistics().getHisto1D("2glifetime")->GetXaxis()->SetTitle("Positronium lifetime");
 getStatistics().getHisto1D("2glifetime")->GetYaxis()->SetTitle("Counts");
 //
 getStatistics().createHistogram(
                                 new TH1F("NhitPerEvent", "Number of hits per event",
                                          11, -1.0, 10.0));
 getStatistics().getHisto1D("NhitPerEvent")->GetXaxis()->SetTitle("Number of hits");
 getStatistics().getHisto1D("NhitPerEvent")->GetYaxis()->SetTitle("Counts");
}
