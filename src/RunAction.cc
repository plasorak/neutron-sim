//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
#include "RunAction.hh"
#include "Run.hh"
#include "RunMessenger.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>


RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim, ActionType action_type)
  : fDetector(det)
  , fPrimary(prim)
  , fActionType(action_type)
{
 fHistoManager = new HistoManager();
 fRunMessenger = new RunMessenger(this);
}


RunAction::~RunAction()
{
 delete fHistoManager;
 delete fRunMessenger;
}


G4Run* RunAction::GenerateRun()
{
  switch (fActionType) {
    case ActionType::kPrimaryInteraction:
      fRun = new RunPrimaryInteraction(fDetector);
      return fRun;
    case ActionType::kCaptureDistance:
      fRun = new RunCaptureDistance(fDetector);
      return fRun;
    default:
      G4ExceptionDescription msg;
      msg << "Unknown action type: " << static_cast<int>(fActionType);
      G4Exception("RunAction::GenerateRun()", "RunAction001", FatalException, msg);
      return nullptr;
  }

}


void RunAction::BeginOfRunAction(const G4Run*)
{
  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();

  // keep run condition
  if (fPrimary) {
    G4ParticleDefinition* particle
      = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }

  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if (isMaster) fRun->EndOfRun(fPrint);

  //save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }

  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetPrintFlag(G4bool flag)
{
  fPrint = flag;
}

