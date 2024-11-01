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
/// \file Hadr03.cc
/// \brief Main program of the hadronic/Hadr03 example

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "G4ParticleHPManager.hh"
#include "CLI11.hpp"

int main(int argc, char** argv) {

  CLI::App app{"Neutron studies"};
  argv = app.ensure_utf8(argv);

  std::string filename = "";
  app.add_option("-f,--file", filename, "The G4 mac file")->required();

  int n_thread = 1;
  app.add_option("-t,--n-thread", n_thread, "The number of thread to run with");

  ActionType action_type = kFirstInteraction;
  std::map<std::string, ActionType> map{
    {"FirstInteraction", ActionType::kFirstInteraction},
    {"NeutronInteractionDistance", ActionType::kNeutronInteractionDistance},
  };

  app.add_option("-s,--study-type", action_type, "Which study to run")
    ->transform(CLI::CheckedTransformer(map, CLI::ignore_case));

  CLI11_PARSE(app, argc, argv);

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  //construct the run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(n_thread);

  //set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);

  PhysicsList* phys = new PhysicsList;
  runManager->SetUserInitialization(phys);
  runManager->SetUserInitialization(
    new ActionInitialization(det,action_type)
  );

  // Replaced HP environmental variables with C++ calls
  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes( false );
  G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState( true );
  G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation( true );
  G4ParticleHPManager::GetInstance()->SetNeglectDoppler( false );
  G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( true );
  G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel( false );
  G4ParticleHPManager::GetInstance()->SetUseNRESP71Model( false );

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();


  //batch mode
  std::string command = "/control/execute ";
  UImanager->ApplyCommand(command+filename);

  delete runManager;
}