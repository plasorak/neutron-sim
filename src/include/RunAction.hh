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
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#pragma once

#include "G4UserRunAction.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>
#include "ActionInitialization.hh"

class DetectorConstruction;
class Run;
class RunMessenger;
class PrimaryGeneratorAction;
class HistoManager;
class G4Run;

class RunAction : public G4UserRunAction {

  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*, ActionType);
   ~RunAction() override;

  public:
    G4Run* GenerateRun() override;
    void BeginOfRunAction(const G4Run*) override;
    void   EndOfRunAction(const G4Run*) override;

    void SetPrintFlag(G4bool);

  private:
    DetectorConstruction*      fDetector     = nullptr;
    PrimaryGeneratorAction*    fPrimary      = nullptr;
    Run*                       fRun          = nullptr;
    HistoManager*              fHistoManager = nullptr;
    RunMessenger*              fRunMessenger = nullptr;

    ActionType fActionType = ActionType::kPrimaryInteraction;
    G4bool   fPrint = true;      //optional printing
};

