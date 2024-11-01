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
/// \file Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4ProcessTable.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicProcess.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"

#include <iostream>
#include <map>
#include <numeric>
#include <iterator>

RunFirstInteraction::RunFirstInteraction(DetectorConstruction* det)
: fDetector(det) {

  for (G4int i=0; i<3; i++)
    fPbalance[i] = 0.;

  fPbalance[1] = DBL_MAX;
}


void RunFirstInteraction::CountProcesses(G4VProcess* process) {

  if (process == nullptr)
    return;

  G4String procName = process->GetProcessName();

  auto it = fProcCounter.find(procName);

  if (it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  } else {
    fProcCounter[procName]++;
  }

}


void RunFirstInteraction::SumTrack(G4double trackl, G4VProcess* process) {

  if (process == nullptr)
    return;

  G4String procName = process->GetProcessName();
  auto it = fSumTrack.find(procName);

  if (it == fSumTrack.end()) {
    fSumTrack[procName] = trackl;
    fSumTrack2[procName] = trackl*trackl;
  }  else {
    fSumTrack[procName] += trackl;
    fSumTrack2[procName] += trackl*trackl;

  }

}


void RunFirstInteraction::CountNuclearChannel(G4String name, G4double Q) {

  auto it = fNuclChannelMap.find(name);

  if (it == fNuclChannelMap.end()) {
    fNuclChannelMap[name] = NuclChannel(1, Q);
  } else {
    NuclChannel& data = it->second;
    data.fCount++;
    data.fQ += Q;
  }
}


void RunFirstInteraction::ParticleCount(G4String name, G4double Ekin) {

  auto it = fParticleDataMap.find(name);

  if (it == fParticleDataMap.end()) {
    fParticleDataMap[name] = ParticleData(1, Ekin, Ekin, Ekin);
  } else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin;
  }
}


void RunFirstInteraction::Merge(const G4Run* run) {

  auto const* localRun = static_cast<const RunFirstInteraction*>(run);

  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  for (auto const& st: localRun->fSumTrack){
    auto procName = st.first;
    auto localCount = st.second;
    if (fSumTrack.find(procName) == fSumTrack.end()) fSumTrack[procName] = localCount;
    else                                             fSumTrack[procName] += localCount;
  }

  for (auto const& st: localRun->fSumTrack2){
    auto procName = st.first;
    auto localCount = st.second;
    if (fSumTrack2.find(procName) == fSumTrack2.end()) fSumTrack2[procName] = localCount;
    else                                               fSumTrack2[procName] += localCount;
  }

  for (auto const& itp: localRun->fProcCounter){
    auto procName = itp.first;
    auto localCount = itp.second;
    if (fProcCounter.find(procName) == fProcCounter.end()) fProcCounter[procName] = localCount;
    else                                                   fProcCounter[procName] += localCount;
  }

  for (auto const& itc : localRun->fNuclChannelMap){
    G4String name = itc.first;
    const NuclChannel& localData = itc.second;

    if (fNuclChannelMap.find(name) == fNuclChannelMap.end()) {
      fNuclChannelMap[name] = NuclChannel(localData.fCount, localData.fQ);
    } else {
      NuclChannel& data = fNuclChannelMap[name];
      data.fCount += localData.fCount;
      data.fQ     += localData.fQ;
    }
  }

  for (auto const& itn: localRun->fParticleDataMap) {
    G4String name = itn.first;
    const ParticleData& localData = itn.second;

    if (fParticleDataMap.find(name) == fParticleDataMap.end()) {
      fParticleDataMap[name] = ParticleData(
        localData.fCount,
        localData.fEmean,
        localData.fEmin,
        localData.fEmax
      );
    } else {
      ParticleData& data = fParticleDataMap[name];
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin)
        data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax)
        data.fEmax = emax;
    }
  }

  G4Run::Merge(run);
}

void RunFirstInteraction::EndOfRun(G4bool print) {

  auto e = std::string(G4BestUnit(fEkin,"Energy"));
  while (e.find(" ") != std::string::npos)
    e = e.replace(e.find(" "), 1, "");

  std::ofstream data_file("neutron_"+std::string(e)+".txt");

  //run condition
  const G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();

  G4String Particle = fParticle->GetParticleName();
  data_file << "The run is " << numberOfEvent << " "<< Particle << " of "
            << G4BestUnit(fEkin,"Energy") << " through "
            << G4BestUnit(fDetector->GetSize(),"Length") << " of "
            << material->GetName() << " (density: "
            << G4BestUnit(density,"Volumic Mass") << ")\n";

  if (numberOfEvent == 0) return;

  //frequency of processes
  data_file << "\nProcess calls frequency:\n";
  G4int survive = 0;
  for (auto const& it: fProcCounter){
     G4String procName = it.first;
     G4int count = it.second;
     data_file << procName << "," << count << "\n";
     if (procName == "Transportation") survive = count;
  }

  if (survive > 0) {
    data_file << "Nb of incident particles surviving after "
              << G4BestUnit(fDetector->GetSize(),"Length") << " of "
              << material->GetName() << " : " << survive << "\n";
  }

  //compute mean free path and related quantities
  data_file << "\nMean free paths / Cross sections:\n";

  auto total_interactions = std::accumulate(
    std::begin(fProcCounter),
    std::end(fProcCounter),
    0,
    [] (G4int value, const std::map<G4String,G4int>::value_type& p) { return value + p.second; }
  );

  for (auto const& st: fSumTrack) {
    auto proc = st.first;
    G4double MeanFreePath = st.second * total_interactions / fProcCounter[proc] / fProcCounter[proc];
    G4double MeanTrack2 = fSumTrack2[proc] * total_interactions / fProcCounter[proc] / fProcCounter[proc];
    G4double rms = std::sqrt(std::fabs(MeanTrack2 - MeanFreePath*MeanFreePath));

    G4double CrossSection = 0.0;

    if (MeanFreePath > 0.0)
      CrossSection = 1. / MeanFreePath;

    G4double massicMFP = MeanFreePath*density;
    G4double massicCS  = 0.0;

    if (massicMFP > 0.0)
      massicCS = 1. / massicMFP;

    data_file << st.first << ":\n"
              << "MeanFreePath: " << G4BestUnit(MeanFreePath,"Length") << "+-" << G4BestUnit(rms,"Length")
              << ", massic: " << G4BestUnit(massicMFP, "Mass/Surface")
              << "\nCrossSection: " << CrossSection*cm << " cm^-1 "
              << ", massic: " << G4BestUnit(massicCS, "Surface/Mass")
              << "\nTotalCount: " << fProcCounter[proc] << "\n\n";

  }

  data_file << "\nVerification: crossSections from G4HadronicProcessStore\n";

  G4ProcessTable* processTable  = G4ProcessTable::GetProcessTable();
  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  G4double sumc1 = 0.0, sumc2 = 0.0;

  const G4Element* element = (material->GetNumberOfElements() == 1)? material->GetElement(0) : nullptr;

  for (auto const& it: fProcCounter) {
    G4String procName = it.first;
    const G4VProcess* process = processTable->FindProcess(procName, fParticle);
    PrintXS(process, material, element, store, density, sumc1, sumc2, data_file);
  }

  if (sumc1 > 0.0) {
    data_file << "Total" << " = " << G4BestUnit(sumc1, "Surface/Mass") << " ";
    if (sumc2 > 0.0) {
      data_file << G4BestUnit(sumc2, "Surface");
    }
    data_file << "\n";
  } else {
    data_file << " not available\n";
  }

  data_file << "\nList of nuclear reactions:\n";
  for (auto const& ic: fNuclChannelMap){
    G4String name    = ic.first;
    NuclChannel data = ic.second;
    G4int count = data.fCount;
    G4double Q  = data.fQ/count;
    data_file << name
              << " (" << G4BestUnit(Q, "Energy") << "), "
              << count << "\n";
  }

  data_file << "\nList of generated particles:\n";

  for (auto const& itn: fParticleDataMap) {
    G4String name = itn.first;

    if (name != "neutron")
      continue;

    ParticleData data = itn.second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;
    data_file << name
              << "," << count
              << "," << G4BestUnit(eMean, "Energy")
              << "," << G4BestUnit(eMin, "Energy")
              << "," << G4BestUnit(eMax, "Energy")
              << "\n";
  }

  fProcCounter.clear();
  fNuclChannelMap.clear();
  fParticleDataMap.clear();

  data_file.close();
}

void RunFirstInteraction::PrintXS(
  const G4VProcess* proc,
  const G4Material* mat, const G4Element* elm,
  G4HadronicProcessStore* store, G4double density,
  G4double& sum1, G4double& sum2, std::ofstream& buffer) {

  if (proc == nullptr)
    return;

  G4double xs1 = store->GetCrossSectionPerVolume(fParticle, fEkin, proc, mat);
  G4double massSigma = xs1/density;
  sum1 += massSigma;

  if (nullptr != elm) {
    G4double xs2 = store->GetCrossSectionPerAtom(fParticle, fEkin, proc, elm, mat);
    sum2 += xs2;
    buffer << proc->GetProcessName() << " = " << G4BestUnit(massSigma, "Surface/Mass") << " " << G4BestUnit(xs2, "Surface") << "\n";
  } else {
    buffer << proc->GetProcessName() << " = " << G4BestUnit(massSigma, "Surface/Mass") << "\n";
  }
}

  void RunNeutronInteractionDistance::SaveDistance(const G4StepPoint*, const G4VProcess*, const std::string){}
  void RunNeutronInteractionDistance::Merge(const G4Run*) {}
  void RunNeutronInteractionDistance::EndOfRun(G4bool) {}
