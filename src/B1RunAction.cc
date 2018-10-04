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
// $Id: B1RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "vector"
#include "string"

#include <TH1F.h>
#include <TH2D.h>
#include "TFile.h"
#include <sstream>

#include "G4ThreeVector.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)

  {

   //MY ADDITION
    i = 0;

  // add new units for dose
  //


  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager 
 G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  //Initiate the number of relevant events in the run
  i = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B1RunAction::Histogram(G4int n)
{
  G4int j;
  //std::string Title  = "E-" +  std::to_string(n) +"-Cu.root" ;
  //TH1D *h1 = new TH1D ("h1", Title.c_str(), 100,0., 3.1);

  TH1D *h1 = new TH1D ("h1","10^6 from Cu", 100,0., 3.1);
  for (j = 0; j < i; j++)
    {
      h1->Fill(Evector[j]);
    }
  //FIGURE OUT WHY THIS DOESN'T WORK!
  //TFile* File = new TFile("E-" +  std::to_string(n) +"-Cu.root", "recreate");
  TFile* File = new TFile("TestSingleH3.root ", "recreate");
  h1->Write();

  delete h1;
  //File->Close();
}

void B1RunAction::HistogramClassified(G4int n)
{
  TH1D *hCompton = new TH1D ("hCompton", "10^6 Gammas producing Compton", 100,0., 3.1);
  TH1D *hPairProd = new TH1D ("hPairProd", "10^6 Gammas producing Pair Production", 100,0., 3.1);
  TH1D *hPhot = new TH1D ("hPhot", "10^6 Gammas producing Photoelectric Effect", 100,0., 3.1);
  TH1D *hFirstGamma = new TH1D ("hFirstGamma", "E distribution of the first gamma arriving to the Xenon", 100,0., 3.1);
  for (std::vector<G4double>::iterator it = EComptVector.begin(); it != EComptVector.end(); it++)
    {
      hCompton->Fill(*it);
    }

  for (std::vector<G4double>::iterator it = EPairProdVector.begin(); it != EPairProdVector.end(); it++)
    {
      hPairProd->Fill(*it);
    }
  
  for (std::vector<G4double>::iterator it = EPhotoVector.begin(); it != EPhotoVector.end(); it++)
    {
      hPhot->Fill(*it);
    }

    for (std::vector<G4double>::iterator it = EFirstGammaVector.begin(); it != EFirstGammaVector.end(); it++)
    {
      hFirstGamma->Fill(*it);
    }
  TFile* File = new TFile("TestSeparatedH3.root", "recreate");
  
  hCompton->Write();
  hPairProd->Write();
  hPhot->Write();
  hFirstGamma->Write();

  delete hCompton;
  delete hPairProd;
  delete hPhot;
  delete hFirstGamma;
  //File->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::HistogramPosition()
{
  TH2D *hXY = new TH2D ("hXY", "XY positions for the first E deposition in Xenon", 1000,-3300., 3300., 1000,-3300., 3300.);
  TH2D *hXZ = new TH2D ("hXZ", "XZ positions for the first E deposition in Xenon", 1000,-3300., 3300., 1000,-3300., 3300.);
  TH2D *hYZ = new TH2D ("hYZ", "YZ positions for the first E deposition in Xenon", 1000,-3300., 3300., 1000,-3300., 3300.);

  for (std::vector<G4ThreeVector>::iterator it = PositionFirstGammaVector.begin(); it != PositionFirstGammaVector.end(); it++)
    {
      hXY->Fill(it->getX(), it->getY());
      hXZ->Fill(it->getX(), it->getZ());
      hYZ->Fill(it->getY(), it->getZ());
    }
  
  TFile* File = new TFile("TestPositionH3.root", "recreate");
  
  hXY->Write();
  hXZ->Write();
  hYZ->Write();

  delete hXY;
  delete hXZ;
  delete hYZ;
  //File->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B1RunAction::SaveFile(G4int n)
{
  std::ofstream f;
  G4int j;

  f.open ("E-" + std::to_string(n) + "-CuTest.txt");
  
  
  for (j = 0; j < i; j++)
    {
      f << Evector[j] << std:: endl;
    }
  f.close();
  
}



void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();

  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
  /*
  G4int j = 0;
  
  G4cout << "EVENT                    E deposited (MeV) " << G4endl;

  for (j = 0; j < i; j++)
    {
       G4cout << j << "                          " << Evector[j] << G4endl;
    }
  G4cout << "Number of relevant events " << i << G4endl;
  */
 
  B1RunAction::HistogramClassified(nofEvents);
  B1RunAction::Histogram(nofEvents);
  B1RunAction::HistogramPosition();

  G4double totalE = 0.;

  for (std::vector<G4double>::iterator it = Evector.begin(); it != Evector.end(); it++)
    {
      totalE += *it;
    }

  // B1RunAction::SaveFile(nofEvents);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::Eevent(G4double edep)
{
  Evector.push_back(edep);
  i++;
}

void B1RunAction::Eclassified(G4double edep, G4String ProcessName)
{
  if(ProcessName.compareTo("compt") == 0) {EComptVector.push_back(edep);}

  else if(ProcessName.compareTo("ePairProd") == 0) {EPairProdVector.push_back(edep);}

  else if(ProcessName.compareTo("phot") == 0)
    {
      EPhotoVector.push_back(edep);
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EFirstGamma (G4double EGamma)
{
  EFirstGammaVector.push_back(EGamma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::PositionFirstGamma (G4ThreeVector position)
{
  PositionFirstGammaVector.push_back(position);
}
