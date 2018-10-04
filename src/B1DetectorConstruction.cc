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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "G4Transform3D.hh"

#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;


  //Dymensions of the smallest cylinder (Xenon)
  
  G4double Xe_rmin =  0.*m, Xe_rmax = 3.*m;
  G4double Xe_hz = 3.*m;
  G4double Xe_phimin = 0.*deg, Xe_phimax = 360.*deg;

  G4double Poly_thick = 0.03*m; //thickness Polyethylene
  G4double Cu_thick = 0.12*m;  //thickness Copper
  G4double tank2_thick = 0.03*m;
  G4double water_thick = 3*m;
  G4double tank1_thick = 0.03*m;

  G4double world_thick = 0.01*m;
  

  //     
  // World
  //
  G4double world_sizeXY =  2*(Xe_rmax + Poly_thick + Cu_thick + tank2_thick + water_thick + tank1_thick + world_thick);
  G4double world_sizeZ  =  2*(Xe_hz + Poly_thick + Cu_thick + tank2_thick + water_thick + tank1_thick + world_thick);
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  
  //     
  // 1. Steel tank
  //


  //Steel definition
  G4double density;
  G4int ncomponents;
  G4double fractionmass;
  G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
  StainlessSteel->AddMaterial(nist->FindOrBuildMaterial("G4_C"), fractionmass=0.001);
  StainlessSteel->AddMaterial(nist->FindOrBuildMaterial("G4_Si"), fractionmass=0.007);
  StainlessSteel->AddMaterial(nist->FindOrBuildMaterial("G4_Cr"), fractionmass=0.18);
  StainlessSteel->AddMaterial(nist->FindOrBuildMaterial("G4_Mn"), fractionmass=0.01);
  StainlessSteel->AddMaterial(nist->FindOrBuildMaterial("G4_Fe"), fractionmass=0.712);
  StainlessSteel->AddMaterial(nist->FindOrBuildMaterial("G4_Ni"), fractionmass=0.09);


  //Material of the cointainer 
  G4Material* tank1_mat = StainlessSteel;

  //Position
  //Not included in the rest of the sections, because all cilinders are centered in 0
  G4ThreeVector pos = G4ThreeVector(0, 0, 0);
  
  //Logical Volume       
  G4double tank1_rmin =  0.*m, tank1_rmax = Xe_rmax + Poly_thick + Cu_thick + tank2_thick + water_thick + tank1_thick;
  G4double tank1_hz = Xe_hz + Poly_thick + Cu_thick + tank2_thick + water_thick + tank1_thick;
  G4double tank1_phimin = 0.*deg, tank1_phimax = 360.*deg;
  G4Tubs* solidTank1 =    
    new G4Tubs("Steel", 
    tank1_rmin, tank1_rmax, tank1_hz,
    tank1_phimin, tank1_phimax);
                      
  G4LogicalVolume* logicTank1 =                         
    new G4LogicalVolume(solidTank1,         //its solid
                        tank1_mat,          //its material
                        "Tank1");           //its name
  //Physical Volume
               
  new G4PVPlacement(0,                       //no rotation
		    pos,                    //its position
                    logicTank1,             //its logical volume
                    "Tank1",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // 2. Water
  //

  //Material 
  G4Material* water_mat = nist->FindOrBuildMaterial("G4_WATER"); 
    
  //Logical Volume
  G4double water_rmin =  0.*m, water_rmax = Xe_rmax + Poly_thick + Cu_thick + tank2_thick + water_thick;
  G4double water_hz = Xe_hz + Poly_thick + Cu_thick + tank2_thick + water_thick;
  G4double water_phimin = 0.*deg, water_phimax = 360.*deg;
  G4Tubs* solidWater =    
    new G4Tubs("Water", 
    water_rmin,water_rmax, water_hz,
    water_phimin, water_phimax);
                      
  G4LogicalVolume* logicWater =                         
    new G4LogicalVolume(solidWater,         //its solid
                       water_mat,          //its material
                        "Water");           //its name

  //Physical Volume
               
  new G4PVPlacement(0,                       //no rotation
		    pos,                    //its position
                    logicWater,             //its logical volume
                    "Water",                //its name
                    logicTank1,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // 3. Second tank of steel. Pressure Vessel
  //

  //Material 
  G4Material* tank2_mat = StainlessSteel; 
    
  //Logical Volume
  G4double tank2_rmin =  0.*m, tank2_rmax = Xe_rmax + Poly_thick + Cu_thick + tank2_thick;
  G4double tank2_hz = Xe_hz + Poly_thick + Cu_thick + tank2_thick;
  G4double tank2_phimin = 0.*deg, tank2_phimax = 360.*deg;
  G4Tubs* solidTank2 =    
    new G4Tubs("Tank2", 
    tank2_rmin,tank2_rmax, tank2_hz,
    tank2_phimin, tank2_phimax);
                      
  G4LogicalVolume* logicTank2 =                         
    new G4LogicalVolume(solidTank2,         //its solid
                       tank2_mat,          //its material
                        "Tank2");           //its name

  //Physical Volume
               
  new G4PVPlacement(0,                       //no rotation
		    pos,                    //its position
                    logicTank2,             //its logical volume
                    "Tank2",                //its name
                    logicWater,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // 4. Copper
  //

  //Material 
  G4Material* Cu_mat = nist->FindOrBuildMaterial("G4_Cu");
    
  //Logical Volume
  G4double Cu_rmin =  0.*m, Cu_rmax = Xe_rmax + Poly_thick + Cu_thick;
  G4double Cu_hz = Xe_hz + Poly_thick + Cu_thick;
  G4double Cu_phimin = 0.*deg, Cu_phimax = 360.*deg;
  G4Tubs* solidCu =    
    new G4Tubs("Copper", 
    Cu_rmin,Cu_rmax, Cu_hz,
    Cu_phimin, Cu_phimax);
                      
  G4LogicalVolume* logicCu =                         
    new G4LogicalVolume(solidCu,         //its solid
                       Cu_mat,          //its material
                        "Copper");           //its name

  //Physical Volume
               
  new G4PVPlacement(0,                       //no rotation
		    pos,                    //its position
                    logicCu,             //its logical volume
                    "Copper",                //its name
                    logicTank2,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //
  // 5. Polyethylene
  //


  //Material 
  G4Material* Poly_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    
  //Logical Volume

  G4double Poly_rmin =  0.*m, Poly_rmax = Xe_rmax +  Poly_thick;
  G4double Poly_hz = Xe_hz + Poly_thick;
  G4double Poly_phimin = 0.*deg, Poly_phimax = 360.*deg;
  G4Tubs* solidPoly =    
    new G4Tubs("Polyethylene", 
    Poly_rmin,Poly_rmax, Poly_hz,
    Poly_phimin, Poly_phimax);
                      
  G4LogicalVolume* logicPoly =                         
    new G4LogicalVolume(solidPoly,         //its solid
                       Poly_mat,          //its material
                        "Polyethylene");           //its name

  //Physical Volume
               
  new G4PVPlacement(0,                       //no rotation
		    pos,                    //its position
                    logicPoly,             //its logical volume
                    "Polyethylene",                //its name
                    logicCu,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  


  //
  // 6. Xenon TPC
  //

  // Definition of Xe136

  G4double a; //mass of mole
  G4double z; //atomic number
  G4Material* GXe  =
    new G4Material("Xenon136",z=54., a = 135.907227*g/mole, density = 0.08663*g/cm3);

  
  //Material 
  G4Material* Xe_mat = GXe;

  
  //Logical Volume       

  G4Tubs* solidXe =    
    new G4Tubs("Xenon", 
    Xe_rmin,Xe_rmax, Xe_hz,
    Xe_phimin, Xe_phimax);
                      
  G4LogicalVolume* logicXe =                         
    new G4LogicalVolume(solidXe,         //its solid
                       Xe_mat,          //its material
                        "Xenon");           //its name

  //Physical Volume
               
  new G4PVPlacement(0,                       //no rotation
		    pos,                    //its position
                    logicXe,             //its logical volume
                    "Xenon",                //its name
                    logicPoly,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  


  // Set Xenon as scoring volume
  //
  fScoringVolume = logicXe ;

  
  //
  //always return the physical World
  //
  return physWorld;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
