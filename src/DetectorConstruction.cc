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
// DetectorConstruction.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "CalorimeterSD.hh"
#include "CheckScatteringSD.hh"
#include "DetectorMessenger.hh"
#include "StandardDetectorSD.hh"
#include "StepRecordSD.hh"
#include "TrackingDetectorSD.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4SDManager.hh"

#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include "G4VisAttributes.hh"

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4String conf) : G4VUserDetectorConstruction(), fConfig(conf)
{
    if (fConfig != "prad" && fConfig != "drad" && fConfig != "test")
        fConfig = "prad";

    fVisAtts.clear();

    fWorldSizeXY = 150.0 * cm;
    fWorldSizeZ = 800.0 * cm;

    fTargetCenter = -455.0 * cm; // PRad survey
    fTargetR = 14.5 * cm;;
    fTargetHalfL = 2.75 * cm;
    fTargetMat = "D2Gas";
    fTargetDensityRatio = 1.0;

    fRecoilDetNSeg = 20;
    fRecoilDetCenter = fTargetCenter;
    fRecoilDetR = 13.5 * cm;
    fRecoilDetHalfL = 2.6 * cm;
    fRecoilDetL1Thickness = 200 * um;
    fRecoilDetL2Thickness = 300 * um;

    fDownChamberCenter = fTargetCenter + 74.0 * mm + 71.0 * cm / 2.0;
    fVacBoxCenter = fTargetCenter + 74.0 * mm + 71.0 * cm + 425.17 * cm / 2.0;

    fGEMCenter[0] = 217.5 * cm;
    fGEMCenter[1] = ( 257.5 + 0.0) * cm;

    fSciPlaneCenter = 262.5 * cm;

    fCrystalSurf = (295.0 + 0.0) * cm;

    fExtDensityRatio = 1.0;

    if (fConfig == "drad") {
        fTargetSDOn = false;
        fRecoilDetSDOn = true;
        fGEMSDOn = true;
        fSciPlaneSDOn = true;
        fHyCalSDOn = true;
        fVirtualSDOn = false;
    } else {
        fTargetSDOn = false;
        fRecoilDetSDOn = false;
        fGEMSDOn = true;
        fSciPlaneSDOn = false;
        fHyCalSDOn = true;
        fVirtualSDOn = false;
    }

    fAttenuationLG = 0.0;
    fReflectanceLG = 1.0;

    detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();

    // Define volumes
    if (fConfig == "drad")
        return DefineDRadVolumes();
    else if (fConfig == "test")
        return DefineTestVolumes();
    else
        return DefinePRadVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    if (fConfig == "drad")
        DefineDRadSDs();
    else if (fConfig == "test")
        DefineTestSDs();
    else
        DefinePRadSDs();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    G4String symbol;
    G4int z, n;
    G4double a;
    G4double density;
    G4int ncomponents, natoms;
    G4double fractionmass;

    G4NistManager *pNM = G4NistManager::Instance();

    // Define elements from NIST material table
    G4Element *H  = pNM->FindOrBuildElement(z = 1);
    G4Element *He = pNM->FindOrBuildElement(z = 2);
    G4Element *C  = pNM->FindOrBuildElement(z = 6);
    G4Element *N  = pNM->FindOrBuildElement(z = 7);
    G4Element *O  = pNM->FindOrBuildElement(z = 8);
    G4Element *F  = pNM->FindOrBuildElement(z = 9);
    G4Element *Na = pNM->FindOrBuildElement(z = 11);
    G4Element *Al = pNM->FindOrBuildElement(z = 13);
    G4Element *Si = pNM->FindOrBuildElement(z = 14);
    G4Element *P  = pNM->FindOrBuildElement(z = 15);
    G4Element *S  = pNM->FindOrBuildElement(z = 16);
    G4Element *Ar = pNM->FindOrBuildElement(z = 18);
    G4Element *K  = pNM->FindOrBuildElement(z = 19);
    G4Element *Cr = pNM->FindOrBuildElement(z = 24);
    G4Element *Mn = pNM->FindOrBuildElement(z = 25);
    G4Element *Fe = pNM->FindOrBuildElement(z = 26);
    G4Element *Ni = pNM->FindOrBuildElement(z = 28);
    G4Element *Cu = pNM->FindOrBuildElement(z = 29);
    G4Element *Zn = pNM->FindOrBuildElement(z = 30);
    G4Element *As = pNM->FindOrBuildElement(z = 33);
    G4Element *Ta = pNM->FindOrBuildElement(z = 73);
    G4Element *W  = pNM->FindOrBuildElement(z = 74);
    G4Element *Pb = pNM->FindOrBuildElement(z = 82);
    
    //G4Material* H2LiquidG4 = pNM->FindOrBuildMaterial("G4_lH2", false, false);
    //G4cout<<H2Liquid->GetName()<<" "<<H2Liquid->GetDensity()<<" "<<H2Liquid->GetTemperature()<<" "<<H2Liquid->GetPressure()<<" "<<H2Liquid->GetRadlen()<<" "<<H2Liquid->GetNumberOfElements()<<G4endl;
    

    // Define isotopes
    G4Isotope *H2 = new G4Isotope("H2", z = 1, n = 2, a = 2.0141 * g / mole);
    G4Element *D = new G4Element("Deuterium", symbol = "D", ncomponents = 1);
    D->AddIsotope(H2, 1.0);

    // Define materials

    // Space Vacuum
    G4Material *Galaxy = new G4Material("Galaxy", density = universe_mean_density, ncomponents = 1, kStateGas, 0.1 * kelvin, 1.0e-19 * pascal);
    Galaxy->AddElement(H, fractionmass = 1.0);
    fVisAtts[Galaxy->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Air
    G4Material *Air = new G4Material("Air", density = 1.292 * mg / cm3, ncomponents = 2);
    Air->AddElement(N, fractionmass = 0.7);
    Air->AddElement(O, fractionmass = 0.3);
    fVisAtts[Air->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Air vacuum of 1.e-6 torr at room temperature, 1 atmosphere = 760 torr
    G4Material *Vacuum = new G4Material("Vacuum", density = 1.0e-6 / 760.0 * 1.292 * mg / cm3, ncomponents = 1, kStateGas, STP_Temperature, 1.0e-6 / 760.0 * atmosphere);
    Vacuum->AddMaterial(Air, fractionmass = 1.0);
    fVisAtts[Vacuum->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Hydrogen Gas (T = 19.5 K, P = 470 mTorr)
    G4Material *H2Gas = new G4Material("H2Gas", density = fTargetDensityRatio * 0.47 / 760.0 * 273.15 / 19.5 * 0.08988 * mg / cm3, ncomponents = 1, kStateGas, 19.5 * kelvin, fTargetDensityRatio * 0.47 / 760.0 * atmosphere);
    H2Gas->AddElement(H, natoms = 2);
    fVisAtts[H2Gas->GetName()] = new G4VisAttributes(G4Colour::Cyan());
    
    // Hydrogen liquid
    G4Material *H2Liquid = new G4Material("H2Liquid", density = fTargetDensityRatio * 70.80 * mg / cm3, ncomponents = 1, kStateLiquid, 20.0 * kelvin, fTargetDensityRatio * 1.135 * atmosphere);
    H2Liquid->AddElement(H, natoms = 2);
    fVisAtts[H2Liquid->GetName()] = new G4VisAttributes(G4Colour::Cyan());
    fVisAtts[H2Liquid->GetName()]->SetForceSolid(true);
    // Deuteron Gas
    G4Material *D2Gas = new G4Material("D2Gas", density = fTargetDensityRatio * 0.47 / 760.0 * 273.15 / 19.5 * 0.1796 * mg / cm3, ncomponents = 1, kStateGas, 19.5 * kelvin, fTargetDensityRatio * 0.47 / 760.0 * atmosphere);
    D2Gas->AddElement(D, natoms = 2);
    fVisAtts[D2Gas->GetName()] = new G4VisAttributes(G4Colour::Cyan());

    // Copper C101
    G4Material *Copper = new G4Material("Copper", density = fExtDensityRatio * 8.92 * g / cm3, ncomponents = 1);
    Copper->AddElement(Cu, natoms = 1);
    fVisAtts[Copper->GetName()] = new G4VisAttributes(G4Colour::Brown());
    G4Material *Copper0d2 = new G4Material("Copper0.2", Copper->GetDensity() * 0.2, Copper);
    fVisAtts[Copper0d2->GetName()] = new G4VisAttributes(G4Colour::Brown());
    G4Material *Copper0d75 = new G4Material("Copper0.75", Copper->GetDensity() * 0.75, Copper);
    fVisAtts[Copper0d75->GetName()] = new G4VisAttributes(G4Colour::Brown());
    G4Material *Copper0d8 = new G4Material("Copper0.8", Copper->GetDensity() * 0.8, Copper);
    fVisAtts[Copper0d8->GetName()] = new G4VisAttributes(G4Colour::Brown());

    // Kapton
    G4Material *Kapton = new G4Material("Kapton", density = fExtDensityRatio * 1.42 * g / cm3, ncomponents = 4);
    Kapton->AddElement(H, fractionmass = 0.0273);
    Kapton->AddElement(C, fractionmass = 0.7213);
    Kapton->AddElement(N, fractionmass = 0.0765);
    Kapton->AddElement(O, fractionmass = 0.1749);
    fVisAtts[Kapton->GetName()] = new G4VisAttributes(G4Colour::Brown());
    G4Material *Kapton0d2 = new G4Material("Kapton0.2", Kapton->GetDensity() * 0.2, Kapton);
    fVisAtts[Kapton0d2->GetName()] = new G4VisAttributes(G4Colour::Brown());
    G4Material *Kapton0d8 = new G4Material("Kapton0.8", Kapton->GetDensity() * 0.8, Kapton);
    fVisAtts[Kapton0d8->GetName()] = new G4VisAttributes(G4Colour::Brown());

    // Silicon
    G4Material *Silicon = new G4Material("Silicon", density = 2.329 * g / cm3, ncomponents = 1);
    Silicon->AddElement(Si, natoms = 1);
    fVisAtts[Silicon->GetName()] = new G4VisAttributes(G4Colour::Green());

    // Aluminum
    G4Material *Aluminum = new G4Material("Aluminum", density = fExtDensityRatio * 2.700 * g / cm3, ncomponents = 1);
    Aluminum->AddElement(Al, natoms = 1);
    fVisAtts[Aluminum->GetName()] = new G4VisAttributes(G4Colour::Grey());
    fVisAtts[Aluminum->GetName()]->SetForceSolid(true);
    
    //Tantalum
    G4Material *Tantalum = new G4Material("Tantalum", density = fExtDensityRatio * 16.65 * g / cm3, ncomponents = 1);
    Tantalum->AddElement(Ta, natoms = 1);
    fVisAtts[Tantalum->GetName()] = new G4VisAttributes(G4Colour::Red());
    
    // Tedlar
    G4Material *Tedlar = new G4Material("Tedlar", density = 1.545 * g / cm3, ncomponents = 3);
    Tedlar->AddElement(H, natoms = 3);
    Tedlar->AddElement(C, natoms = 2);
    Tedlar->AddElement(F, natoms = 1);
    fVisAtts[Tedlar->GetName()] = new G4VisAttributes(G4Colour::Grey());

    // Stainless Steel
    G4Material *SSteel = new G4Material("SSteel", density = fExtDensityRatio * 7.9 * g / cm3, ncomponents = 9);
    SSteel->AddElement(C, fractionmass = 0.0007);
    SSteel->AddElement(Si, fractionmass = 0.01);
    SSteel->AddElement(Mn, fractionmass = 0.02);
    SSteel->AddElement(Ni, fractionmass = 0.09);
    SSteel->AddElement(P, fractionmass = 0.00045);
    SSteel->AddElement(S, fractionmass = 0.00015);
    SSteel->AddElement(Cr, fractionmass = 0.18);
    SSteel->AddElement(N, fractionmass = 0.0011);
    SSteel->AddElement(Fe, fractionmass = 0.6976);
    fVisAtts[SSteel->GetName()] = new G4VisAttributes(G4Colour::Grey());
    fVisAtts[SSteel->GetName()]->SetForceSolid(true);

    // Nickel
    G4Material *Nickel = new G4Material("Nickel", density = fExtDensityRatio * 8.908 * g / cm3, ncomponents = 1);
    Nickel->AddElement(Ni, natoms = 1);
    fVisAtts[Nickel->GetName()] = new G4VisAttributes(G4Colour::Black());

    // GEM Frame G10
    G4Material *NemaG10 = new G4Material("NemaG10", density = fExtDensityRatio * 1.700 * g / cm3, ncomponents = 4);
    NemaG10->AddElement(Si, natoms = 1);
    NemaG10->AddElement(O, natoms = 2);
    NemaG10->AddElement(C, natoms = 3);
    NemaG10->AddElement(H, natoms = 3);
    fVisAtts[NemaG10->GetName()] = new G4VisAttributes(G4Colour::Brown());

    // Ar/CO2 Gas
    G4Material *CO2 = new G4Material("CO2", density = fExtDensityRatio * 1.842e-3 * g / cm3, ncomponents = 2);
    CO2->AddElement(C, natoms = 1);
    CO2->AddElement(O, natoms = 2);
    G4Material *ArCO2 = new G4Material("ArCO2", density = fExtDensityRatio * 1.715e-3 * g / cm3, ncomponents = 2);
    ArCO2->AddElement(Ar, fractionmass = 0.7);
    ArCO2->AddMaterial(CO2, fractionmass = 0.3);
    fVisAtts[ArCO2->GetName()] = new G4VisAttributes(G4Colour::Yellow());

    // He Gas
    G4Material *HeGas = new G4Material("HeGas", density = fExtDensityRatio * 0.1786e-3 * g / cm3, ncomponents = 1);
    HeGas->AddElement(He, natoms = 1);
    fVisAtts[HeGas->GetName()] = new G4VisAttributes(G4Colour::Cyan());

    // Scintillator EJ204
    G4Material *EJ204 = new G4Material("EJ204", density = fExtDensityRatio * 1.032 * g / cm3, ncomponents = 2);
    EJ204->AddElement(H, natoms = 521);
    EJ204->AddElement(C, natoms = 474);
    fVisAtts[EJ204->GetName()] = new G4VisAttributes(G4Colour::Green());

    // Rohacell 31 IG
    G4Material *Rohacell = new G4Material("Rohacell", density = fExtDensityRatio * 0.023 * g / cm3, ncomponents = 3);
    Rohacell->AddElement(C, natoms = 5);
    Rohacell->AddElement(H, natoms = 8);
    Rohacell->AddElement(O, natoms = 2);
    fVisAtts[Rohacell->GetName()] = new G4VisAttributes(G4Colour::Grey());

    // Tungsten
    G4Material *Tungsten = new G4Material("Tungsten", density = 19.25 * g / cm3, ncomponents = 1);
    Tungsten->AddElement(W, natoms = 1);
    fVisAtts[Tungsten->GetName()] = new G4VisAttributes(G4Colour::Black());
    
    //Lead
    G4Material *Lead = new G4Material("Lead", density = 11.35 * g / cm3, ncomponents = 1);
    Lead->AddElement(Pb, natoms = 1);
    fVisAtts[Lead->GetName()] = new G4VisAttributes(G4Colour::Red());

    // Polyester (3M VM-2000 reflector)
    G4Material *Polyester = new G4Material("Polyester", density = 1.37 * g / cm3, ncomponents = 3);
    Polyester->AddElement(C, natoms = 10);
    Polyester->AddElement(H, natoms = 8);
    Polyester->AddElement(O, natoms = 4);
    fVisAtts[Polyester->GetName()] = new G4VisAttributes(G4VisAttributes::Invisible);

    // Brass
    G4Material *Brass = new G4Material("Brass", density = 8.53 * g / cm3, ncomponents = 2);
    Brass->AddElement(Cu, fractionmass = 0.7);
    Brass->AddElement(Zn, fractionmass = 0.3);
    fVisAtts[Brass->GetName()] = new G4VisAttributes(G4Color::Brown());

    // PbWO4 Crystal
    G4Material *PbWO4 = new G4Material("PbWO4", density = 8.280 * g / cm3, ncomponents = 3);
    PbWO4->AddElement(Pb, natoms = 1);
    PbWO4->AddElement(W, natoms = 1);
    PbWO4->AddElement(O, natoms = 4);
    fVisAtts[PbWO4->GetName()] = new G4VisAttributes(G4Colour::Blue());

    // Silica
    G4Material *SiO2 = new G4Material("SiO2", density = 2.200 * g / cm3, ncomponents = 2);
    SiO2->AddElement(Si, natoms = 1);
    SiO2->AddElement(O, natoms = 2);
    fVisAtts[SiO2->GetName()] = new G4VisAttributes(G4Colour::Green());

    // Lead Glass
    G4Material *PbO = new G4Material("PbO", density = 9.530 * g / cm3, ncomponents = 2);
    PbO->AddElement(Pb, natoms = 1);
    PbO->AddElement(O, natoms = 1);

    G4Material *K2O = new G4Material("K2O", density = 2.320 * g / cm3, ncomponents = 2);
    K2O->AddElement(K, natoms = 2);
    K2O->AddElement(O, natoms = 1);

    G4Material *Na2O = new G4Material("Na2O", density = 2.270 * g / cm3, ncomponents = 2);
    Na2O->AddElement(Na, natoms = 2);
    Na2O->AddElement(O, natoms = 1);

    G4Material *As2O3 = new G4Material("As2O3", density = 3.738 * g / cm3, ncomponents = 2);
    As2O3->AddElement(As, natoms = 2);
    As2O3->AddElement(O, natoms = 3);

    G4Material *PbGlass = new G4Material("PbGlass", density = 3.86 * g / cm3, ncomponents = 5);
    PbGlass->AddMaterial(PbO, fractionmass = 0.5080);
    PbGlass->AddMaterial(SiO2, fractionmass = 0.4170);
    PbGlass->AddMaterial(K2O, fractionmass = 0.0422);
    PbGlass->AddMaterial(Na2O, fractionmass = 0.0278);
    PbGlass->AddMaterial(As2O3, fractionmass = 0.0050);
    fVisAtts[PbGlass->GetName()] = new G4VisAttributes(G4Colour::Blue());

    // Virtual Detector Material
    G4Material *VirtualDetM = new G4Material("VirtualDetM", density = universe_mean_density, ncomponents = 1, kStateGas, 0.1 * kelvin, 1.0e-19 * pascal);
    VirtualDetM->AddElement(H, fractionmass = 1.0);
    fVisAtts[VirtualDetM->GetName()] = new G4VisAttributes(G4Colour::Cyan());
    
    //Viton
    G4Material *Viton = new G4Material("Viton", density = 2.5 * g / cm3, ncomponents = 3); // HFP, VF2
    Viton->AddElement(C, natoms = 5);
    Viton->AddElement(H, natoms = 2);
    Viton->AddElement(F, natoms = 8);
    fVisAtts[Viton->GetName()] = new G4VisAttributes(G4Colour::Red());

    // Print out material table
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefinePRadVolumes()
{
    G4Material *DefaultM = G4Material::GetMaterial("Galaxy");
    G4Material *UCollimatorM = G4Material::GetMaterial("Nickel");
    G4Material *VacuumTubeM = G4Material::GetMaterial("SSteel");
    G4Material *VirtualDetM = G4Material::GetMaterial("VirtualDetM");
    //G4Material *HeBagM = G4Material::GetMaterial("HeGas");

    G4double TargetR = fTargetR;
    G4double TargetHalfL = fTargetHalfL;
    G4Material *TargetM = G4Material::GetMaterial("H2Gas");
    if (fTargetMat == "LH2"){
        TargetM = G4Material::GetMaterial("H2Liquid");
        G4cout<<"using "<<fTargetHalfL*2.<<" mm long liquid hydrogen target, with radius "<<fTargetR<<" mm"<<G4endl;
    }
    else if (fTargetMat == "Ta") {
        TargetM = G4Material::GetMaterial("Tantalum");
        G4cout<<"using "<<fTargetHalfL*2.<<" mm long Tantalum target, with radius "<<fTargetR<<" mm"<<G4endl;
    }
    else{
        G4cout<<"using "<<fTargetHalfL*2.<<" mm long gas hydrogen target, with radius "<<fTargetR<<" mm"<<G4endl;
    }
    // World
    G4VSolid *solidWorld = new G4Box("WorldS", fWorldSizeXY, fWorldSizeXY, fWorldSizeZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, DefaultM, "WorldLV");
    G4VPhysicalVolume *physiWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0);

    // Target
    // Target Container
    G4VSolid *solidTargetCon = new G4Box("TargetContainerS", TargetR + 10.0 * mm, TargetR + 10.0 * mm, TargetHalfL + 1.0 * mm);
    G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, DefaultM, "TargetContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fTargetCenter), logicTargetCon, "Target Container", logicWorld, false, 0);

    // Target material
    G4VSolid *solidTarget = new G4Tubs("TargetS", 0, TargetR, TargetHalfL, 0, twopi);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetM, "TargetLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);

    if (fTargetMat == "LH2"){
        G4Material *TargetCellM = G4Material::GetMaterial("Aluminum");
        G4Tubs *CellTube1 = new G4Tubs("CellTube1", 0, TargetR + 2 * mm, TargetHalfL + 0.125 * mm, 0, twopi);
        G4Tubs *CellTube2 = new G4Tubs("CellTube2", 0, TargetR + 0.0001 * mm, TargetHalfL + 0.0001 * mm, 0, twopi);
        G4SubtractionSolid * solidCell = new G4SubtractionSolid("TargetCellS", CellTube1, CellTube2);
        G4LogicalVolume* logicCell = new G4LogicalVolume(solidCell, TargetCellM, "TargetCellLV");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);
    }
    else if (fTargetMat == "Ta"){
        //do nothing here
    }
    else{
        G4Material *TargetCellM = G4Material::GetMaterial("Copper");
        G4Material *TargetWindowM = G4Material::GetMaterial("Kapton");
        // Target cell
        G4double CellXY = 3.5 * cm;
        G4Box *CellBox = new G4Box("CellBox", CellXY, CellXY, TargetHalfL);
        G4Tubs *CellTube = new G4Tubs("CellTube", 0, TargetR, TargetHalfL + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidCell = new G4SubtractionSolid("TargetCellS", CellBox, CellTube);
        G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TargetCellM, "TargetCellLV");
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

        // Target cell windows
        G4double CellApertureR = 2.0 * mm;
        G4double CellWinThickness = 7.5 * um;
        G4Box *CellWinBox = new G4Box("CellWinBox", CellXY, CellXY, CellWinThickness / 2.0);
        G4Tubs *CellWinTube = new G4Tubs("CellWinTube", 0, CellApertureR, CellWinThickness + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidCellWin = new G4SubtractionSolid("TargetWindowS", CellWinBox, CellWinTube);
        G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, TargetWindowM, "TargetWindowLV");
        new G4PVPlacement(0, G4ThreeVector(0, 0, -TargetHalfL - CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 0);
        new G4PVPlacement(0, G4ThreeVector(0, 0, +TargetHalfL + CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 1);
    }
    // Upstream collimator
    // Dimension from PRad beam line drawing (search PRad in JLab drawing database)
    G4double UCollimatorHalfL = 11.8 * 2.54 / 2.0 * cm;
    G4double UCollimatorOR = 3.9 * 2.54 / 2.0 * cm;
    G4double UCollimatorIR = 12.7 * mm / 2.0;
    G4double UCollimatorCenter = fTargetCenter - 2.03 * m + UCollimatorHalfL;
    G4VSolid *solidUCollimator = new G4Tubs("UCollimatorS", UCollimatorIR, UCollimatorOR, UCollimatorHalfL, 0, twopi);
    G4LogicalVolume *logicUCollimator = new G4LogicalVolume(solidUCollimator, UCollimatorM, "UCollimatorLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, UCollimatorCenter), logicUCollimator, "Upstream Collimator", logicWorld, false, 0);

    // Upstream beam pipe
    // Dimension from PRad target drawings and PRad beam line drawings
    G4double UBeamPipeIR = 1.87 * 2.54 / 2.0 * cm;
    G4double UBeamPipeOR = 2.0 * 2.54 / 2.0 * cm;
    G4double UBeamPipeHalfL = 0.79 * m;
    G4double UBeamPipeOffset = 0.105 * m;
    G4double UBeamPipeCenter = fTargetCenter - UBeamPipeOffset - UBeamPipeHalfL - 4 * cm; // subtract 4cm for target length
    G4VSolid *solidUBeamPipe = new G4Tubs("UBeamPipeS", UBeamPipeIR, UBeamPipeOR, UBeamPipeHalfL, 0, 2 * pi);
    G4LogicalVolume *logicUBeamPipe = new G4LogicalVolume(solidUBeamPipe, VacuumTubeM, "UBeamPipeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, UBeamPipeCenter), logicUBeamPipe, "Upstream Beam Pipe", logicWorld, false, 0);

    AddVaccumBox(logicWorld);

    // Center of two GEM should be at -3000.0 + 89.0 + (5226.16 + 5186.45) / 2 + 4.6525 = 2299.9575 mm // (5226.16 + 5186.45) / 2 from Weizhi
    //fGEMCenter[0] = 229.99575 * cm;
    //AddGEM(logicWorld, 0, false);

    AddGEM(logicWorld, 0, false);
    AddGEM(logicWorld, 1, false);
    
    // He bag (Only He gas for now)
    //G4Box *HeBagBox = new G4Box("HeBagBox", 1.0 * m, 1.0 * m, (fGEMCenter[1] - fGEMCenter[0] - 5.65 * cm) / 2.0 - 1*cm);
    //G4Tubs *HeBagTube = new G4Tubs("HeBagTube", 0, 22.0 * mm, (fGEMCenter[1] - fGEMCenter[0] - 5.65 * cm + 1.0 * mm) / 2.0 - 1*cm, 0, twopi);
    //G4SubtractionSolid *solidHeBag = new G4SubtractionSolid("HeBagS", HeBagBox, HeBagTube);
    //G4LogicalVolume *logicHeBag = new G4LogicalVolume(solidHeBag, HeBagM, "HeBagLV");
    //new G4PVPlacement(0, G4ThreeVector(0, 0, (fGEMCenter[0] + fGEMCenter[1]) / 2.0), logicHeBag, "He Bag", logicWorld, false, 0);
    
    // The crystal surface should be at -3000.0 + 89.0 + 5646.15 = 2735.15 mm // 5646.15 from Weizhi
    //fCrystalSurf = 273.515 * cm; // Surface of the PWO
    
    AddHyCal(logicWorld);
    
    // Virtual Detector
    /*G4double VirtualDetR = 50.0 * cm;
    G4double VirtualDetZ = 0.1 * mm;
    G4VSolid *solidVirtualDet = new G4Tubs("VirtualDetS", 0. , VirtualDetR, VirtualDetZ / 2.0, 0, twopi);
    G4LogicalVolume *logicVirtualDet = new G4LogicalVolume(solidVirtualDet, VirtualDetM, "VirtualDetLV");
    //new G4PVPlacement(0, G4ThreeVector(0, 0, fGEMCenter[0] - 3.5*cm), logicVirtualDet, "Virtual Detector", logicWorld, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, fTargetCenter + TargetHalfL + 10 * mm), logicVirtualDet, "Virtual Detector", logicWorld, false, 0);
    */
    G4LogicalVolumeStore *pLogicalVolume = G4LogicalVolumeStore::GetInstance();

    for (unsigned long i = 0; i < pLogicalVolume->size(); i++)
        (*pLogicalVolume)[i]->SetVisAttributes(fVisAtts[(*pLogicalVolume)[i]->GetMaterial()->GetName()]);

    // Always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefinePRadSDs()
{
    if (fTargetSDOn) {
        CheckScatteringSD *TargetSD = new CheckScatteringSD("TargetSD", "TG");
        G4SDManager::GetSDMpointer()->AddNewDetector(TargetSD);
        SetSensitiveDetector("TargetLV", TargetSD);
    }

    if (fGEMSDOn) {
        TrackingDetectorSD *GEMSD = new TrackingDetectorSD("GEMSD", "GEM");
        G4SDManager::GetSDMpointer()->AddNewDetector(GEMSD);
        SetSensitiveDetector("GEM0CathodeLV", GEMSD);
        SetSensitiveDetector("GEM1CathodeLV", GEMSD);
    }

    if (fHyCalSDOn) {
        CalorimeterSD *HyCalSD = new CalorimeterSD("HyCalSD", "HC", "database/pwo_attenuation.dat");
        HyCalSD->SetAttenuationLG(fAttenuationLG);
        G4SDManager::GetSDMpointer()->AddNewDetector(HyCalSD);

        for (int i = 0; i < 1152; i++)
            SetSensitiveDetector(Form("PbWO4Absorber%04dLV", i), HyCalSD);

        for (int i = 0; i < 576; i++)
            SetSensitiveDetector(Form("PbGlassAbsorber%04dLV", i), HyCalSD);
    }

    if (fVirtualSDOn) {
        StandardDetectorSD *VirtualSD = new StandardDetectorSD("VirtualSD", "VD");
        G4SDManager::GetSDMpointer()->AddNewDetector(VirtualSD);
        SetSensitiveDetector("VirtualDetLV", VirtualSD);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineDRadVolumes()
{
    G4Material *DefaultM = G4Material::GetMaterial("Galaxy");
    G4Material *TargetM = G4Material::GetMaterial("D2Gas");

    if (fTargetMat == "H2Gas" || fTargetMat == "hydrogen")
        TargetM = G4Material::GetMaterial("H2Gas");

    G4Material *TargetCellM = G4Material::GetMaterial("Kapton");
    G4Material *TargetWindowM = G4Material::GetMaterial("Kapton");
    G4Material *RecoilDetectorM = G4Material::GetMaterial("Silicon");
    G4Material *RecoilDetCoverM = G4Material::GetMaterial("SiO2");
    G4Material *HeBagM = G4Material::GetMaterial("HeGas");
    G4Material *ScintillatorPlaneM = G4Material::GetMaterial("EJ204");

    // World
    G4VSolid *solidWorld = new G4Box("WorldS", fWorldSizeXY, fWorldSizeXY, fWorldSizeZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, DefaultM, "WorldLV");
    G4VPhysicalVolume *physiWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0);

    // Target
    G4double CellXY = 15.0 * cm;

    // Target Container
    G4VSolid *solidTargetCon = new G4Box("TargetContainerS", CellXY + 0.1 * cm, CellXY + 0.1 * cm, fTargetHalfL + 0.1 * cm);
    G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, DefaultM, "TargetContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fTargetCenter), logicTargetCon, "Target Container", logicWorld, false, 0);

    // Target material
    G4VSolid *solidTarget = new G4Tubs("TargetS", 0, fTargetR, fTargetHalfL, 0, twopi);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetM, "TargetLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);

    // Target cell
    G4Box *CellBox = new G4Box("CellBox", CellXY, CellXY, fTargetHalfL);
    G4Tubs *CellTube = new G4Tubs("CellTube", 0, fTargetR, fTargetHalfL + 1.0 * mm, 0, twopi);
    G4SubtractionSolid *solidCell = new G4SubtractionSolid("TargetCellS", CellBox, CellTube);
    G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TargetCellM, "TargetCellLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

    // Target cell windows
    G4double CellApertureR = 2.0 * mm;
    G4double CellWinThickness = 7.5 * um;
    G4Box *CellWinBox = new G4Box("CellWinBox", CellXY, CellXY, CellWinThickness / 2.0);
    G4Tubs *CellWinTube = new G4Tubs("CellWinTube", 0, CellApertureR, CellWinThickness + 1.0 * mm, 0, twopi);
    G4SubtractionSolid *solidCellWin = new G4SubtractionSolid("TargetWindowS", CellWinBox, CellWinTube);
    G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, TargetWindowM, "TargetWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -fTargetHalfL - CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, +fTargetHalfL + CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 1);

    // Recoil detector
    G4double RecoilDetCenter = fRecoilDetCenter - fTargetCenter;
    G4double CoverThickness = 0.5 * um;
    G4double RecoilDetAng = twopi / fRecoilDetNSeg;
    G4double RecoilDetL2OR = fRecoilDetR * cos(RecoilDetAng / 2.0);
    G4double RecoilDetL2IRC = RecoilDetL2OR - fRecoilDetL2Thickness + CoverThickness;
    G4double RecoilDetL2IR = RecoilDetL2OR - fRecoilDetL2Thickness;
    G4double RecoilDetL1OR = RecoilDetL2IR;
    G4double RecoilDetL1IRC = RecoilDetL1OR - fRecoilDetL1Thickness + CoverThickness;
    G4double RecoilDetL1IR = RecoilDetL1OR - fRecoilDetL1Thickness;

    G4double rInnerRDL2[] = {RecoilDetL2IRC, RecoilDetL2IRC};
    G4double rOuterRDL2[] = {RecoilDetL2OR, RecoilDetL2OR};
    G4double zPlaneRDL2[] = {-fRecoilDetHalfL, fRecoilDetHalfL};
    G4VSolid *solidRecoilDet2 = new G4Polyhedra("RecoilDet2S", 0, twopi, fRecoilDetNSeg, 2, zPlaneRDL2, rInnerRDL2, rOuterRDL2);

    G4double rInnerRDL2Cover[] = {RecoilDetL2IR, RecoilDetL2IR};
    G4double rOuterRDL2Cover[] = {RecoilDetL2IRC, RecoilDetL2IRC};
    G4double zPlaneRDL2Cover[] = {-fRecoilDetHalfL, fRecoilDetHalfL};
    G4VSolid *solidRecoilDet2Cover = new G4Polyhedra("RecoilDet2CoverS", 0, twopi, fRecoilDetNSeg, 2, zPlaneRDL2Cover, rInnerRDL2Cover, rOuterRDL2Cover);

    G4double rInnerRDL1[] = {RecoilDetL1IRC, RecoilDetL1IRC};
    G4double rOuterRDL1[] = {RecoilDetL1OR, RecoilDetL1OR};
    G4double zPlaneRDL1[] = {-fRecoilDetHalfL, fRecoilDetHalfL};
    G4VSolid *solidRecoilDet1 = new G4Polyhedra("RecoilDet1S", 0, twopi, fRecoilDetNSeg, 2, zPlaneRDL1, rInnerRDL1, rOuterRDL1);

    G4double rInnerRDL1Cover[] = {RecoilDetL1IR, RecoilDetL1IR};
    G4double rOuterRDL1Cover[] = {RecoilDetL1IRC, RecoilDetL1IRC};
    G4double zPlaneRDL1Cover[] = {-fRecoilDetHalfL, fRecoilDetHalfL};
    G4VSolid *solidRecoilDet1Cover = new G4Polyhedra("RecoilDet1CoverS", 0, twopi, fRecoilDetNSeg, 2, zPlaneRDL1Cover, rInnerRDL1Cover, rOuterRDL1Cover);

    G4LogicalVolume *logicRecoilDet1 = new G4LogicalVolume(solidRecoilDet1, RecoilDetectorM, "RecoilDet1LV");
    G4LogicalVolume *logicRecoilDet1Cover = new G4LogicalVolume(solidRecoilDet1Cover, RecoilDetCoverM, "RecoilDet1CoverLV");
    G4LogicalVolume *logicRecoilDet2 = new G4LogicalVolume(solidRecoilDet2, RecoilDetectorM, "RecoilDet2LV");
    G4LogicalVolume *logicRecoilDet2Cover = new G4LogicalVolume(solidRecoilDet2Cover, RecoilDetCoverM, "RecoilDet2CoverLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, RecoilDetCenter), logicRecoilDet1Cover, "Recoil Detector 1 Cover", logicTarget, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, RecoilDetCenter), logicRecoilDet1, "Recoil Detector 1", logicTarget, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, RecoilDetCenter), logicRecoilDet2Cover, "Recoil Detector 2 Cover", logicTarget, false, 1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, RecoilDetCenter), logicRecoilDet2, "Recoil Detector 2", logicTarget, false, 1);

    AddVaccumBox(logicWorld);

    AddGEM(logicWorld, 0, true);
    AddGEM(logicWorld, 1, false);

    // He bag (Only He gas for now)
    G4Box *HeBagBox = new G4Box("HeBagBox", 1.0 * m, 1.0 * m, (fGEMCenter[1] - fGEMCenter[0] - 5.65 * cm) / 2.0);
    G4Tubs *HeBagTube = new G4Tubs("HeBagTube", 0, 22.0 * mm, (fGEMCenter[1] - fGEMCenter[0] - 5.65 * cm + 1.0 * mm) / 2.0, 0, twopi);
    G4SubtractionSolid *solidHeBag = new G4SubtractionSolid("HeBagS", HeBagBox, HeBagTube);
    G4LogicalVolume *logicHeBag = new G4LogicalVolume(solidHeBag, HeBagM, "HeBagLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, (fGEMCenter[0] + fGEMCenter[1]) / 2.0), logicHeBag, "He Bag", logicWorld, false, 0);

    // Scintillator plane
    G4double SciPlaneThickness = 5.0 * mm;
    G4double SciPlaneHalfX = 75.0 * cm;
    G4double SciPlaneHalfY = 75.0 * cm;
    G4Box *SciPlaneBox = new G4Box("ScintillatorPlaneBox", SciPlaneHalfX, SciPlaneHalfY, SciPlaneThickness / 2.0);
    G4Tubs *SciPlaneTube = new G4Tubs("ScintillatorPlaneTube", 0, 22.0 * mm, (SciPlaneThickness + 1.0 * mm) / 2.0, 0, twopi);
    G4SubtractionSolid *solidSciPlane = new G4SubtractionSolid("ScintillatorPlaneS", SciPlaneBox, SciPlaneTube);
    G4LogicalVolume *logicSciPlane = new G4LogicalVolume(solidSciPlane, ScintillatorPlaneM, "ScintillatorPlaneLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fSciPlaneCenter), logicSciPlane, "Scintillator Plane", logicWorld, false, 0);

    AddHyCal(logicWorld);

    G4LogicalVolumeStore *pLogicalVolume = G4LogicalVolumeStore::GetInstance();

    for (unsigned long i = 0; i < pLogicalVolume->size(); i++)
        (*pLogicalVolume)[i]->SetVisAttributes(fVisAtts[(*pLogicalVolume)[i]->GetMaterial()->GetName()]);

    // Always return the physical World
    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineDRadSDs()
{
    if (fRecoilDetSDOn) {
        TrackingDetectorSD *RecoilDetSD = new TrackingDetectorSD("RecoilDetectorSD", "RD");
        G4SDManager::GetSDMpointer()->AddNewDetector(RecoilDetSD);
        SetSensitiveDetector("RecoilDet1LV", RecoilDetSD);
        SetSensitiveDetector("RecoilDet2LV", RecoilDetSD);
    }

    if (fGEMSDOn) {
        TrackingDetectorSD *GEMSD = new TrackingDetectorSD("GEMSD", "GEM");
        G4SDManager::GetSDMpointer()->AddNewDetector(GEMSD);
        SetSensitiveDetector("GEM0CathodeLV", GEMSD);
        SetSensitiveDetector("GEM1CathodeLV", GEMSD);
    }

    if (fSciPlaneSDOn) {
        StandardDetectorSD *SciPlaneSD = new StandardDetectorSD("ScintillatorPlaneSD", "SP");
        G4SDManager::GetSDMpointer()->AddNewDetector(SciPlaneSD);
        SetSensitiveDetector("ScintillatorPlaneLV", SciPlaneSD);
    }

    if (fHyCalSDOn) {
        CalorimeterSD *HyCalSD = new CalorimeterSD("HyCalSD", "HC", "database/pwo_attenuation.dat");
        HyCalSD->SetAttenuationLG(fAttenuationLG);
        G4SDManager::GetSDMpointer()->AddNewDetector(HyCalSD);

        for (int i = 0; i < 1152; i++)
            SetSensitiveDetector(Form("PbWO4Absorber%04dLV", i), HyCalSD);

        for (int i = 0; i < 576; i++)
            SetSensitiveDetector(Form("PbGlassAbsorber%04dLV", i), HyCalSD);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineTestVolumes()
{
    G4Material *DefaultM = G4Material::GetMaterial("Galaxy");
    G4Material *TargetM = G4Material::GetMaterial("H2Gas");
    G4Material *TargetCellM = G4Material::GetMaterial("Copper");
    G4Material *TargetWindowM = G4Material::GetMaterial("Kapton");
    G4Material *VirtualDetM = G4Material::GetMaterial("VirtualDetM");

    // World
    G4double WorldSizeXY = 10.0 * cm;
    G4double WorldSizeZ = 100.0 * cm;
    G4VSolid *solidWorld = new G4Box("WorldS", WorldSizeXY, WorldSizeXY, WorldSizeZ);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, DefaultM, "WorldLV");
    G4VPhysicalVolume *physiWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0);

    // Target
    // Target Container
    G4VSolid *solidTargetCon = new G4Box("TargetContainerS", 3.5 * cm, 3.5 * cm, 2.1 * cm);
    G4LogicalVolume *logicTargetCon = new G4LogicalVolume(solidTargetCon, DefaultM, "TargetContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTargetCon, "Target Container", logicWorld, false, 0);

    // Target material
    G4double TargetR = 25.0 * mm;
    G4double TargetHalfL = 20.0 * mm;
    G4VSolid *solidTarget = new G4Tubs("TargetS", 0, TargetR, TargetHalfL, 0, twopi);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, TargetM, "TargetLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTarget, "Target Material", logicTargetCon, false, 0);
    //logicTarget->SetUserLimits(new G4UserLimits(0.5 * mm));

    // Target cell
    G4double CellXY = 3.5 * cm;
    G4Box *CellBox = new G4Box("CellBox", CellXY, CellXY, TargetHalfL);
    G4Tubs *CellTube = new G4Tubs("CellTube", 0, TargetR, TargetHalfL + 1.0 * mm, 0, twopi);
    G4SubtractionSolid *solidCell = new G4SubtractionSolid("TargetCellS", CellBox, CellTube);
    G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell, TargetCellM, "TargetCellLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicCell, "Target Cell", logicTargetCon, false, 0);

    // Target cell windows
    G4double CellApertureR = 2.0 * mm;
    G4double CellWinThickness = 7.5 * um;
    G4Box *CellWinBox = new G4Box("CellWinBox", CellXY, CellXY, CellWinThickness / 2.0);
    G4Tubs *CellWinTube = new G4Tubs("CellWinTube", 0, CellApertureR, CellWinThickness + 1.0 * mm, 0, twopi);
    G4SubtractionSolid *solidCellWin = new G4SubtractionSolid("TargetWindowS", CellWinBox, CellWinTube);
    G4LogicalVolume *logicCellWin = new G4LogicalVolume(solidCellWin, TargetWindowM, "TargetWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -TargetHalfL - CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, +TargetHalfL + CellWinThickness / 2.0), logicCellWin, "Target Window", logicTargetCon, false, 1);

    // Virtual Detector
    G4double VirtualDetZ = 0.1 * mm;
    G4double VirtualDetL = 99.0 * cm;
    G4double VirtualDetIR = (VirtualDetL - 20.0 * mm) * tan(0.5 / 180.0 * pi);
    G4double VirtualDetOR = (VirtualDetL + 20.0 * mm) * tan(10.0 / 180.0 * pi);
    G4VSolid *solidVirtualDet = new G4Tubs("VirtualDetS", VirtualDetIR, VirtualDetOR, VirtualDetZ, 0, twopi);
    G4LogicalVolume *logicVirtualDet = new G4LogicalVolume(solidVirtualDet, VirtualDetM, "VirtualDetLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, VirtualDetL), logicVirtualDet, "Virtual Detector", logicWorld, false, 0);

    G4LogicalVolumeStore *pLogicalVolume = G4LogicalVolumeStore::GetInstance();

    for (unsigned long i = 0; i < pLogicalVolume->size(); i++)
        (*pLogicalVolume)[i]->SetVisAttributes(fVisAtts[(*pLogicalVolume)[i]->GetMaterial()->GetName()]);

    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineTestSDs()
{
    StepRecordSD *TargetSD = new StepRecordSD("TargetSD", "TG");
    G4SDManager::GetSDMpointer()->AddNewDetector(TargetSD);
    SetSensitiveDetector("TargetLV", TargetSD);

    StandardDetectorSD *VirtualSD = new StandardDetectorSD("VirtualSD", "VD");
    G4SDManager::GetSDMpointer()->AddNewDetector(VirtualSD);
    SetSensitiveDetector("VirtualDetLV", VirtualSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddVaccumBox(G4LogicalVolume *mother)
{
    G4Material *ChamberM = G4Material::GetMaterial("Aluminum");
    G4Material *ChamberWindowM = G4Material::GetMaterial("Kapton");
    G4Material *VacuumBoxM = G4Material::GetMaterial("Aluminum");
    G4Material *VacuumTubeM = G4Material::GetMaterial("Aluminum");
    G4Material *FlangeM1 = G4Material::GetMaterial("Aluminum");
    G4Material *FlangeM2 = G4Material::GetMaterial("SSteel");
    G4Material *FlangeM3 = G4Material::GetMaterial("Viton");

    // Target chamber
    // For now, only built the downstream chamber with window
    // The downstream chamber window should locate at -3000.0 + 89.0 + 74.0  = -2837.0 mm
    // The length of the downstream chamber is 381.7 mm
    // The total length of the downstream chamber and the tube in total is 710.0 mm
    // Here the downstream chamber and the tube are built together to be the new downstream chamber.
    // So the center of this geometry should be at -2837.0 + 710.0 / 2 = -2482.0 mm
    G4double shiftCenter = fTargetCenter + (1600.00)*mm;
    fDownChamberCenter = shiftCenter + 74.0 * mm + 71.0 * cm / 2.0;
    G4double DownChamberHalfL = 71.0 / 2.0 * cm;
    G4double DownChamberUR = 8.00 * cm;

    //downstream beam pipe
    G4double DownstreamBeamPipeHalfZ = (fDownChamberCenter - DownChamberHalfL - (fTargetCenter + 20*cm))/2.;
    G4double DownstreamBeamPipeCenter = fTargetCenter + 20*cm + DownstreamBeamPipeHalfZ;
    G4VSolid *solidDownStreamBeamPipe = new G4Tubs("DownStreamBeamPipe", 17.30*cm, 17.78*cm, DownstreamBeamPipeHalfZ *0.999, 0, twopi);
    G4LogicalVolume *logicDownStreamBeamPipe = new G4LogicalVolume(solidDownStreamBeamPipe, ChamberM, "DownStreamBeamPipeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, DownstreamBeamPipeCenter), logicDownStreamBeamPipe, "Down Stream Beam Pipe", mother, false, 0);

    // Downstream chamber
    //G4double rInnerDC[] = {7.56 * cm, 7.56 * cm, 7.56 * cm, 7.56 * cm, 17.30 * cm, 17.30 * cm};
    //G4double rOuterDC[] = {8.00 * cm, 8.00 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm};
    //G4double zPlaneDC[] = {0,         32.83 * cm, 32.83 * cm, 35.37 * cm, 35.37 * cm, 71.00 * cm};
    //G4VSolid *solidDownChamber = new G4Polycone("DownstreamChamberS", 0, twopi, 6, zPlaneDC, rInnerDC, rOuterDC);
    //G4LogicalVolume *logicDownChamber = new G4LogicalVolume(solidDownChamber, ChamberM, "DownstreamChamberLV");
    //new G4PVPlacement(0, G4ThreeVector(0, 0, fDownChamberCenter - DownChamberHalfL), logicDownChamber, "Downstream Chamber", mother, false, 0);

    // Downstream chamber
    //G4double rInnerDC[] = {7.56 * cm, 7.56 * cm, 7.56 * cm, 7.56 * cm, 17.30 * cm, 17.30 * cm};
    //G4double rOuterDC[] = {8.00 * cm, 8.00 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm};
    G4double rInnerDC[] = {17.30 * cm, 17.30 * cm,  17.30 * cm,  17.30 * cm,  17.30 * cm, 17.30 * cm};
    G4double rOuterDC[] = {17.78 * cm, 17.78 * cm,  17.78 * cm, 17.78 * cm, 17.78 * cm, 17.78 * cm};
    G4double zPlaneDC[] = {0,         32.83 * cm, 32.83 * cm, 35.37 * cm, 35.37 * cm, 71.00 * cm};
    G4VSolid *solidDownChamber = new G4Polycone("DownstreamChamberS", 0, twopi, 6, zPlaneDC, rInnerDC, rOuterDC);
    G4LogicalVolume *logicDownChamber = new G4LogicalVolume(solidDownChamber, ChamberM, "DownstreamChamberLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fDownChamberCenter - DownChamberHalfL), logicDownChamber, "Downstream Chamber", mother, false, 0);

    // Downstream chamber window
    G4double DownChamberApertureR = 22.8 * mm;
    G4double DownChamberWinThickness = 7.5 * um;
    G4Tubs *solidDownChamberWin = new G4Tubs("DownstreamChamberWindowS", DownChamberApertureR, DownChamberUR, DownChamberWinThickness / 2.0, 0, twopi);
    G4LogicalVolume *logicDownChamberWin = new G4LogicalVolume(solidDownChamberWin, ChamberWindowM, "DownstreamChamberWindowLV");
    //new G4PVPlacement(0, G4ThreeVector(0, 0, fDownChamberCenter - DownChamberHalfL - DownChamberWinThickness / 2.0), logicDownChamberWin, "Downstream Chamber Window", mother, false, 0);

    // Vacuum box
    // The length of the vacuum box is 4251.7 mm
    // So the center of this geometry should be at -3000.0 + 89.0 + 74.0 + 710.0 + 2125.85 = -1.15 mm
    fVacBoxCenter = shiftCenter + 74.0 * mm + 71.0 * cm + 425.17 * cm / 2.0;
    G4double VacBoxHalfL = 425.17 * cm / 2.0;
    G4double VacBoxMaxR = 78.11 * cm;
    G4double rInner2[] = {17.30 * cm, 17.30 * cm, 50.17 * cm, 50.17 * cm, 78.11 * cm, 78.11 * cm};
    G4double rOuter2[] = {17.78 * cm, 17.78 * cm, 50.80 * cm, 50.80 * cm, 78.74 * cm, 78.74 * cm};
    G4double zPlane2[] = {0, 6.8 * cm, 17.6 * cm, 215.3 * cm, 229.5 * cm, 425.17 * cm};
    G4VSolid *solidVacBox = new G4Polycone("VacuumBoxS", 0, twopi, 6, zPlane2, rInner2, rOuter2);
    G4LogicalVolume *logicVacBox = new G4LogicalVolume(solidVacBox, VacuumBoxM, "VacuumBoxLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fVacBoxCenter - VacBoxHalfL), logicVacBox, "Vacuum Box", mother, false, 0);

    // Vacuum box window
    G4double VacBoxWinFlangeOffset = 3.81 * cm;
    G4double ArcDistance = 5.59 * cm;
    G4double ArcEndR = (ArcDistance * ArcDistance + VacBoxMaxR * VacBoxMaxR) / (2 * ArcDistance);
    G4double ArcEndThickness = 1.0 * mm;
    G4double VacBoxWinApertureR = 3.0 * cm;
    G4VSolid *solidVacBoxWin = new G4Sphere("VacuumBoxWindowS", ArcEndR - ArcEndThickness, ArcEndR, 0, twopi, pi - asin(VacBoxMaxR / ArcEndR), asin(VacBoxMaxR / ArcEndR) - asin((VacBoxWinApertureR + 0.1 * mm) / ArcEndR));
    G4LogicalVolume *logicVacBoxWin = new G4LogicalVolume(solidVacBoxWin, VacuumBoxM, "VacuumBoxWindowLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fVacBoxCenter + VacBoxHalfL + ArcEndR - ArcDistance - VacBoxWinFlangeOffset), logicVacBoxWin, "Vacuum Box Window", mother, false, 0);

    // Vacuum Tube
    G4double VacTubeOR = 1.75 * cm; //previous 1.9 cm, 1.375 in = 1.75cm
    G4double VacTubeIR = VacTubeOR - 0.0889 * cm; // 0.049 in = 0.12446 cm from Eugene, 0.035 in = 0.0889 cm for Al pipe
    G4double VacTubeL = fWorldSizeZ - 10.0 * cm - fVacBoxCenter - VacBoxHalfL + ArcDistance + VacBoxWinFlangeOffset - (34.6+4.0+1.5)*2.*mm; //The last is the lenth of flange
    G4VSolid *solidVacTube = new G4Tubs("VacuumTubeS", VacTubeIR, VacTubeOR, VacTubeL / 2.0, 0, twopi);
    G4LogicalVolume *logicVacTube = new G4LogicalVolume(solidVacTube, VacuumTubeM, "VacuumTubeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fWorldSizeZ - 10.0 * cm - VacTubeL / 2.0), logicVacTube, "Vacuum Tube", mother, false, 0);

    //Flange
    G4double FlangeCenter = fVacBoxCenter + VacBoxHalfL - ArcDistance - VacBoxWinFlangeOffset + (34.6*mm + 4.0*mm + 1.5*mm);//The last is the half lenth of flange
    //center ring
    G4double RingID = 40.0 * mm;
    G4double RingMD = RingID + 0.5 * mm; //Middle
    G4double RingOD = RingMD + 6.85 * mm;
    G4double SteelRingL = 8.0 * mm;
    G4double RingL = 3.0 * mm;
    G4VSolid* solidSteelRing = new G4Tubs("SteelRingS", RingID*0.5, RingMD*0.5, SteelRingL*0.5, 0, twopi);
    G4LogicalVolume* logicSteelRing = new G4LogicalVolume(solidSteelRing, FlangeM2, "StellRingLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter), logicSteelRing, "StellRing", mother, false, 0);
    G4VSolid* SolidRing = new G4Tubs("RingS", RingMD*0.5, RingOD*0.5, RingL*0.5, 0, twopi);
    G4LogicalVolume* logicRing = new G4LogicalVolume(SolidRing, FlangeM3, "RingLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter), logicRing, "Ring", mother, false, 0);

    //flange tube
    G4double FlangeOD = 55.0 * mm;
    G4double FlangeID = 34.8 * mm;
    G4double FlangeL = 4.0 * mm;
    G4VSolid* solidFlangeTub = new G4Tubs("FlangeTube", FlangeID*0.5, FlangeOD*0.5, FlangeL*0.5, 0, twopi);

    G4double BlankOD = 41.5 * mm;
    G4double BlankL = 2.55 * mm;
    G4VSolid* solidBlank = new G4Tubs("BlankS", 0., BlankOD*0.5, BlankL*0.5, 0, twopi);

    G4SubtractionSolid* solidFlange1 = new G4SubtractionSolid("FlangeS1", solidFlangeTub, solidBlank, 0, G4ThreeVector(0, 0, -0.5*FlangeL+0.5*BlankL-0.05*mm));
    G4LogicalVolume* logicFlange1 = new G4LogicalVolume(solidFlange1, FlangeM1, "Flange1LV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter+RingL*0.5+FlangeL*0.5+0.1*mm), logicFlange1, "Flange1", mother, false, 0);

    G4SubtractionSolid* solidFlange2 = new G4SubtractionSolid("FlangeS2", solidFlangeTub, solidBlank, 0, G4ThreeVector(0, 0, 0.5*FlangeL-0.5*BlankL+0.05*mm));
    G4LogicalVolume* logicFlange2 = new G4LogicalVolume(solidFlange2, FlangeM1, "Flange2LV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter-RingL*0.5-FlangeL*0.5-0.1*mm), logicFlange2, "Flange2", mother, false, 0);

    G4double FTubeOD = 38.1 * mm;
    G4double FTubeID = 34.8 * mm;
    G4double FTubeL = 34.6 * mm;
    G4VSolid* solidFTube = new G4Tubs("FTubeS", FTubeID*0.5, FTubeOD*0.5, FTubeL*0.5, 0, twopi);
    G4LogicalVolume* logicFTube = new G4LogicalVolume(solidFTube, FlangeM1, "FTubeLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter+RingL*0.5+FlangeL+FTubeL*0.5+0.2*mm), logicFTube, "FTube1", mother, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter-RingL*0.5-FlangeL-FTubeL*0.5-0.2*mm), logicFTube, "FTube2", mother, false, 0);

    //Hing Clamp on flange
    G4double HingTubeID = 49.0 * mm;
    G4double HingTubeOD = 67.82 * mm;
    G4double HingTubeL = 15.20 * mm;
    G4double HingBlankTubeOD = FlangeOD + 1.0 * mm;
    G4double HingBlankTubeL = RingL + FlangeL*2. + 1.*mm;
    G4VSolid* solidHingTube = new G4Tubs("HingTubeS", HingTubeID*0.5, HingTubeOD*0.5, HingTubeL*0.5, 0, twopi);
    G4VSolid* solidHingBlankTube = new G4Tubs("HingBlankTubeS", 0., HingBlankTubeOD*0.5, HingBlankTubeL*0.5, 0, twopi);
    G4SubtractionSolid* solidHingClamp = new G4SubtractionSolid("HingClampS", solidHingTube, solidHingBlankTube, 0, G4ThreeVector(0, 0, 0));
    G4LogicalVolume* logicHingClamp = new G4LogicalVolume(solidHingClamp, FlangeM1, "HingClampLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, FlangeCenter), logicHingClamp, "HingClamp", mother, false, 0);

    G4VSolid* solidClampBox = new G4Box("ClampBoxS", 24.*mm*0.5, 7.*mm*0.5, 15.2*mm*0.5);
    G4LogicalVolume* logicClampBox = new G4LogicalVolume(solidClampBox, FlangeM1, "ClampBoXLV");
    G4VSolid* solidClampBox2 = new G4Box("ClampBoxS2", 40.*mm*0.5, 7.*mm*0.5, 15.2*mm*0.5);
    G4LogicalVolume* logicClampBox2 = new G4LogicalVolume(solidClampBox2, FlangeM1, "ClampBoX2LV");
    new G4PVPlacement(0, G4ThreeVector(0, 7.*mm*0.5+HingTubeOD*0.5+1.0*mm, FlangeCenter), logicClampBox, "ClampBox1", mother, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, -7.*mm*0.5-HingTubeOD*0.5-1.0*mm, FlangeCenter), logicClampBox2, "ClampBox2", mother, false, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddGEM(G4LogicalVolume *mother, int layerid, bool culess)
{
    G4Material *DefaultM = G4Material::GetMaterial("Galaxy");
    G4Material *GEMFrameM = G4Material::GetMaterial("NemaG10");
    G4Material *GEMGasM = G4Material::GetMaterial("ArCO2");
    G4Material *GEMFoilM = G4Material::GetMaterial("Kapton");
    G4Material *GEMFoil0d2M = G4Material::GetMaterial("Kapton0.2");
    G4Material *GEMFoil0d8M = G4Material::GetMaterial("Kapton0.8");
    G4Material *GEMCuM = G4Material::GetMaterial("Copper");
    G4Material *GEMCu0d2M = G4Material::GetMaterial("Copper0.2");
    G4Material *GEMCu0d75M = G4Material::GetMaterial("Copper0.75");
    G4Material *GEMCu0d8M = G4Material::GetMaterial("Copper0.8");
    G4Material *GEMGlueM = G4Material::GetMaterial("Kapton"); // TODO: Add actual Glue material

    // GEM
    G4double GEMCenter = fGEMCenter[layerid];
    G4double GEMGap = 3.971 * cm; // Gap between two GEM // 3.971 from Weizhi
    G4double GEMHalfX = 55.04 * cm / 2.0;
    G4double GEMHalfY = 122.88 * cm / 2.0;
    G4double GEMHalfT = (15.0 * mm + 455.0 * um) / 2.0; // 2 * 25 + 5 + 50 (win) + 6 * 5 + 3 * 50 (foil) + 5 + 5 + 50 + 50 + 60 (readout)

    if (culess) GEMHalfT = (15.0 * mm + 410.0 * um) / 2.0; // 2 * 25 + 50 (win) + 3 * 50 (foil) + 50 + 50 + 60 (readout)

    G4double GEMSpacerWh = 0.3 * mm / 2.0;
    G4double GEMSpacerWv = 0.5 * mm / 2.0;
    G4double GEMSpacerT = (2.0 - 0.1) * mm;
    G4double GEMHoleR = 2.2 * cm;
    G4double GEMCenterHalfXY = 7.4 * cm / 2.0;
    G4double GEMFrameWidth = 1.5 * cm;
    G4double GEMCenterOffset = GEMHalfX + GEMFrameWidth - GEMCenterHalfXY;

    // GEM Container
    G4Box *GEMConBox = new G4Box(Form("GEM%dConBox", layerid), 1.0 * m, 1.0 * m, (GEMGap + 2.0 * GEMHalfT + 1.0 * mm) / 2.0);
    G4Tubs *GEMConTube = new G4Tubs(Form("GEM%dConTube", layerid), 0, GEMHoleR, (GEMGap + 2.0 * GEMHalfT + 1.0 * mm) / 2.0 + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEMCon = new G4SubtractionSolid(Form("GEM%dContainerS", layerid), GEMConBox, GEMConTube);
    G4LogicalVolume *logicGEMCon = new G4LogicalVolume(solidGEMCon, DefaultM, Form("GEM%dContainerLV", layerid));
    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMCenter), logicGEMCon, Form("GEM %d Container", layerid), mother, false, 2 * layerid);

    // GEM
    G4Box *GEMBox = new G4Box(Form("GEM%dBox", layerid), GEMHalfX + GEMFrameWidth, GEMHalfY + GEMFrameWidth * 2.0, GEMHalfT);
    G4Tubs *GEMTube = new G4Tubs(Form("GEM%dTube", layerid), 0, GEMHoleR, GEMHalfT + 0.1 * mm, 0, twopi);
    G4SubtractionSolid *solidGEM = new G4SubtractionSolid(Form("GEM%dS", layerid), GEMBox, GEMTube, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEM = new G4LogicalVolume(solidGEM, DefaultM, Form("GEM%dLV", layerid));
    new G4PVPlacement(0, G4ThreeVector(GEMCenterOffset, 0, GEMGap / 2.0), logicGEM, Form("GEM %d L", layerid), logicGEMCon, false, 0);
    G4RotationMatrix rmGEM;
    rmGEM.rotateZ(180.0 * deg);
    new G4PVPlacement(G4Transform3D(rmGEM, G4ThreeVector(-GEMCenterOffset, 0, -GEMGap / 2.0)), logicGEM, Form("GEM %d R", layerid), logicGEMCon, false, 1);

    // GEM Gas
    G4Box *GEMGasBox = new G4Box(Form("GEM%dGasBox", layerid), GEMHalfX, GEMHalfY, GEMHalfT);
    G4Box *GEMSubBox = new G4Box(Form("GEM%dSubBox", layerid), GEMCenterHalfXY, GEMCenterHalfXY, GEMHalfT + 0.1 * mm);
    G4SubtractionSolid *solidGEMGas = new G4SubtractionSolid(Form("GEM%dGasS", layerid), GEMGasBox, GEMSubBox, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMGas = new G4LogicalVolume(solidGEMGas, GEMGasM, Form("GEM%dGasLV", layerid));
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMGas, Form("GEM %d Gas", layerid), logicGEM, false, 0);

    // GEM Frame
    G4Box *GEMFrameBox1 = new G4Box(Form("GEM%dFrameBox1", layerid), GEMHalfX + GEMFrameWidth, GEMHalfY + GEMFrameWidth * 2.0, GEMHalfT);
    G4Box *GEMFrameBox2 = new G4Box(Form("GEM%dFrameBox2", layerid), GEMHalfX, GEMHalfY, GEMHalfT + 0.1 * mm);
    G4SubtractionSolid *solidGEMFrame = new G4SubtractionSolid(Form("GEM%dFrameS", layerid), GEMFrameBox1, GEMFrameBox2);
    G4LogicalVolume *logicGEMFrame = new G4LogicalVolume(solidGEMFrame, GEMFrameM, Form("GEM%dFrameLV", layerid));
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGEMFrame, Form("GEM %d Frame", layerid), logicGEM, false, 0);

    G4Box *GEMPipeBox = new G4Box(Form("GEM%dPipeBox", layerid), GEMCenterHalfXY - GEMFrameWidth / 2.0, GEMCenterHalfXY, GEMHalfT);
    G4SubtractionSolid *solidGEMPipe = new G4SubtractionSolid(Form("GEM%dPipeS", layerid), GEMPipeBox, GEMTube, 0, G4ThreeVector(-GEMFrameWidth / 2.0, 0, 0));
    G4LogicalVolume *logicGEMPipe = new G4LogicalVolume(solidGEMPipe, GEMFrameM, Form("GEM%dPipeLV", layerid));
    new G4PVPlacement(0, G4ThreeVector(-GEMCenterOffset + GEMFrameWidth / 2.0, 0, 0), logicGEMPipe, Form("GEM %d Pipe", layerid), logicGEM, false, 0);

    // GEM Spacer
    G4double GEMSpacerOffset = 91.3 * mm;
    G4Box *GEMHSpacerBox1 = new G4Box(Form("GEM%dHSpacerBox1", layerid), GEMHalfX, GEMSpacerWh, GEMSpacerT / 2.0);
    G4Box *GEMHSpacerBox2 = new G4Box(Form("GEM%dHSpacerBox2", layerid), GEMHalfX - GEMCenterHalfXY + GEMFrameWidth / 2.0, GEMSpacerWh, GEMSpacerT / 2.0);
    G4Box *GEMVSpacerBox = new G4Box(Form("GEM%dVSpacerBox", layerid), GEMSpacerWv, GEMHalfY, GEMSpacerT / 2.0);
    G4UnionSolid *GEMSpacerPiece1 = new G4UnionSolid(Form("GEM%dSpacerPiece1", layerid), GEMVSpacerBox, GEMHSpacerBox1, 0, G4ThreeVector(GEMSpacerOffset, 204.0 * mm, 0));
    G4UnionSolid *GEMSpacerPiece2 = new G4UnionSolid(Form("GEM%dSpacerPiece2", layerid), GEMSpacerPiece1, GEMHSpacerBox1, 0, G4ThreeVector(GEMSpacerOffset, 409.3 * mm, 0));
    G4UnionSolid *GEMSpacerPiece3 = new G4UnionSolid(Form("GEM%dSpacerPiece3", layerid), GEMSpacerPiece2, GEMHSpacerBox2, 0, G4ThreeVector(GEMSpacerOffset + GEMCenterHalfXY - GEMFrameWidth / 2.0, 0, 0));
    G4UnionSolid *GEMSpacerPiece4 = new G4UnionSolid(Form("GEM%dSpacerPiece4", layerid), GEMSpacerPiece3, GEMHSpacerBox1, 0, G4ThreeVector(GEMSpacerOffset, -204.0 * mm, 0));
    G4UnionSolid *GEMSpacerPiece5 = new G4UnionSolid(Form("GEM%dSpacerPiece5", layerid), GEMSpacerPiece4, GEMHSpacerBox1, 0, G4ThreeVector(GEMSpacerOffset, -409.3 * mm, 0));
    G4UnionSolid *solidGEMSpacer = new G4UnionSolid(Form("GEM%dSpacerS", layerid), GEMSpacerPiece5, GEMVSpacerBox, 0, G4ThreeVector(GEMSpacerOffset * 2.0, 0, 0));
    G4LogicalVolume *logicGEMSpacer = new G4LogicalVolume(solidGEMSpacer, GEMFrameM, Form("GEM%dSpacerLV", layerid));

    // GEM Foil
    G4double GEMWinT = 25.0 * um;
    G4double GEMFoilT = 50.0 * um;
    G4double GEMCuT = 5.0 * um;
    G4double GEMGlueT = 60.0 * um;

    G4Box *GEMWinBox = new G4Box(Form("GEM%dWinBox", layerid), GEMHalfX, GEMHalfY, GEMWinT / 2.0);
    G4SubtractionSolid *solidGEMWin = new G4SubtractionSolid(Form("GEM%dWinS", layerid), GEMWinBox, GEMSubBox, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMWin = new G4LogicalVolume(solidGEMWin, GEMFoilM, Form("GEM%dWinLV", layerid));

    G4Box *GEMFoilBox = new G4Box(Form("GEM%dFoilBox", layerid), GEMHalfX, GEMHalfY, GEMFoilT / 2.0);
    G4SubtractionSolid *solidGEMFoil = new G4SubtractionSolid(Form("GEM%dFoilS", layerid), GEMFoilBox, GEMSubBox, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMFoil = new G4LogicalVolume(solidGEMFoil, GEMFoil0d8M, Form("GEM%dFoilLV", layerid));
    G4LogicalVolume *logicGEMFoil80 = new G4LogicalVolume(solidGEMFoil, GEMFoil0d2M, Form("GEM%dFoil80LV", layerid));
    G4LogicalVolume *logicGEMFoil350 = new G4LogicalVolume(solidGEMFoil, GEMFoilM, Form("GEM%dFoil350LV", layerid));
    G4LogicalVolume *logicGEMCathode = new G4LogicalVolume(solidGEMFoil, GEMFoilM, Form("GEM%dCathodeLV", layerid));

    G4Box *GEMCuBox = new G4Box(Form("GEM%dCuBox", layerid), GEMHalfX, GEMHalfY, GEMCuT / 2.0);
    G4SubtractionSolid *solidGEMCu = new G4SubtractionSolid(Form("GEM%dCuS", layerid), GEMCuBox, GEMSubBox, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMCu = new G4LogicalVolume(solidGEMCu, GEMCu0d8M, Form("GEM%dCuLV", layerid));
    G4LogicalVolume *logicGEMCu80 = new G4LogicalVolume(solidGEMCu, GEMCu0d2M, Form("GEM%dCu80LV", layerid));
    G4LogicalVolume *logicGEMCu350 = new G4LogicalVolume(solidGEMCu, GEMCu0d75M, Form("GEM%dCu350LV", layerid));
    G4LogicalVolume *logicGEMCathodeCu = new G4LogicalVolume(solidGEMCu, GEMCuM, Form("GEM%dCathodeCuLV", layerid));

    G4Box *GEMGlueBox = new G4Box(Form("GEM%dGlueBox", layerid), GEMHalfX, GEMHalfY, GEMGlueT / 2.0);
    G4SubtractionSolid *solidGEMGlue = new G4SubtractionSolid(Form("GEM%dGlueS", layerid), GEMGlueBox, GEMSubBox, 0, G4ThreeVector(-GEMCenterOffset, 0, 0));
    G4LogicalVolume *logicGEMGlue = new G4LogicalVolume(solidGEMGlue, GEMGlueM, Form("GEM%dGlueLV", layerid));

    G4double zoff = -GEMHalfT;

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMWinT / 2.0), logicGEMWin, Form("GEM %d Window", layerid), logicGEMGas, false, 0);
    zoff += GEMWinT;

    new G4PVPlacement(0, G4ThreeVector(-GEMSpacerOffset, 0, zoff + GEMSpacerT / 2.0), logicGEMSpacer, Form("GEM %d Spacer", layerid), logicGEMGas, false, 0);
    zoff += 3.0 * mm;

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMFoilT / 2.0), logicGEMCathode, Form("GEM %d Cathode", layerid), logicGEMGas, false, 0);
    zoff += GEMFoilT;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCathodeCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 0);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(-GEMSpacerOffset, 0, zoff + GEMSpacerT / 2.0), logicGEMSpacer, Form("GEM %d Spacer", layerid), logicGEMGas, false, 1);
    zoff += 3.0 * mm;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 1);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMFoilT / 2.0), logicGEMFoil, Form("GEM %d Foil", layerid), logicGEMGas, false, 0);
    zoff += GEMFoilT;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 2);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(-GEMSpacerOffset, 0, zoff + GEMSpacerT / 2.0), logicGEMSpacer, Form("GEM %d Spacer", layerid), logicGEMGas, false, 2);
    zoff += 2.0 * mm;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 3);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMFoilT / 2.0), logicGEMFoil, Form("GEM %d Foil", layerid), logicGEMGas, false, 1);
    zoff += GEMFoilT;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 4);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(-GEMSpacerOffset, 0, zoff + GEMSpacerT / 2.0), logicGEMSpacer, Form("GEM %d Spacer", layerid), logicGEMGas, false, 3);
    zoff += 2.0 * mm;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 5);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMFoilT / 2.0), logicGEMFoil, Form("GEM %d Foil", layerid), logicGEMGas, false, 2);
    zoff += GEMFoilT;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu, Form("GEM %d Copper", layerid), logicGEMGas, false, 6);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(-GEMSpacerOffset, 0, zoff + GEMSpacerT / 2.0 + 0.01 * mm), logicGEMSpacer, Form("GEM %d Spacer", layerid), logicGEMGas, false, 4);
    zoff += 2.0 * mm;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu80, Form("GEM %d Copper", layerid), logicGEMGas, false, 7);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMFoilT / 2.0), logicGEMFoil80, Form("GEM %d Foil", layerid), logicGEMGas, false, 3);
    zoff += GEMFoilT;

    if (!culess) {
        new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMCuT / 2.0), logicGEMCu350, Form("GEM %d Copper", layerid), logicGEMGas, false, 8);
        zoff += GEMCuT;
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMFoilT / 2.0), logicGEMFoil350, Form("GEM %d Foil", layerid), logicGEMGas, false, 4);
    zoff += GEMFoilT;
    new G4PVPlacement(0, G4ThreeVector(0, 0, zoff + GEMGlueT / 2.0), logicGEMGlue, Form("GEM %d Glue", layerid), logicGEMGas, false, 0);
    zoff += GEMGlueT;

    new G4PVPlacement(0, G4ThreeVector(0, 0, GEMHalfT - GEMWinT / 2.0), logicGEMWin, Form("GEM %d Window", layerid), logicGEMGas, false, 1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::AddHyCal(G4LogicalVolume *mother)
{
    G4Material *HyCalConM = G4Material::GetMaterial("Air");
    G4Material *HyCalBoxM = G4Material::GetMaterial("Rohacell");
    G4Material *HyCalBoxWinM = G4Material::GetMaterial("Tedlar");
    G4Material *CollimatorM = G4Material::GetMaterial("Tungsten");
    G4Material *PbWO4ModuleM = G4Material::GetMaterial("PbWO4");
    G4Material *PbGlassModuleM = G4Material::GetMaterial("PbGlass");
    G4Material *TedlarTapeM = G4Material::GetMaterial("Tedlar");
    G4Material *ReflectorM = G4Material::GetMaterial("Polyester");
    G4Material *PlateM = G4Material::GetMaterial("Brass");
    G4Material *LeadM = G4Material::GetMaterial("Lead");
    G4Material *VirtualDetM = G4Material::GetMaterial("VirtualDetM");
    
    std::ifstream dimension_file;
    dimension_file.open("database/hycal_module_shuffled.dat");
    G4double pwo[1152][4], lg[576][4];

    for (int i = 0; i < 576; i++)
        dimension_file >> lg[i][0] >> lg[i][1] >> lg[i][2] >> lg[i][3];

    for (int i = 0; i < 1152; i++)
        dimension_file >> pwo[i][0] >> pwo[i][1] >> pwo[i][2] >> pwo[i][3];

    dimension_file.close();

    G4double MaxStep = 1.0 * mm;

    // HyCal
    G4double CrystalL = 18.0 * cm;
    G4double PbGlassL = 45.0 * cm;
    G4double CrystalDiffL = 9.73 * cm; // according to last survey (april 2017)
    G4double CrystalCenter = fCrystalSurf + CrystalL / 2.0;
    G4double PbGlassCenter = fCrystalSurf - CrystalDiffL + PbGlassL / 2.0;

    // HyCal container
    G4double HyCalBoxCenter = PbGlassCenter - 9.0 * cm + 30.0 * cm; // Check
    G4Box *HyCalConNoHole = new G4Box("HyCalConNoHole", 80.0 * cm, 80.0 * cm, 65.0 * cm);
    G4Tubs *HyCalConHole = new G4Tubs("HyCalConHole", 0, 31.75 * mm, 66.0 * cm, 0, twopi);
    G4SubtractionSolid *solidHyCalCon = new G4SubtractionSolid("HyCalConS", HyCalConNoHole, HyCalConHole);
    G4LogicalVolume *logicHyCalCon = new G4LogicalVolume(solidHyCalCon, HyCalConM, "HyCalConLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, HyCalBoxCenter), logicHyCalCon, "HyCal Container", mother, false, 0);
    logicHyCalCon->SetUserLimits(new G4UserLimits(MaxStep));

    // HyCal box
    G4Box *HyCalBoxOuter = new G4Box("HyCalBoxOuter", 72.54 * cm, 72.54 * cm, 62.54 * cm);
    G4Box *HyCalBoxInner = new G4Box("HyCalBoxInner", 70.0 * cm, 70.0 * cm, 60.0 * cm);
    G4SubtractionSolid *HyCalBoxNoHole = new G4SubtractionSolid("HyCalBoxNoHole", HyCalBoxOuter, HyCalBoxInner);
    G4Tubs *HyCalBoxHole = new G4Tubs("HyCalBoxHole", 0, 3.175 * cm, 65.0 * cm, 0, twopi);
    G4SubtractionSolid *solidHyCalBox = new G4SubtractionSolid("HyCalBoxS", HyCalBoxNoHole, HyCalBoxHole);
    G4LogicalVolume *logicHyCalBox = new G4LogicalVolume(solidHyCalBox, HyCalBoxM, "HyCalBoxLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicHyCalBox, "HyCal Box", logicHyCalCon, false, 0);

    // HyCal box window
    G4VSolid *solidHyCalBoxWin = new G4Tubs("HyCalBoxWinS", 1.90 * cm, 5.08 * cm, 19.0 * um, 0, twopi);
    G4LogicalVolume *logicHyCalBoxWin = new G4LogicalVolume(solidHyCalBoxWin, HyCalBoxWinM, "HyCalBoxWinLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -62.54 * cm - 19.0 * um), logicHyCalBoxWin, "HyCal Box Window", logicHyCalCon, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -60.00 * cm + 19.0 * um), logicHyCalBoxWin, "HyCal Box Window", logicHyCalCon, false, 1);

    // Lead glass module container
    G4double PlateThickness = 25.4 * 0.08 * mm; // 0.08 in
    G4double PbGlassPlateHoleR = 25.4 * 0.5 * mm; // 0.5 in
    G4double ReflectorT = 25.0 * um;
    G4double StripThickness = 25.4 * 0.001 * mm; // 0.001 in

    G4Box *PbGlassConBox = new G4Box("PbGlassConBox", 58.199 * cm, 58.165 * cm, PbGlassL / 2.0 + PlateThickness);
    G4Box *PbGlassConHole = new G4Box("PbGlassConHole", 35.309 * cm, 35.275 * cm, PbGlassL / 2.0 + PlateThickness + 1.0 * mm);
    G4SubtractionSolid *solidPbGlassCon = new G4SubtractionSolid("PbGlassModuleContainerS", PbGlassConBox, PbGlassConHole);
    G4LogicalVolume *logicPbGlassCon = new G4LogicalVolume(solidPbGlassCon, HyCalConM, "PbGlassModuleContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, PbGlassCenter - HyCalBoxCenter), logicPbGlassCon, "PbGlass Module Container", logicHyCalCon, false, 0);

    for (int i = 0; i < 576; i++) {
        G4double PbGlassX = lg[i][0] - 0.15 * mm;
        G4double PbGlassY = lg[i][1] - 0.15 * mm;

        // Lead glass module
        G4VSolid *solidPbGlassModule = new G4Box(Form("PbGlassModule%04dS", i), lg[i][0] * mm / 2.0, lg[i][1] * mm / 2.0, PbGlassL / 2.0 + PlateThickness);
        G4LogicalVolume *logicPbGlassModule = new G4LogicalVolume(solidPbGlassModule, HyCalConM, Form("PbGlassModule%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(lg[i][2] * mm, lg[i][3] * mm, 0), logicPbGlassModule, "PbGlass Module", logicPbGlassCon, false, i);

        G4VSolid *solidPbGlassReflector = new G4Box(Form("PbGlassReflector%04dS", i), PbGlassX / 2.0 + ReflectorT, PbGlassY / 2.0 + ReflectorT, PbGlassL / 2.0);
        G4LogicalVolume *logicPbGlassReflector = new G4LogicalVolume(solidPbGlassReflector, ReflectorM, Form("PbGlassReflector%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicPbGlassReflector, "Reflector", logicPbGlassModule, false, 0);

        G4VSolid *solidPbGlassAbsorber = new G4Box(Form("PbGlassAbsorber%04dS", i), PbGlassX / 2.0, PbGlassY / 2.0, PbGlassL / 2.0);
        G4LogicalVolume *logicPbGlassAbsorber = new G4LogicalVolume(solidPbGlassAbsorber, PbGlassModuleM, Form("PbGlassAbsorber%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicPbGlassAbsorber, "PbGlass Absorber", logicPbGlassReflector, false, 0);

        G4Box *PbGlassPlateBox = new G4Box("PbGlassPlateBox", PbGlassX / 2.0, PbGlassY / 2.0, PlateThickness / 2.0);
        G4Tubs *PbGlassPlateHole = new G4Tubs("PbGlassPlateHole", 0, PbGlassPlateHoleR, PlateThickness / 2.0 + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidPbGlassPlate = new G4SubtractionSolid(Form("PbGlassPlate%04dS", i), PbGlassPlateBox, PbGlassPlateHole);
        G4LogicalVolume *logicPbGlassPlate = new G4LogicalVolume(solidPbGlassPlate, PlateM, Form("PbGlassPlate%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, - PbGlassL / 2.0 - PlateThickness / 2.0), logicPbGlassPlate, "Brass Plate", logicPbGlassModule, false, 0);


        if ((lg[i][2] > -353.09 && lg[i][3] > 352.75) || (lg[i][2] < 353.09 && lg[i][3] < -352.75)) {
            G4double StripWidth = PbGlassX;

            G4VSolid *solidPbGlassStrip = new G4Box(Form("PbGlassStrip%04dS", i), StripWidth / 2.0, StripThickness / 2.0,  PbGlassL / 2.0);
            G4LogicalVolume *logicPbGlassStrip = new G4LogicalVolume(solidPbGlassStrip, PlateM, Form("PbGlassStrip%04dLV", i));
            new G4PVPlacement(0, G4ThreeVector(0, (lg[i][1] * mm - StripThickness) / 2.0, 0), logicPbGlassStrip, "Brass Strip", logicPbGlassModule, false, 0);
            new G4PVPlacement(0, G4ThreeVector(0, (-lg[i][1] * mm + StripThickness) / 2.0, 0), logicPbGlassStrip, "Brass Strip", logicPbGlassModule, false, 1);
        } else {
            G4double StripWidth = PbGlassY;

            G4VSolid *solidPbGlassStrip = new G4Box(Form("PbGlassStrip%04dS", i), StripThickness / 2.0, StripWidth / 2.0, PbGlassL / 2.0);
            G4LogicalVolume *logicPbGlassStrip = new G4LogicalVolume(solidPbGlassStrip, PlateM, Form("PbGlassStrip%04dLV", i));
            new G4PVPlacement(0, G4ThreeVector((lg[i][0] * mm - StripThickness) / 2.0, 0, 0), logicPbGlassStrip, "Brass Strip", logicPbGlassModule, false, 0);
            new G4PVPlacement(0, G4ThreeVector((-lg[i][0] * mm + StripThickness) / 2.0, 0, 0), logicPbGlassStrip, "Brass Strip", logicPbGlassModule, false, 1);
        }
    }

    // PbWO4 module container
    G4double CrystalPlateHoleR = 25.4 * 0.25 * mm; // 0.25 in
    G4double TedlarTapeT = 38.1 * um;
    ReflectorT = 63.0 * um;

    G4Box *PbWO4ConBox = new G4Box("PbWO4ConBox", 35.309 * cm, 35.275 * cm, CrystalL / 2.0 + PlateThickness);
    G4Box *PbWO4ConHole = new G4Box("PbWO4ConHole", 2.0 * cm, 2.0 * cm, CrystalL / 2.0 + PlateThickness + 1.0 * mm);
    G4SubtractionSolid *solidPbWO4Con = new G4SubtractionSolid("PbWO4ModuleContainerS", PbWO4ConBox, PbWO4ConHole);
    G4LogicalVolume *logicPbWO4Con = new G4LogicalVolume(solidPbWO4Con, HyCalConM, "PbWO4ModuleContainerLV");
    new G4PVPlacement(0, G4ThreeVector(0, 0, CrystalCenter - HyCalBoxCenter), logicPbWO4Con, "PbWO4 Module Container", logicHyCalCon, false, 576);

    for (int i = 0; i < 1152; i++) {
        G4double CrystalX = pwo[i][0] - 0.27 * mm;
        G4double CrystalY = pwo[i][1] - 0.25 * mm;
        G4double StripWidth = 25.4 * 0.5 * mm; // 0.5 in

        // PbWO4 module
        G4VSolid *solidPbWO4Module = new G4Box(Form("PbWO4Module%04dS", i), pwo[i][0] * mm / 2.0, pwo[i][1] * mm / 2.0, CrystalL / 2.0 + PlateThickness);
        G4LogicalVolume *logicPbWO4Module = new G4LogicalVolume(solidPbWO4Module, HyCalConM, Form("PbWO4Module%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(pwo[i][2] * mm, pwo[i][3] * mm, 0), logicPbWO4Module, "PbGlass Module", logicPbWO4Con, false, i);

        G4VSolid *solidPbWO4TedlarTape = new G4Box(Form("PbWO4TedlarTape%04dS", i), CrystalX / 2.0 + TedlarTapeT + ReflectorT, CrystalY / 2.0 + TedlarTapeT + ReflectorT, CrystalL / 2.0);
        G4LogicalVolume *logicPbWO4TedlarTape = new G4LogicalVolume(solidPbWO4TedlarTape, TedlarTapeM, Form("PbWO4TedlarTape%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicPbWO4TedlarTape, "Tedlar Tape", logicPbWO4Module, false, 0);

        G4VSolid *solidPbWO4Reflector = new G4Box(Form("PbWO4Reflector%04dS", i), CrystalX / 2.0 + ReflectorT, CrystalY / 2.0 + ReflectorT, CrystalL / 2.0);
        G4LogicalVolume *logicPbWO4Reflector = new G4LogicalVolume(solidPbWO4Reflector, ReflectorM, Form("PbWO4Reflector%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicPbWO4Reflector, "Reflector", logicPbWO4TedlarTape, false, 0);

        G4VSolid *solidPbWO4Absorber = new G4Box(Form("PbWO4Absorber%04dS", i), CrystalX / 2.0, CrystalY / 2.0, CrystalL / 2.0);
        G4LogicalVolume *logicPbWO4Absorber = new G4LogicalVolume(solidPbWO4Absorber, PbWO4ModuleM, Form("PbWO4Absorber%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicPbWO4Absorber, "PbWO4 Absorber", logicPbWO4Reflector, false, 0);

        G4Box *PbWO4PlateBox = new G4Box("PbWO4PlateBox", CrystalX / 2.0, CrystalY / 2.0, PlateThickness / 2.0);
        G4Tubs *PbWO4PlateHole = new G4Tubs("PbWO4PlateHole", 0, CrystalPlateHoleR, PlateThickness / 2.0 + 1.0 * mm, 0, twopi);
        G4SubtractionSolid *solidPbWO4Plate = new G4SubtractionSolid(Form("PbWO4Plate%04dS", i), PbWO4PlateBox, PbWO4PlateHole);
        G4LogicalVolume *logicPbWO4Plate = new G4LogicalVolume(solidPbWO4Plate, PlateM, Form("PbWO4Plate%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector(0, 0, -(CrystalL + PlateThickness) / 2.0), logicPbWO4Plate, "Brass Plate", logicPbWO4Module, false, 0);

        G4VSolid *solidPbWO4Strip = new G4Box(Form("PbWO4Strip%04dS", i), StripThickness / 2.0, StripWidth / 2.0, CrystalL / 2.0);
        G4LogicalVolume *logicPbWO4Strip = new G4LogicalVolume(solidPbWO4Strip, PlateM, Form("PbWO4Strip%04dLV", i));
        new G4PVPlacement(0, G4ThreeVector((pwo[i][0] * mm - StripThickness) / 2.0, 0, 0), logicPbWO4Strip, "Brass Strip", logicPbWO4Module, false, 0);
        new G4PVPlacement(0, G4ThreeVector((-pwo[i][0] * mm + StripThickness) / 2.0, 0, 0), logicPbWO4Strip, "Brass Strip", logicPbWO4Module, false, 1);
    }

   // Collimator
    G4Box *CollimatorBox = new G4Box("CollimatorBox", 4.07 * cm, 4.07 * cm, 3.02 * cm);
    G4Tubs *CollimatorTube = new G4Tubs("CollimatorTube", 0, 1.95 * cm, 3.1 * cm, 0, twopi);
    G4SubtractionSolid *solidCollimator = new G4SubtractionSolid("CollimatorS", CollimatorBox, CollimatorTube);
    G4LogicalVolume *logicCollimator = new G4LogicalVolume(solidCollimator, CollimatorM, "CollimatorLV");
    G4RotationMatrix rmColl;
    rmColl.rotateZ(-0.0 * deg);
    new G4PVPlacement(G4Transform3D(rmColl, G4ThreeVector(0, 0, fCrystalSurf - PlateThickness - 3.1 * cm - HyCalBoxCenter-0.05*cm)), logicCollimator, "Collimator", logicHyCalCon, false, 0);
    
    //Virtual Det 
    G4double VirtualDetZ = 0.01 * mm;
    //G4VSolid *solidVirtualDet1 = new G4Box("VirtualDetS1", 600*mm, 600*mm, VirtualDetZ);
    G4VSolid *solidVirtualDet1 = new G4Box("VirtualDetS1", 350*mm, 350*mm, VirtualDetZ);
    G4VSolid *solidVirtualDet2 = new G4Tubs("VirtualDet2", 0, 1.95 * cm, 3.1 * cm, 0, twopi);
    
    G4SubtractionSolid *solidVirtualDetSub = new G4SubtractionSolid("solidVirtualDetSub", solidVirtualDet1, solidVirtualDet2);
    G4LogicalVolume *logicVirtualDet = new G4LogicalVolume(solidVirtualDetSub, VirtualDetM, "VirtualDetLV");

    
    //new G4PVPlacement(0, G4ThreeVector(0., 0, CrystalCenter - HyCalBoxCenter - CrystalL/2. - 11*cm ), logicVirtualDet, "Virtual Detector", logicHyCalCon, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0., 0, CrystalCenter - HyCalBoxCenter - CrystalL/2. - 0.3*cm ), logicVirtualDet, "Virtual Detector", logicHyCalCon, false, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
