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
// TrackingMessenger.cc
// Developer : Chao Gu
// History:
//   Apr 2017, C. Gu, Original version.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingMessenger.hh"

#include "TrackingAction.hh"

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingMessenger::TrackingMessenger(TrackingAction *act) : G4UImessenger(), Action(act)
{
    TrackingDir = new G4UIdirectory("/pradsim/tracking/");
    TrackingDir->SetGuidance("Tracking control");

    NoSecondaryCmd = new G4UIcmdWithABool("/pradsim/tracking/nosecondary", this);
    NoSecondaryCmd->SetGuidance("Turn on/off secondary particles");
    NoSecondaryCmd->SetParameterName("nosecondary", false);

    SaveTrackInfoCmd = new G4UIcmdWithABool("/pradsim/tracking/savetrackinfo", this);
    SaveTrackInfoCmd->SetGuidance("Turn on/off saving tracking infomation");
    SaveTrackInfoCmd->SetParameterName("savetrackinfo", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingMessenger::~TrackingMessenger()
{
    delete NoSecondaryCmd;
    delete SaveTrackInfoCmd;
    delete TrackingDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if (command == NoSecondaryCmd)
        Action->SetNoSecondary(NoSecondaryCmd->GetNewBoolValue(newValue));

    if (command == SaveTrackInfoCmd)
        Action->SetSaveTrackInfo(SaveTrackInfoCmd->GetNewBoolValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
