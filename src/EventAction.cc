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
// EventAction.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "EventMessenger.hh"
#include "GlobalVars.hh"
#include "RootTree.hh"
#include "StandardHit.hh"
#include "TrackingAction.hh"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TTree.h"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4UserEventAction.hh"

#include "G4String.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(G4String conf) : G4UserEventAction(), fEventID(0), fPrintModulo(1000), fOnlyRecordHits(false)
{
    //fCollName = "HCColl";
    fCollName = "VDColl";

    if (conf == "test")
        fCollName = "VDColl";

    Register(gRootTree->GetTree());

    eventMessenger = new EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
    delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *evt)
{
    fEventID = evt->GetEventID();

    if (fEventID % fPrintModulo == 0)
        G4cout << "\n---> Begin of event: " << fEventID << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *evt)
{
    if (fOnlyRecordHits) {
        G4HCofThisEvent *HCE = evt->GetHCofThisEvent();

        G4int nHC = HCE->GetNumberOfCollections();
	bool WriteFlag = false;
        for (G4int i = 0; i < nHC; i++) {
            G4String ColName = HCE->GetHC(i)->GetName();
	    StandardHitsCollection *HitColl = (StandardHitsCollection *) HCE->GetHC(i);
	    G4int nHits = HitColl->entries();
	    if (nHits > 0){
                if (ColName == "HCColl" || ColName == "GEMColl")  { // Hard-coded detector name in DetectorConstruction.cc
                    WriteFlag = true;
                }
	    }
        }
        if (WriteFlag) gRootTree->FillTree();

    } else
        gRootTree->FillTree();

    TrackingAction *theTrackingAction = (TrackingAction *)G4RunManager::GetRunManager()->GetUserTrackingAction();
    theTrackingAction->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::Register(TTree *tree)
{
    tree->Branch("EventID", &fEventID, "EventID/I");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
