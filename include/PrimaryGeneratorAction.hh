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
// PrimaryGeneratorAction.hh
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Rewrite with class PrimaryGenerator.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4String.hh"
#include "G4ThreeVector.hh"

class G4VPrimaryGenerator;

class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(G4String conf);
    virtual ~PrimaryGeneratorAction();

    void GeneratePrimaries(G4Event *);

    inline void SetGunType(G4String val);
    inline void SetEventType(G4String evtype);
    inline void SetRecoilParticle(G4String val);

    inline void SetBeamEnergy(G4double val);
    inline void SetPosition(G4ThreeVector xyz);
    inline void SetTheta(G4double theta);
    inline void SetPhi(G4double phi);
    inline void SetThetaRange(G4double lo, G4double hi);

    inline void SetEnpRange(G4double lo, G4double hi);

    inline void SetEventFile(G4String path);
    inline void SetPileUpProfile(G4String profile);
    inline void SetTargetProfile(G4String profile);
    
    inline void SetPointPID(G4String val);

private:
    G4String fConfig;

    G4String fGunType;
    G4String fEventType;
    G4String fRecoilParticle;
    G4String fPointPID;

    double fE;
    double fX, fY, fZ;
    double fTheta, fPhi;
    double fThetaLo, fThetaHi;

    double fEnpLo, fEnpHi;

    G4String fEventFile;
    G4String fTargetProfile;
    G4String fPileUpProfile;

    G4VPrimaryGenerator *fPrimaryGenerator;

    PrimaryGeneratorMessenger *gunMessenger; // pointer to the messenger
};

inline void PrimaryGeneratorAction::SetGunType(G4String val)
{
    fGunType = val;
}

inline void PrimaryGeneratorAction::SetPointPID(G4String val)
{
    fPointPID = val;
}

inline void PrimaryGeneratorAction::SetEventType(G4String val)
{
    fEventType = val;
}

inline void PrimaryGeneratorAction::SetRecoilParticle(G4String val)
{
    fRecoilParticle = val;
}

inline void PrimaryGeneratorAction::SetBeamEnergy(G4double val)
{
    fE = val;
}

inline void PrimaryGeneratorAction::SetPosition(G4ThreeVector xyz)
{
    fX = xyz.x();
    fY = xyz.y();
    fZ = xyz.z();
}

inline void PrimaryGeneratorAction::SetTheta(G4double theta)
{
    fTheta = theta;
}

inline void PrimaryGeneratorAction::SetPhi(G4double phi)
{
    fPhi = phi;
}

inline void PrimaryGeneratorAction::SetThetaRange(G4double lo, G4double hi)
{
    if (lo > -9999) fThetaLo = lo;

    if (hi > -9999) fThetaHi = hi;
}

inline void PrimaryGeneratorAction::SetEnpRange(G4double lo, G4double hi)
{
    if (lo > -9999) fEnpLo = lo;

    if (hi > -9999) fEnpHi = hi;
}

inline void PrimaryGeneratorAction::SetEventFile(G4String path)
{
    fEventFile = path;
}

inline void PrimaryGeneratorAction::SetPileUpProfile(G4String profile)
{
    fPileUpProfile = profile;
}

inline void PrimaryGeneratorAction::SetTargetProfile(G4String profile)
{
    fTargetProfile = profile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
