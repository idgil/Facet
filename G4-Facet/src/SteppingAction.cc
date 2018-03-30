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
// -------------------------------------------------------------------
// $Id: SteppingAction.cc,v 1.1 2012-09-21 mraine Exp $
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Analysis.hh"

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "stddef.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* RuAct,EventAction* event)
: G4UserSteppingAction(), runAction(RuAct), eventAction(event)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
     G4double flagParticle=0.;
     G4double energydep = 0.;
     G4double kinenergy = 0.;
     G4double PartType = 0.;
     G4double momentumx = 0.;
     G4double momentumy = 0.;
     G4double momentumz= 0.;
     G4double momentumxy2= 0.;
     G4double momentum2 = 0.;
     G4double angz = 0.;
     G4double angphi = 0.;
 //    G4int colimflag = 0;
     G4double edep = 0.;
     G4int detNo = -1.;
     G4int lineNo = -1;
     G4int layerNo = -1;
     G4double pi  = 3.14159265358979323846;

// In which volume point is located

  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetPhysicalVolume();
  G4String name = volume->GetName();

  G4StepPoint* point1 = step->GetPreStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();

 // G4int CopyNo = volume->GetCopyNo();
 // G4String MatName = volume->GetLogicalVolume()->GetMaterial()->GetName();
 // G4cout << "name =  " << "  " << name << "  copy = " << CopyNo << "  Metrial =" << MatName<< G4endl;

// Kill outgoing particles
 //  if (name == "Shield") step->GetTrack()->SetTrackStatus(fStopAndKill);

  const G4Track* track = step->GetTrack();
  if ( track->GetParentID() == 0 )
      {
      PartType = 1.; //Primary
      // Check if the particle just borned and is not a secondary
      if ( track->GetCurrentStepNumber() == 1)
         {

      energydep = step->GetTotalEnergyDeposit();
      kinenergy = step->GetTrack()->GetKineticEnergy();
      kinenergy = energydep+kinenergy;
      eventAction->AddInitSpectrum(kinenergy); //primary energy in current event
      //G4cout << "energy =  " << "  " << kinenergy << G4endl;

// Start angle of primary particle
      G4ThreeVector momentum = track->GetMomentumDirection();
   //   momentumx = momentum.x();
   //   momentumy = momentum.y();
   //   momentumz = momentum.z();
   //   momentum2 = sqrt(pow(momentumx,2)+pow(momentumy,2)+pow(momentumz,2));
   //   angz = acos(momentumz/momentum2)*180./pi;
      angz = momentum.theta()*180/pi;
      angphi= momentum.phi()*180/pi+180;
      eventAction->SetStartAngle(angz,angphi);
  //   G4cout << "angle =  " << "  " << momentum.phi()*180/pi+180<< G4endl;
         }


      }
   // If the secondary - find the kind of particle and its energy loss in detector
  else
     {
        if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "e-")      	  flagParticle = 1;
        if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "proton")  	  flagParticle = 2;
        if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "GenericIon")    flagParticle = 3;
        if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "gamma")    flagParticle = 4;
        if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ->GetParticleName() == "neutron")    flagParticle = 5;
    }
        edep =0.;
        edep = step->GetTotalEnergyDeposit();

 // Detector number with signal
    if (name == "pixel_PV" && edep > 0.)
    {

     detNo= step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
  //   G4cout << "pointx =  " << "  " << point1->GetPosition().getX()*mm<< G4endl;
  //   G4cout << "pointy =  " << "  " << point1->GetPosition().getY()*mm<< G4endl;
  //   G4cout << "pointz =  " << "  " << point1->GetPosition().getZ()*mm<< G4endl;
    G4VPhysicalVolume* mother = touch1->GetVolume(1);
 //      G4VPhysicalVolume* grandmother = touch1->GetVolume(3);
     lineNo = mother->GetCopyNo();
  //   layerNo = grandmother->GetCopyNo();
     G4String layername = touch1->GetVolume(3)->GetName();
     if (layername == "CCD1_PV") layerNo = 0;
     if (layername == "CCD2_PV") layerNo = 1;
     if (layername == "CCD3_PV") layerNo = 2;

      eventAction->AddDetEdep(detNo,lineNo, layerNo,edep);

       if (point1->GetStepStatus() == fGeomBoundary || track->GetCurrentStepNumber() == 1) // if particle entered new volume or just borned in volume
       {
           G4ThreeVector momentum = track->GetMomentumDirection();
           angz=momentum.theta()*180/pi;
           angphi=momentum.phi()*180/pi+180;
//           eventAction->SetAngle(angz, angphi);
           eventAction->AddAngle(detNo,lineNo, layerNo,angz, angphi);
       }
    }

}    
