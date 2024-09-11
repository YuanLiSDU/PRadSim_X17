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
// PRadSim.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite with ROOT support.
//   Mar 2017, C. Gu, Add DRad configuration.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysListEmModified.hh"
#include "PhysListPureEm.hh"
#include "RootTree.hh"
#include "SteppingVerbose.hh"

#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"

#include "G4BuilderType.hh"
#include "G4ios.hh"
#include "G4String.hh"
#include "Randomize.hh"

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <string>
#include <time.h>

#ifdef G4VIS_USE
    #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
    #include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootTree *gRootTree = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void usage(int, char **argv)
{
    printf("usage: %s [options] [MACRO_NAME]\n", argv[0]);
    printf("  -c, --conf=prad          Set configuration\n");
    printf("  -s, --seed=1             Set random seed\n");
    printf("  -p, --physics=FTFP_BERT  Set physics list\n");
    printf("  -h, --help               Print usage\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc, char **argv)
{
    std::string conf = "prad";
    std::string physics_list = "FTFP_BERT";
    std::string seed = "random";
    std::string macro;
    macro.clear();

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"conf", required_argument, 0, 'c'},
            {"physics", required_argument, 0, 'p'},
            {"seed", required_argument, 0, 's'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "c:hp:s:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'c':
            conf = optarg;
            break;

        case 'h':
            usage(argc, argv);
            exit(0);
            break;

        case 'p':
            physics_list = optarg;
            break;

        case 's':
            seed = optarg;
            break;

        case '?':
            // getopt_long already printed an error message
            break;

        default:
            usage(argc, argv);
            break;
        }
    }

    if (optind + 1 == argc)
        macro = argv[optind++];

    // Initialize the random engine
    G4Random::setTheEngine(new CLHEP::Ranlux64Engine);

    if (seed == "random") {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        long randseed =  ts.tv_sec * 1000000000L + ts.tv_nsec;
        std::cout << "Using random seed " << randseed << std::endl;
        G4Random::setTheSeed(randseed);
    } else
        G4Random::setTheSeed(stol(seed));

    // Initialize output root tree
    std::ifstream data_file("output/file.output");
    std::string file_name;
    int run_number = 1;

    if (data_file.is_open()) {
        std::string line;

        while (std::getline(data_file, line)) {
            if (line.at(0) == '#')
                continue;

            std::stringstream iss(line);
            iss >> file_name >> run_number;
            break;
        }

        data_file.close();
    }

    if (file_name.empty()) file_name = "simrun";

    std::ofstream data_file_o("output/file.output");
    data_file_o << file_name << "  " << run_number + 1;
    data_file_o.close();

    std::string path = "output/" + file_name + "_" + std::to_string(run_number) + ".root";
    char outf[256];
    strcpy(outf, path.c_str());

    gRootTree = new RootTree(outf);

    // Construct the default run manager
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    G4RunManager *runManager = new G4RunManager;

    // Set physics list
    bool pure_em = false;
    bool extra_em = false;
    bool local = false;

    if (physics_list.size() >= 2 && physics_list.substr(0, 2) == "EM") {
        pure_em = true;

        std::size_t found = physics_list.find_last_of("_");

        if (found != std::string::npos && physics_list.substr(found + 1) == "EXTRA") {
            extra_em = true;
            physics_list = physics_list.substr(0, found);
        }
    }

    if (physics_list.size() > 6) {
        std::size_t found = physics_list.find_last_of("_");

        if (found != std::string::npos && physics_list.substr(found + 1) == "LOCAL") {
            local = true;
            physics_list = physics_list.substr(0, found);
        }
    }

    G4VModularPhysicsList *physicsList = NULL;

    if (pure_em)
        physicsList = new PhysListPureEm(physics_list, extra_em, 1);
    else {
        G4PhysListFactory factory;

        if (!factory.IsReferencePhysList(physics_list)) {
            std::cout << "Physics list " << physics_list << " is not available in PhysListFactory." << std::endl;
            physics_list = "FTFP_BERT";
        }

        physicsList = factory.GetReferencePhysList(physics_list);
    }

    physicsList->RegisterPhysics(new G4StepLimiterPhysics());

    if (local) {
        physicsList->RemovePhysics(bElectromagnetic);
        physicsList->RegisterPhysics(new PhysListEmModified());
    }

    runManager->SetUserInitialization(physicsList);

    // Set mandatory initialization classes
    DetectorConstruction *detector = new DetectorConstruction(conf);
    runManager->SetUserInitialization(detector);

    ActionInitialization *action = new ActionInitialization(conf);
    runManager->SetUserInitialization(action);

    // Get the pointer to the User Interface manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (!macro.empty()) {
        // batch mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macro);
    } else {
        // interactive mode : define UI session
#ifdef G4VIS_USE
        // Initialize visualization
        G4VisManager *visManager = new G4VisExecutive;
        visManager->Initialize();
	std::cout << "init " << physics_list << " !!!" << std::endl;
#endif

#ifdef G4UI_USE
        G4UIExecutive *ui = new G4UIExecutive(argc, argv);
std::cout << "init UI" << physics_list << " !!!" << std::endl;
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute init_vis.mac");
#endif

        if (ui->IsGUI())
            UImanager->ApplyCommand("/control/execute gui.mac");

        ui->SessionStart();
        delete ui;
#endif

#ifdef G4VIS_USE
        delete visManager;
#endif
    }

    // Job termination
    delete runManager;
    delete gRootTree;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
