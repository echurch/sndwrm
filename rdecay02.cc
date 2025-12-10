#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4HadronicParameters.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "Shielding.hh"
#include "G4OpticalPhysics.hh"

#include <cstdlib>
#include <ctime>

// Usage helper
namespace {
void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << "   app [-m macro ] [-u UIsession] [-r seed] [-p physList] " << G4endl;
}
}

int main(int argc, char** argv)
{
    // --------------------------------------------------------------------
    // Banner
    // --------------------------------------------------------------------
    G4cout << G4endl
    << "===============================================================" << G4endl
    << "           Underground LAr Detectors – GEANT4 App              " << G4endl
    << "===============================================================" << G4endl
    << G4endl;

    // --------------------------------------------------------------------
    // Parse arguments
    // --------------------------------------------------------------------
    G4UIExecutive* ui = nullptr;
    if (argc == 1) ui = new G4UIExecutive(argc, argv);

    G4String macro;
    G4String session;
    G4String physListName;
    G4long seed = time(nullptr);

    for (int i = 1; i < argc; i += 2) {
        if (G4String(argv[i]) == "-m" && i + 1 < argc) {
            macro = argv[i + 1];
        } else if (G4String(argv[i]) == "-u" && i + 1 < argc) {
            session = argv[i + 1];
        } else if (G4String(argv[i]) == "-r" && i + 1 < argc) {
            seed = std::atoi(argv[i + 1]);
        } else if (G4String(argv[i]) == "-p" && i + 1 < argc) {
            physListName = argv[i + 1];
        } else {
            PrintUsage();
            return 1;
        }
    }

    // --------------------------------------------------------------------
    // Random engine
    // --------------------------------------------------------------------
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4Random::setTheSeed(seed);

    // --------------------------------------------------------------------
    // RunManager (Serial)
    // --------------------------------------------------------------------
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

    // --------------------------------------------------------------------
    // Detector
    // --------------------------------------------------------------------
    auto* det = new DetectorConstruction;
    runManager->SetUserInitialization(det);

    // --------------------------------------------------------------------
    // Physics list selection
    // --------------------------------------------------------------------
    G4HadronicParameters::Instance()->SetTypeTablePT("njoy");
    G4HadronicParameters::Instance()->SetEnableNUDEX(true);

    G4PhysListFactory factory;
    G4VModularPhysicsList* physlist = nullptr;

    // If not provided, check PHYSLIST env
    if (!physListName.size()) {
        const char* envP = std::getenv("PHYSLIST");
        if (envP) physListName = envP;
    }

    // If still empty, fallback
    if (!physListName.size()) {
        physListName = "Shielding";
    }

    // If name exists in factory, use reference list
    if (factory.IsReferencePhysList(physListName)) {
        physlist = factory.GetReferencePhysList(physListName);
    } else {
        // Otherwise fallback to your previous Shielding + OpticalPhysics setup
        auto* shielding = new Shielding(1, "HP", "", true);
        shielding->RegisterPhysics(new G4OpticalPhysics());
        physlist = shielding;
    }

    runManager->SetUserInitialization(physlist);

    // --------------------------------------------------------------------
    // User actions
    // --------------------------------------------------------------------
    runManager->SetUserInitialization(new ActionInitialization(det));

    // Initialize
    runManager->Initialize();

    // --------------------------------------------------------------------
    // Visualization
    // --------------------------------------------------------------------
    G4VisManager* visManager = nullptr;
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (ui) {
        visManager = new G4VisExecutive;
        visManager->Initialize();

        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    } else {
        UImanager->ApplyCommand("/control/execute " + macro);
    }

    // Cleanup
    delete visManager;
    delete runManager;

    return 0;
}
