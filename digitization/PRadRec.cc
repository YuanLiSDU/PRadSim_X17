//
// PRadRec.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PRadDataHandler.h"
#include "PRadHyCalSystem.h"
//#include "PRadDensityCorrection.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TIterator.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct HyCalModule{
    double x;
    double y;
    double sizex;
    double sizey;
    int id;
    HyCalModule(double x1, double y1, double sx, double sy, int id1) : x(x1), y(y1), sizex(sx), sizey(sy), id(id1) {}
    HyCalModule() {}
    bool Contain(double x1, double y1){
        if (fabs(x1-x) <= sizex/2. && fabs(y1-y) <= sizey/2.) return true;
        else return false;
    }
};


void usage(int, char **argv)
{
    printf("usage: %s [options] FILE_NAME\n", argv[0]);
    printf("  -h, --help                 Print usage\n");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//**********nonlin correction***********//
#include "ConfigParser.h"
#define T_BLOCKS 2156
double nonlinConst[T_BLOCKS][3];
int    smearMode;
int    doSShapeCorr = 0;
TH2D* moduleECorrHist[T_BLOCKS];
vector<HyCalModule> moduleList;
TFile* moduleECorrFile;

void LoadConst(double e);
double GetModuleECorr(double x, double y);
Double_t EnergyCorrect(Double_t energy, Short_t cid);
//**************************************//

int GetPrimExID(string & s){
    if (s[0] == 'W' || s[0] == 'G'){
        int id = stoi(s.substr(1, s.length() - 1));
        if (s[0] == 'W') id += 1000;
        return id;
    }else{
        return -1;
    }
}

bool SortForR(HyCalModule const & a, HyCalModule const & b){
    return (a.x*a.x + a.y*a.y) < (b.x*b.x + b.y*b.y);
}

int main(int argc, char **argv)
{
    std::string filename;
    double ei = 2142;
    smearMode = 0;
    while (1) {
        static struct option long_options[] = {
            {"help",  no_argument, 0, 'h'},
            {"smear", required_argument, 0, 's'},
            {"energy",  required_argument, 0, 'e'},
	    {"s_shape_corr", required_argument, 0, 'd'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
            
        case 'e':
            ei = atof(optarg);

            if (ei < 1.0) {
                usage(argc, argv);
                exit(1);
            }

            break;
            
        case 's':
            smearMode = stoi(optarg);
	    cout<<"smear mode: "<<smearMode<<endl;
            break;
            
	case 'd':
	    if (stoi(optarg) != 0) doSShapeCorr = 1;
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
        filename = argv[optind++];
    else {
        usage(argc, argv);
        exit(1);
    }

    // simulation data is more like raw evio data with HyCal information only,
    // so we only need hycal system to connected to the handler
    PRadDataHandler *handler = new PRadDataHandler();
    PRadHyCalSystem *hycal = new PRadHyCalSystem("config/hycal.conf");

    int runNumber = 1288;
    string corrFile = "./database/density_params_sim_1GeV.txt";
    if (ei>1500.) { runNumber = 1443; corrFile = "./database/density_params_sim_2GeV.txt"; }

    //PRadDensityCorrection *corrector = NULL;
    //if (doDenCorr) corrector = new PRadDensityCorrection(hycal, runNumber, corrFile);

    handler->SetHyCalSystem(hycal);
    handler->ReadFromEvio(filename);

    std::string tf = filename;
    size_t ppos = tf.rfind('.', tf.length());

    if (ppos != std::string::npos)
        tf = tf.substr(0, ppos);

    std::string outf = tf + "_rec.root";

    TFile *f = new TFile(outf.c_str(), "RECREATE");
    TTree *t = new TTree("T", "Reconstructed Sim Results");
    
    //**********nonlin correction***************//
    LoadConst(ei);
    //******************************************//

    int N;
    double E[100], X[100], Y[100], Z[100]; // maximum number of clusters, 100 is enough
    UInt_t Flag[100];
    Short_t NModule[100];
    Short_t CID[100];
    // retrieve part of the cluster information
    t->Branch("HC.N", &N, "HC.N/I");
    t->Branch("HC.X", X, "HC.X[HC.N]/D");
    t->Branch("HC.Y", Y, "HC.Y[HC.N]/D");
    t->Branch("HC.Z", Z, "HC.Z[HC.N]/D");
    t->Branch("HC.P", E, "HC.P[HC.N]/D");
    t->Branch("HC.CID", CID, "HC.CID[HC.N]/S");
    t->Branch("HC.Flag", Flag, "HC.Flag[HC.N]/i");
    t->Branch("HC.NModule", NModule, "HC.NModule[HC.N]/S");
    

    int i = 1;

    for (auto &event : handler->GetEventData()) {
        if (i % 1000 == 0)
            std::cout << i << " events processed" << std::endl;

        hycal->Reconstruct(event);
        auto &hits = hycal->GetDetector()->GetHits();

        N = (int)hits.size();




        for (int j = 0; j < (int)hits.size(); ++j) {


            X[j] = hits[j].x;
            Y[j] = hits[j].y;
            Z[j] = 5646.12 - 3000.0 + 88.9 + hits[j].z;
            Flag[j] = hits[j].flag;
            NModule[j] = hits[j].nblocks;
            CID[j] = hits[j].cid;
            E[j] = EnergyCorrect(hits[j].E, hits[j].cid);
            double sShapeFactor = 1.;
            if (doSShapeCorr) sShapeFactor = GetModuleECorr(X[j], Y[j]);
            E[j] /= sShapeFactor;
        }

        t->Fill();

        i++;
    }

    f->cd();
    t->Write();
    f->Close();

    return 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

//**********nonlin correction***********//
void LoadConst(double e)
{
    ConfigParser parser;
    string fileName = "./database/calibration/2GeV_nonlin_const_LA_1.25.txt";

    if (e<1500) fileName = "./database/calibration/1GeV_nonlin_const_LA_1.125.txt";
    
    cout<<"using calibration file "<<fileName<<endl;
    
    if (!parser.ReadFile(fileName)){
        std::cout<<"cannot find mc calibration file"<<std::endl;
        exit(0);
    }
    
    int count = 0;
    while(parser.ParseLine()){
        for (int i=0; i<3; i++)   nonlinConst[count][i] = parser.TakeFirst().Double();
        count++;
    }
    //parser.CloseFile();

    if (!doSShapeCorr) return;

    cout<<"loading s shape correction file"<<endl;

    for (int i=0; i<T_BLOCKS; i++) moduleECorrHist[i] = 0;

    moduleList.clear();
    
    if (!parser.ReadFile("database/hycal_module.txt")){
        cout<<"cannot find hycal_module.txt"<<endl;
        exit(0);
    }
    
    string name, type;
    double sx, sy, sz, x, y, z;
    while (parser.ParseLine()){
        parser>>name>>type>>sx>>sy>>sz>>x>>y>>z;
        
        int id = GetPrimExID(name);
        if (id < 0) continue;
        
        moduleList.push_back(HyCalModule(x, y, sx, sy, id));
    }
    
    std::sort(moduleList.begin(), moduleList.end(), SortForR);

    cout<<e<<" "<<smearMode<<endl;

    if (e > 1500){
        if (smearMode == 0){
	    moduleECorrFile = new TFile("./database/calibration/2143MeV_ep_sim_total_LA_depth_LG_1.25_PWO_none_sCorr_10x10_redo.root", "READ");    
	}
	else{
	    moduleECorrFile = new TFile("./database/calibration/2143MeV_ee_sim_total_LA_depth_LG_1.25_PWO_none_sCorr_10x10.root", "READ");
	} 
    }
    else{
	if (smearMode == 0){
	    moduleECorrFile = new TFile("./database/calibration/1101MeV_ep_sim_total_LA_depth_LG_1.125_PWO_none_sCorr_10x10.root", "READ");
	}
	else{
	    moduleECorrFile = new TFile("./database/calibration/1101MeV_ee_sim_total_LA_depth_LG_1.125_PWO_none_sCorr_10x10.root", "READ");
	}
    }

    TIter next(moduleECorrFile->GetListOfKeys());
    TKey *thisKey;
    while ( (thisKey = (TKey*)next()) ){
        TClass *thisClass = gROOT->GetClass(thisKey->GetClassName());
        if (!thisClass->InheritsFrom("TDirectory")) continue;
        TDirectory *thisDir = (TDirectory*)thisKey->ReadObj();
        string name1 = thisKey->GetName();
        int id = stoi(name1);
        if (smearMode == 0){
            moduleECorrHist[id-1] = (TH2D*)thisDir->Get(Form("avg_ep_ratio_%04d", id));
        }else{
	    moduleECorrHist[id-1] = (TH2D*)thisDir->Get(Form("avg_ee_ratio_%04d", id));
	}
    }

}

double GetModuleECorr(double x, double y)
{
    float tx = 0;
    float ty = 0;
    int id = -1;
    
    for (unsigned int i=0; i<moduleList.size(); i++){
        if (moduleList[i].Contain(x, y)) {
            tx = (x - moduleList[i].x)/moduleList[i].sizex;
            ty = (y - moduleList[i].y)/moduleList[i].sizey;
            id = moduleList[i].id;
            break;
        }
    }
    if (id <= 0) return 1.;
    
    if (moduleECorrHist[id-1] == 0) {
        cout<<"module ECorr hist does not exist for "<<id<<" "<<endl;
        return 1.;
    }
    
    int binx = moduleECorrHist[id-1]->GetXaxis()->FindBin(tx);
    int biny = moduleECorrHist[id-1]->GetYaxis()->FindBin(ty);
    
    double factor = moduleECorrHist[id-1]->GetBinContent(binx, biny);

    
    if (factor < 0.7 || factor > 1.3) {cout<<"unusual ECorr factor for "<<id<<" "<<factor<<" "<<binx<<" "<<biny<<endl; return 1.; }
    return factor;
}


double EnergyCorrect(Double_t energy, Short_t cid)
{

    if (cid <= 0) return energy;

    double factor = (nonlinConst[cid-1][0] + nonlinConst[cid-1][1]/sqrt(energy/1000.) + nonlinConst[cid-1][2]/(energy/1000.));
    if (fabs(factor-1.) > 0.1) return energy;
    
    return energy*factor; 
}
//***************************************//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
