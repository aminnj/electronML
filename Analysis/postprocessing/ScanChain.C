#pragma GCC diagnostic ignored "-Wsign-compare"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1.h"
#include "TChain.h"

#include "CMS3.h"
#include "c2numpy.h"

using namespace std;
using namespace tas;


int ScanChain(TChain *ch){

    c2numpy_writer writer;
    c2numpy_init(&writer, "output", 90000000);

    c2numpy_addcolumn(&writer, "ismatch", C2NUMPY_INTC);

    c2numpy_addcolumn(&writer, "matchtype", C2NUMPY_INTC);

    c2numpy_addcolumn(&writer, "mva", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "pt", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "eta", C2NUMPY_FLOAT32);

    c2numpy_addcolumn(&writer, "ele_kfhits", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_oldsigmaietaieta", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_oldsigmaiphiiphi", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_oldcircularity", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_oldr9", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_scletawidth", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_sclphiwidth", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_oldhe", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_psEoverEraw", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_kfchi2", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_chi2_hits", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_fbrem", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_ep", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_eelepout", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_IoEmIop", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_deltaetain", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_deltaphiin", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_deltaetaseed", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "scl_eta", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_gsfhits", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_expectedMissingInnerHits", C2NUMPY_FLOAT32);
    c2numpy_addcolumn(&writer, "ele_convVtxFitProbability", C2NUMPY_FLOAT32);

    c2numpy_open(&writer);

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);

    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("tree");
        cms3.Init(tree);

        TString filename(currentFile->GetTitle());

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {


            cms3.GetEntry(event);
            nEventsTotal++;

            CMS3::progress(nEventsTotal, nEventsChain);

            for (int icell = 0; icell < rhs_e().size(); icell++) {
                std::cout <<  " rhs_e()[icell]: " << rhs_e()[icell] <<  " rhs_iphi()[icell]: " << rhs_iphi()[icell] <<  " rhs_ieta()[icell]: " << rhs_ieta()[icell] <<  std::endl;
            }

            
            break;

            bool ismatch = -1;
            int matchtype = -1; // -1 = ???, 0 = light fake, 1 = heavy fake, 2 = prompt lepton
            float mva = ele_ID1();
            float pt = ele_pt();
            float eta = ele_eta();
            int theircode = mc_ele_matchedFromCB();
            if (theircode == 0) { // bkg unmatch
                ismatch = false;
                matchtype = 0;
            } else if (theircode == 1) { // sig
                ismatch = true;
                matchtype = 2;
            } else if (theircode == 3) { // bkg from b,c
                ismatch = false;
                matchtype = 1;
            } else {
                std::cout << "shouldn't get here, can't have code=" << theircode << std::endl;
            }

            float ele_kfhits_           = ele_kfhits();
            float ele_oldsigmaietaieta_ = ele_oldsigmaietaieta();
            float ele_oldsigmaiphiiphi_ = ele_oldsigmaiphiiphi();
            float ele_oldcircularity_   = ele_oldcircularity();
            float ele_oldr9_            = ele_oldr9();
            float ele_scletawidth_      = ele_scletawidth();
            float ele_sclphiwidth_      = ele_sclphiwidth();
            float ele_oldhe_            = ele_oldhe();
            float ele_psEoverEraw_      = ele_psEoverEraw();
            float ele_kfchi2_           = ele_kfchi2();
            float ele_chi2_hits_        = ele_gsfchi2();
            float ele_fbrem_            = ele_fbrem();
            float ele_ep_               = ele_ep();
            float ele_eelepout_         = ele_eelepout();
            float ele_IoEmIop_          = ele_IoEmIop();
            float ele_deltaetain_       = ele_deltaetain();
            float ele_deltaphiin_       = ele_deltaphiin();
            float ele_deltaetaseed_     = ele_deltaetaseed();
            float scl_eta_              = scl_eta();
            float ele_gsfhits_                  = ele_gsfhits();
            float ele_expectedMissingInnerHits_ = ele_expected_inner_hits();
            float ele_convVtxFitProbability_    = ele_conversionVertexFitProbability();

            //bindVariables
            if(ele_fbrem_ < -1.) ele_fbrem_ = -1.;
            ele_deltaetain_ = fabs(ele_deltaetain_);
            if(ele_deltaetain_ > 0.06) ele_deltaetain_ = 0.06;
            ele_deltaphiin_ = fabs(ele_deltaphiin_);
            if(ele_deltaphiin_ > 0.6) ele_deltaphiin_ = 0.6;
            if(ele_ep_ > 20.) ele_ep_ = 20.;
            if(ele_eelepout_ > 20.) ele_eelepout_ = 20.;
            ele_deltaetaseed_ = fabs(ele_deltaetaseed_);
            if(ele_deltaetaseed_ > 0.2) ele_deltaetaseed_ = 0.2;
            if(ele_oldcircularity_ < -1.) ele_oldcircularity_ = -1;
            if(ele_oldcircularity_ > 2.) ele_oldcircularity_ = 2.;
            if(ele_oldr9_ > 5) ele_oldr9_ = 5;
            if(ele_chi2_hits_ > 200.) ele_chi2_hits_ = 200;
            if(ele_kfchi2_ > 10.) ele_kfchi2_ = 10.;
            if(std::isnan(ele_oldsigmaiphiiphi_)) ele_oldsigmaiphiiphi_ = 0.;

            c2numpy_intc(&writer, ismatch);
            c2numpy_intc(&writer, matchtype);
            c2numpy_float32(&writer, mva);
            c2numpy_float32(&writer, pt);
            c2numpy_float32(&writer, eta);
            c2numpy_float32(&writer, ele_kfhits_);
            c2numpy_float32(&writer, ele_oldsigmaietaieta_);
            c2numpy_float32(&writer, ele_oldsigmaiphiiphi_);
            c2numpy_float32(&writer, ele_oldcircularity_);
            c2numpy_float32(&writer, ele_oldr9_);
            c2numpy_float32(&writer, ele_scletawidth_);
            c2numpy_float32(&writer, ele_sclphiwidth_);
            c2numpy_float32(&writer, ele_oldhe_);
            c2numpy_float32(&writer, ele_psEoverEraw_);
            c2numpy_float32(&writer, ele_kfchi2_);
            c2numpy_float32(&writer, ele_chi2_hits_);
            c2numpy_float32(&writer, ele_fbrem_);
            c2numpy_float32(&writer, ele_ep_);
            c2numpy_float32(&writer, ele_eelepout_);
            c2numpy_float32(&writer, ele_IoEmIop_);
            c2numpy_float32(&writer, ele_deltaetain_);
            c2numpy_float32(&writer, ele_deltaphiin_);
            c2numpy_float32(&writer, ele_deltaetaseed_);
            c2numpy_float32(&writer, scl_eta_);
            c2numpy_float32(&writer, ele_gsfhits_);
            c2numpy_float32(&writer, ele_expectedMissingInnerHits_);
            c2numpy_float32(&writer, ele_convVtxFitProbability_);

        }//event loop

        delete file;
    }//file loop

    c2numpy_close(&writer);

    return 0;

}

