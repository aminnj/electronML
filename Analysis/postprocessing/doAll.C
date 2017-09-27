{
    gROOT->ProcessLine(".L CMS3.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");
    TChain *ch = new TChain("tree");

    ch->Add("output_12.root");
    ScanChain(ch);
}
