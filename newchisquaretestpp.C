#include "chisquaretest.h"

void newchisquaretestpp(int datatype = 1, bool isbk = 1){

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	TString mcfile;
	TString datafile;
	TString bkfile;

	bkfile = "../ZBoson_18/rootfile/normalized/etacut_eff_file.root";

	datafile = "./pp_data_file_stability.root";
	mcfile = "../ZBoson_18/rootfile/modified_signal.root";

	chisquaretest* ovo = new chisquaretest(mcfile,datafile,bkfile);



	//ovo->bincontentcheck(isbk);
	//ovo->RebinAllpp(4);
	ovo->calculatechisqpp(isbk);
	ovo->plottingandformattingpp(isbk);

}