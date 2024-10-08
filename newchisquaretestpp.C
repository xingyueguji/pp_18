#include "chisquaretest.h"

void newchisquaretestpp(bool iseta = 1, bool isbk = 1)
{

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	TString mcfile;
	TString datafile;
	TString bkfile;

	if (iseta)
	{
		bkfile = "etacut_file.root";
	}
	if (!iseta)
	{
		bkfile = "rawfile.root";
	}

	datafile = "./new_pp_data_file_stability_readonly.root";
	mcfile = "modified_signal_test.root";

	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, iseta);

	// ovo->bincontentcheck(isbk);
	// ovo->RebinAllpp(4);
	ovo->calculatechisqpp(isbk);
	ovo->plottingandformattingpp(iseta, isbk);
}