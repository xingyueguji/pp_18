#include "chisquaretest.h"

void newchisquaretest(int type = 1, bool isbk = 1)
{
	// type 0 == raw , 1 == eta
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	TString mcfile;
	TString datafile;
	TString bkfile;

	if (type == 0)
		bkfile = "rawfile.root";
	if (type == 1)
		bkfile = "etacut_file.root";

	datafile = "data_file.root";
	mcfile = "modified_signal_test.root";

	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, type);

	// ovo->bincontentcheck(isbk);
	ovo->RebinAll(4);
	ovo->calculatechisq(isbk);
	ovo->plottingandformatting(type, isbk);
}