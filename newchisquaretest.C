void newchisquaretest(int type = 1, bool isbk = 0, bool isNew = 1)
{
	// type 0 == raw , 1 == eta
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	//Don't forget to change params in header!!! Like binning ... 

	TString mcfile;
	TString datafile;
	TString bkfile;

	if (type == 0)
		bkfile = "../ZBoson_18/rootfile/normalized/rawfile.root";
	if (type == 1)
		bkfile = "../ZBoson_18/rootfile/normalized/etacut_file.root";

	datafile = "../ZBoson_18/rootfile/data_file.root";
	//datafile = "../ZBoson_18/rootfile/shift_-0.150_smear_0.0075_fixed_modified_signal_100.root";
	

	chisquaretest *ovo = new chisquaretest(datafile, bkfile, type,isNew);

	// ovo->bincontentcheck(isbk);
	ovo->RebinAll(type);

	ovo->calculatechisq(isbk);
	ovo->plottingandformatting(type, isbk);
}