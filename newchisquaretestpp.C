void newchisquaretestpp(bool iseta = 1, bool isbk = 1, bool iseff = 1)
{

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	TString mcfile;
	TString datafile;
	TString bkfile;

	if (iseff)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/etacut_eff_file.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/efffile.root";
		}
	}
	if (!iseff)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/etacut_file.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/rawfile.root";
		}
	}

	datafile = "./new_pp_data_file_stability_readonly.root";

	if (isbk)
	{
		mcfile = "../ZBoson_18/rootfile/new_template_reco_gen_bk.root";
	}
	if (!isbk)
	{
		mcfile = "../ZBoson_18/rootfile/new_template_reco_gen_nobk.root";
	}

	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, iseta, isbk, iseff);

	// ovo->bincontentcheck(isbk);
	// ovo->RebinAllpp(2);
	ovo->calculatechisqpp(isbk);
	ovo->plottingandformattingpp(iseta, isbk, iseff);
}