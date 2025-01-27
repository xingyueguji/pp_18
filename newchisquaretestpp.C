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
	// mcfile = Form("../ZBoson_18/rootfile/shift_%.1f_%.1f_smear_%.1f_%.3f_modified_signal_%i_%i_%i.root", shiftlowbin, shifthighbin, smearlowbin, smearhighbin, nbins_mass_shift, nbins_smear, numberofsamples);
	// mcfile = "../ZBoson_18/rootfile/pp_shift_-0.15_-0.1_smear_0.006_0.012_modified_signal_21_21_100.root";

	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, iseta, iseff);

	// ovo->bincontentcheck(isbk);
	// ovo->RebinAllpp(2);
	ovo->calculatechisqpp(isbk);
	ovo->plottingandformattingpp(iseta, isbk,iseff);
}