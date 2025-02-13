void newchisquaretest(int type = 1, bool iseta = 1)
{
	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acooff
	// type 5 = Nominal_no_bk
	// type 6 = rebinned
	// type 7 = massrange

	bool isbk = true;
	if (type == 5)
	{
		isbk = false;
	}

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();

	// Don't forget to change params in header!!! Like binning ...

	TString mcfile;
	TString datafile;
	TString bkfile;
	TString zoominfile;

	if (type == 1 || type == 5 || type == 6)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/Eta_nominal.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/FA_nominal.root";
		}
	}

	if (type == 2)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/Eta_tnpU.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/FA_tnpU.root";
		}
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_tnpU.root";
	}

	if (type == 3)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/Eta_tnpD.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/FA_tnpD.root";
		}
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_tnpD.root";
	}

	if (type == 4)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/Eta_acooff.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/FA_acooff.root";
		}
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_acooff.root";
	}

	if (type == 7)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/Eta_mass_range.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/FA_mass_range.root";
		}
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_nominal_range.root";
	}

	if (type == 1)
	{
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_nominal.root";
	}

	if (type == 5)
	{
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_nominal_no_bksub.root";
	}

	if (type == 6)
	{
		zoominfile = "../ZBoson_18/rootfile/template_PbPb_zoomin_nominal_binning.root";
	}

	datafile = "../ZBoson_18/rootfile/data_file.root";
	// datafile = "../ZBoson_18/rootfile/shift_-0.150_smear_0.0075_fixed_modified_signal_100.root";

	chisquaretest *ovo = new chisquaretest(datafile, bkfile, type, iseta,zoominfile);

	// ovo->bincontentcheck(isbk);
	ovo->RebinAll(iseta, type);
	ovo->calculatechisq(isbk,0);
	ovo->calculatechisq(isbk,1);
	ovo->plottingandformatting(type, iseta);
}