#include "chisquaretest.h"
void newchisquaretestpp(int type = 1, bool iseta = 1)
{

	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acoon
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
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_tnpU.root";
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
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_tnpD.root";
	}

	if (type == 4)
	{
		if (iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/Eta_acoon.root";
		}
		if (!iseta)
		{
			bkfile = "../ZBoson_18/rootfile/normalized/FA_acoon.root";
		}
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_acoon.root";
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
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_nominal_range.root";
	}

	if (type == 1)
	{
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_nominal.root";
	}

	if (type == 5)
	{
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_nominal_no_bksub.root";
	}

	if (type == 6)
	{
		zoominfile = "../ZBoson_18/rootfile/template_pp_zoomin_nominal_binning.root";
	}

	datafile = "./new_pp_data_file_stability_readonly.root";

	if (type == 7)
	{
		mcfile = "../ZBoson_18/rootfile/template_pp_mass_range.root";
	}
	else if (type != 5)
	{
		mcfile = "../ZBoson_18/rootfile/template_pp_bk.root";
	}
	else
	{
		mcfile = "../ZBoson_18/rootfile/template_pp_nobk.root";
	}

	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, type, iseta,zoominfile);

	// ovo->bincontentcheck(isbk);
	ovo->RebinAllpp(iseta, type);
	ovo->calculatechisqpp(isbk);
	ovo->plottingandformattingpp(type, iseta);
}