#include "chisquaretest.h"
void newchisquaretestpp(int type = 1, int version = 2, bool firsttime = 0)
{

	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = AcoUp
	// type 5 = AcoDown
	// type 6 = Nominal_no_bk
	// type 7 = massrange

	bool isbk = true;

	if (type == 6)
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
	TString prefix;

	if (version == 0)
	{
		prefix = "../ZBoson_18/rootfile/version_0/";
	}
	if (version == 1)
	{
		prefix = "../ZBoson_18/rootfile/version_1/";
	}
	if (version == 2)
	{
		prefix = "../ZBoson_18/rootfile/version_2/";
	}

	if (type == 1 || type == 6)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_nominal.root";
	}

	if (type == 2)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_tnpU.root";
		zoominfile = prefix + "template_pp_zoomin_tnpU.root";
	}

	if (type == 3)
	{

		bkfile = "../ZBoson_18/rootfile/normalized/FA_tnpD.root";
		zoominfile = prefix + "template_pp_zoomin_tnpD.root";
	}

	if (type == 4)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_acoup.root";
		zoominfile = prefix + "template_pp_zoomin_acoup.root";
	}

	if (type == 5)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_acodown.root";
		zoominfile = prefix + "template_pp_zoomin_acodown.root";
	}

	if (type == 7)
	{

		bkfile = "../ZBoson_18/rootfile/normalized/FA_mass_range.root";
		zoominfile = prefix + "template_pp_zoomin_nominal_range.root";
	}

	if (type == 1)
	{
		zoominfile = prefix + "template_pp_zoomin_nominal.root";
	}

	if (type == 6)
	{
		zoominfile = prefix + "template_pp_zoomin_nominal_no_bksub.root";
	}

	datafile = "./new_pp_data_file_stability_readonly.root";

	if (type == 7)
	{
		mcfile = prefix + "template_pp_mass_range.root";
	}
	else if (type != 6)
	{
		mcfile = prefix + "template_pp_bk.root";
	}
	else
	{
		mcfile = prefix + "template_pp_nobk.root";
	}

	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, type, zoominfile, version, firsttime);

	// ovo->bincontentcheck(isbk);
	ovo->RebinAllpp(type);
	ovo->calculatechisqpp(isbk);
	ovo->plottingandformattingpp(type, version);
}