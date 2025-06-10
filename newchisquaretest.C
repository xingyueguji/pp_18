#include "chisquaretest.h"
void newchisquaretest(int type = 1, int version = 0)
{
	// version 0 = legacy, version 1 = direct fill, version 2 = modification
	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acooff
	// type 5 = Nominal_no_bk
	// type 6 = rebinned
	// type 7 = massrange
	// type 8 = HF up
	// type 9 = HF down

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
	TString prefix;
	if (version == 0)
	{
		mcfile = "../ZBoson_18/rootfile/version_0/template_PbPb.root";
		prefix = "../ZBoson_18/rootfile/version_0/";
	}
	if (version == 1)
	{
		mcfile = "../ZBoson_18/rootfile/version_1/template_PbPb.root";
		prefix = "../ZBoson_18/rootfile/version_1/";
	}
	if (version == 2)
	{
		mcfile = "../ZBoson_18/rootfile/version_2/template_PbPb.root";
		prefix = "../ZBoson_18/rootfile/version_2/";
	}

	if (type == 1 || type == 5 || type == 6)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_nominal.root";
	}

	if (type == 2)
	{

		bkfile = "../ZBoson_18/rootfile/normalized/FA_tnpU.root";
		zoominfile = prefix + "template_PbPb_zoomin_tnpU.root";
	}

	if (type == 3)
	{

		bkfile = "../ZBoson_18/rootfile/normalized/FA_tnpD.root";
		zoominfile = prefix + "template_PbPb_zoomin_tnpD.root";
	}

	if (type == 4)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_acoon.root";
		zoominfile = prefix + "template_PbPb_zoomin_acoon.root";
	}

	if (type == 7)
	{

		bkfile = "../ZBoson_18/rootfile/normalized/FA_mass_range.root";
		zoominfile = prefix + "template_PbPb_zoomin_nominal_range.root";
	}

	if (type == 8)
	{

		bkfile = "../ZBoson_18/rootfile/normalized/FA_HF_up.root";
		zoominfile = prefix + "template_PbPb_zoomin_HF_up.root";
	}

	if (type == 9)
	{
		bkfile = "../ZBoson_18/rootfile/normalized/FA_HF_down.root";
		zoominfile = prefix + "template_PbPb_zoomin_HF_down.root";
	}

	if (type == 1)
	{
		zoominfile = prefix + "template_PbPb_zoomin_nominal.root";
	}

	if (type == 5)
	{
		zoominfile = prefix + "template_PbPb_zoomin_nominal_no_bksub.root";
	}

	if (type == 6)
	{
		zoominfile = prefix + "template_PbPb_zoomin_nominal_binning.root";
	}

	datafile = "../ZBoson_18/rootfile/data_file.root";
	chisquaretest *ovo = new chisquaretest(mcfile, datafile, bkfile, type, zoominfile);

	// ovo->bincontentcheck(isbk);
	ovo->RebinAll(type);
	ovo->calculatechisq(isbk, 0);
	ovo->calculatechisq(isbk, 1);
	ovo->plottingandformatting(type);
}