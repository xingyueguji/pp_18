#include "plotting_helper.h"
#include "tdrStyle.C"
#include "CMS_lumi.C"

class chisquaretest : public plotting_helper
{
public:
	chisquaretest();
	chisquaretest(TString s2, TString s3, int type, bool iseta, TString s4);
	chisquaretest(TString s1, TString s2, TString s3, int type, bool iseta, TString s4);
	~chisquaretest();
	Double_t myownfunctionchi2(TH1D *h1, TH1D *h2);
	void calculatechisq(bool isbk, bool iszoomin = 0);
	void calculatechisqpp(bool isbk);
	void plottingandformatting(int type, bool iseta);
	void plottingandformattingpp(int type, bool iseta);
	void bincontentcheck(bool isbk);
	void RebinAll(bool iseta, int type);
	TGraph *RemoveInvalidPoints(TGraph *originalGraph);
	TGraph *CombineGraphsToCircle(TGraph *graph1, TGraph *graph2);
	void RebinAllpp(bool iseta, int x);
	Double_t getuncertainty(TH2D *h_1, int type, Int_t minBinX, Int_t minBinY);
	void getcontour(TH2D *h1, int type, Int_t minBinX, Int_t minBinY, Double_t minBinContent, Double_t *arrayleft, Double_t *arrayright);
	void drawcontour(TGraph *onesig_left, TGraph *onesig_right, TGraph *twosig_left, TGraph *twosig_right, int iteration, bool iseta, bool ispp, TString contourtitle, TString contoursaving);
	std::vector<double> createCustomBinning(
		TH1D *hist,
		double range1_min, double range1_max, int rebin1,
		double range2_min, double range2_max, int rebin2,
		double range3_min, double range3_max, int rebin3,
		double range4_min = -1, double range4_max = -1, int rebin4 = 1);

	Int_t mapthreecases(int type, bool isbk, bool iseff);
	void saveChi2Region(TH2D *hist, int binX_min, int binY_min, int my_case, bool iseta, int iteration, int region_size = 3, bool ispp = true);
	TH1D *reducebin(TH1D *h_old_data, bool iseta, bool isleft);
	TH1D *reducebinpp(TH1D *h_old_data, bool isleft);
	void readlimit(int type, bool iseta, int cent);
	void readlimitpp(int type, bool iseta);

	static const int nbins_mass_shift = 42;
	static const int nbins_smear = 42;
	static const int nbins_cent = 11;

	// those default values will be changed.

	double lowbin_mass_shift = 0.0;
	double highbin_mass_shift = 0.0;
	double lowbin_smear = 0.0;
	double highbin_smear = 0.0;

	double h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
	double h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

	double h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
	double h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;

	TH1D *h_mc_signal[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *h_mc_signal_zoomin[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *h_data[nbins_cent];
	TH1D *h_mc_bk[nbins_cent];
	TH1D *h_data_bksub[nbins_cent];

	TH2D *h_chisquare[nbins_cent];
	TH2D *h_chisquare_zoomin[nbins_cent];

	Double_t dMass_HI[5];
	Double_t dMass_Err_HI[5];

	Double_t dWidth_HI[5];
	Double_t dWidth_Err_HI[5];

	Double_t xposition_HI[5];
	Double_t xposition_err_HI[5];

	Double_t contour_x_left_onesig_HI[nbins_smear];
	Double_t contour_x_right_onesig_HI[nbins_smear];

	Double_t contour_x_left_before_onesig_HI[nbins_smear] = {};
	Double_t contour_x_left_after_onesig_HI[nbins_smear] = {};
	Double_t contour_x_right_before_onesig_HI[nbins_smear] = {};
	Double_t contour_x_right_after_onesig_HI[nbins_smear] = {};

	Double_t contour_x_left_weighted_onesig[nbins_smear] = {};
	Double_t contour_x_right_weighted_onesig[nbins_smear] = {};

	Double_t contour_x_left_before_twosig_HI[nbins_smear] = {};
	Double_t contour_x_left_after_twosig_HI[nbins_smear] = {};
	Double_t contour_x_right_before_twosig_HI[nbins_smear] = {};
	Double_t contour_x_right_after_twosig_HI[nbins_smear] = {};

	Double_t contour_x_left_weighted_twosig[nbins_smear] = {};
	Double_t contour_x_right_weighted_twosig[nbins_smear] = {};

	Double_t localmin_x = 0;
	Double_t localmin_y = 0;

	Double_t contour_y_HI[nbins_cent][nbins_smear] = {};

	Double_t contour_x_left_twosig_HI[6];
	Double_t contour_x_right_twosig_HI[6];

	TH1D *h_mc_signal_pp[nbins_mass_shift][nbins_smear][23];
	TH1D *h_mc_signal_pp_zoomin[nbins_mass_shift][nbins_smear];
	TH1D *h_data_pp[23];
	TH1D *h_mc_bk_pp;
	TH1D *h_data_bksub_pp[23];

	TH2D *h_chisquare_pp[23];
	TH2D *h_chisquare_pp_zoomin;

	Double_t dMass_pp[23];
	Double_t dMass_Err_pp[23];

	Double_t dWidth_pp[23];
	Double_t dWidth_Err_pp[23];

	Double_t xposition[23];
	Double_t xposition_err[23];

	Double_t xposition_pp_inclusive = 1;
	Double_t xposition_pp_inclusive_err = 0;

	TGraphErrors *g_pp_dmass;
	TGraphErrors *g_pp_dwidth;

	TGraphErrors *g_HI_dmass;
	TGraphErrors *g_HI_dwidth;
	TGraph *g_HI_contour_1sig_left;
	TGraph *g_HI_contour_1sig_right;
	TGraph *g_HI_contour_2sig_left;
	TGraph *g_HI_contour_2sig_right;

	double chisquarearray[nbins_mass_shift][nbins_smear][nbins_cent] = {};

	TString mcfilepath;
	TString datafilepath;
	TString bkfilepath;

	TCanvas *c_2d_chisquare_ndf[nbins_cent];
	TCanvas *c_2d_chisquare_ndf_zoomin[nbins_cent];
	TCanvas *c_contour_HI[nbins_cent];
	// What I need: invariant mass comparison data vs mc raw, data vs data bk sub.
	TCanvas *c_data_mc_raw[nbins_cent];
	TCanvas *c_data_data_bk[nbins_cent];

	TCanvas *c_2d_chisquare_ndf_pp[23];
	TCanvas *c_2d_chisquare_ndf_pp_zoomin;
	TCanvas *c_data_mc_raw_pp[23];
	TCanvas *c_data_data_bk_pp[23];

	TFile *mcfile;
	TFile *datafile;
	TFile *bkfile;
	TFile *mcfile_zoomin;

	double eta_mass_shift_array_low[nbins_cent] = {-0.4, -0.25, -0.4, -0.45, 0, 0, 0, 0, 0, 0, -0.28};
	double eta_mass_shift_array_high[nbins_cent] = {-0.08, 0.15, 0.075, -0.05, 0, 0, 0, 0, 0, 0, -0.07};
	double eta_mass_smear_array_low[nbins_cent] = {-0.25, -0.05, -0.15, -0.3, 0, 0, 0, 0, 0, 0, 0.03};
	double eta_mass_smear_array_high[nbins_cent] = {0.3, 0.65, 0.65, 0.5, 0, 0, 0, 0, 0, 0, 0.36};

	double raw_mass_shift_array_low[nbins_cent] = {-0.24, -0.32, -0.32, -0.24, 0, 0, 0, 0, 0, 0, -0.22};
	double raw_mass_shift_array_high[nbins_cent] = {-0.02, -0.06, -0.04, 0.02, 0, 0, 0, 0, 0, 0, -0.08};
	double raw_mass_smear_array_low[nbins_cent] = {0.05, 0.175, -0.1, 0.2, 0, 0, 0, 0, 0, 0, 0.22};
	double raw_mass_smear_array_high[nbins_cent] = {0.4, 0.6, 0.4, 0.7, 0, 0, 0, 0, 0, 0, 0.45};

	double placeholder_mass_shift_array_low_zoomin;
	double placeholder_mass_shift_array_high_zoomin;
	double placeholder_mass_smear_array_low_zoomin;
	double placeholder_mass_smear_array_high_zoomin;

	double placeholder_pp_mass_shift_array_low_zoomin;
	double placeholder_pp_mass_shift_array_high_zoomin;
	double placeholder_pp_mass_smear_array_low_zoomin;
	double placeholder_pp_mass_smear_array_high_zoomin;

	// eta with bk

	double eta_pp_mass_shift_low_with_bk = -0.14;
	double eta_pp_mass_shift_high_with_bk = -0.06;
	double eta_pp_smear_low_with_bk = 0.2;
	double eta_pp_smear_high_with_bk = 0.32;

	// eta without bk

	double eta_pp_mass_shift_low_without_bk = -0.15;
	double eta_pp_mass_shift_high_without_bk = -0.06;
	double eta_pp_smear_low_without_bk = 0.23;
	double eta_pp_smear_high_without_bk = 0.37;

	// eta mass range

	double eta_pp_mass_shift_low_mass_range = -0.14;
	double eta_pp_mass_shift_high_mass_range = -0.07;
	double eta_pp_smear_low_mass_range = 0.13;
	double eta_pp_smear_high_mass_range = 0.25;

	// raw with bk

	double raw_pp_mass_shift_low_with_bk = -0.105;
	double raw_pp_mass_shift_high_with_bk = -0.05;
	double raw_pp_smear_low_with_bk = 0.2;
	double raw_pp_smear_high_with_bk = 0.3;

	// raw without bk

	double raw_pp_mass_shift_low_without_bk = -0.09;
	double raw_pp_mass_shift_high_without_bk = -0.05;
	double raw_pp_smear_low_without_bk = 0.28;
	double raw_pp_smear_high_without_bk = 0.36;

	// raw mass range

	double raw_pp_mass_shift_low_mass_range = -0.1;
	double raw_pp_mass_shift_high_mass_range = -0.05;
	double raw_pp_smear_low_mass_range = 0.16;
	double raw_pp_smear_high_mass_range = 0.22;
};
chisquaretest::chisquaretest()
{
}

chisquaretest::chisquaretest(TString s2, TString s3, int type, bool iseta, TString s4)
{

	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acooff
	// type 5 = Nominal_no_bk
	// type 6 = rebinned
	// type 7 = massrange

	datafilepath = s2;
	bkfilepath = s3;
	TString mcfilepathzoomin = s4;

	mcfile = new TFile("../ZBoson_18/rootfile/template_PbPb.root", "READ");
	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath, "READ");
	mcfile_zoomin = new TFile(mcfilepathzoomin, "READ");

	for (int i = 0; i < 5; i++)
	{
		xposition_HI[i] = i + 1;
		xposition_err_HI[i] = 0;
	}

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		if (!((cent < 4) || (cent == 10)))
			continue;

		if (iseta)
		{
			this->lowbin_mass_shift = eta_mass_shift_array_low[cent];
			this->highbin_mass_shift = eta_mass_shift_array_high[cent];
			this->lowbin_smear = eta_mass_smear_array_low[cent];
			this->highbin_smear = eta_mass_smear_array_high[cent];
		}
		if (!iseta)
		{
			this->lowbin_mass_shift = raw_mass_shift_array_low[cent];
			this->highbin_mass_shift = raw_mass_shift_array_high[cent];
			this->lowbin_smear = raw_mass_smear_array_low[cent];
			this->highbin_smear = raw_mass_smear_array_high[cent];
		}

		this->h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
		this->h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

		this->h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		this->h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;

		// Here's for zoomin

		readlimit(type, iseta, cent);

		double h_low_mass_shift_zoomin = placeholder_mass_shift_array_low_zoomin - ((placeholder_mass_shift_array_high_zoomin - placeholder_mass_shift_array_low_zoomin) / (nbins_mass_shift - 1)) / 2;
		double h_high_mass_shift_zoomin = placeholder_mass_shift_array_high_zoomin + ((placeholder_mass_shift_array_high_zoomin - placeholder_mass_shift_array_low_zoomin) / (nbins_mass_shift - 1)) / 2;
		double h_low_smear_zoomin = placeholder_mass_smear_array_low_zoomin - ((placeholder_mass_smear_array_high_zoomin - placeholder_mass_smear_array_low_zoomin) / (nbins_smear - 1)) / 2;
		double h_high_smear_zoomin = placeholder_mass_smear_array_high_zoomin + ((placeholder_mass_smear_array_high_zoomin - placeholder_mass_smear_array_low_zoomin) / (nbins_smear - 1)) / 2;

		for (int i = 0; i < nbins_smear; i++)
		{
			double bincenter = ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
			contour_y_HI[cent][i] = 2 * (i)*bincenter + lowbin_smear;
		}

		c_2d_chisquare_ndf[cent] = new TCanvas(Form("c_2d_chisquare_ndf_%i", cent), "", 3200, 2400);
		c_2d_chisquare_ndf_zoomin[cent] = new TCanvas(Form("c_2d_chisquare_ndf_zoomin_%i", cent), "", 3200, 2400);
		c_data_mc_raw[cent] = new TCanvas(Form("c_data_mc_raw_%i", cent), "", 800, 800);
		c_data_data_bk[cent] = new TCanvas(Form("c_data_data_bk_%i", cent), "", 800, 600);
		c_contour_HI[cent] = new TCanvas(Form("c_contour_HI_%i", cent), "", 800, 600);

		if (type == 1 || type == 5 || type == 6)
		{
			if (iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("Eta_nominal_%i", cent));
			}
			if (!iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("FA_nominal_%i", cent));
			}
		}

		if (type == 2)
		{
			if (iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("Eta_tnpU_%i", cent));
			}
			if (!iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("FA_tnpU_%i", cent));
			}
		}

		if (type == 3)
		{
			if (iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("Eta_tnpD_%i", cent));
			}
			if (!iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("FA_tnpD_%i", cent));
			}
		}

		if (type == 4)
		{
			if (iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("Eta_AcoOff_%i", cent));
			}
			if (!iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("FA_AcoOff_%i", cent));
			}
		}

		if (type == 7)
		{
			if (iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("Eta_mass_range_%i", cent));
			}
			if (!iseta)
			{
				h_data[cent] = (TH1D *)datafile->Get(Form("FA_mass_range_%i", cent));
			}
		}

		h_mc_bk[cent] = (TH1D *)bkfile->Get(Form("Normalized_mc_bk_%i", cent));
		h_chisquare[cent] = new TH2D(Form("h_chisquare_%i", cent), Form("Cent_%i_%i", this->cenlowlimit[cent], this->cenhighlimit[cent]), nbins_mass_shift, h_low_mass_shift, h_high_mass_shift, nbins_smear, h_low_smear, h_high_smear); // Equation here: half bin to the left and half bin to the right, bin width = range / (21-1)
		h_chisquare_zoomin[cent] = new TH2D(Form("h_chisquare_%i_zoomin", cent), "", nbins_mass_shift, h_low_mass_shift_zoomin, h_high_mass_shift_zoomin, nbins_smear, h_low_smear_zoomin, h_high_smear_zoomin);

		this->areanormalize(h_data[cent]);

		h_data_bksub[cent] = (TH1D *)h_data[cent]->Clone();
		h_data_bksub[cent]->Add(h_mc_bk[cent], -1);

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				if (type == 1 || type == 2 || type == 3 || type == 5 || type == 6)
				{
					if (iseta)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("template_Eta_nominal_%i_%i_%i", shift, smear, cent));
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)mcfile_zoomin->Get(Form("template_Eta_nominal_%i_%i_%i", shift, smear, cent));
					}
					if (!iseta)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("template_FA_nominal_%i_%i_%i", shift, smear, cent));
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)mcfile_zoomin->Get(Form("template_FA_nominal_%i_%i_%i", shift, smear, cent));
					}
				}

				if (type == 4)
				{
					if (iseta)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("template_Eta_acooff_%i_%i_%i", shift, smear, cent));
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)mcfile_zoomin->Get(Form("template_Eta_acooff_%i_%i_%i", shift, smear, cent));
					}
					if (!iseta)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("template_FA_acooff_%i_%i_%i", shift, smear, cent));
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)mcfile_zoomin->Get(Form("template_FA_acooff_%i_%i_%i", shift, smear, cent));
					}
				}

				if (type == 7)
				{
					if (iseta)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("template_Eta_mass_range_%i_%i_%i", shift, smear, cent));
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)mcfile_zoomin->Get(Form("template_Eta_mass_range_%i_%i_%i", shift, smear, cent));
					}
					if (!iseta)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("template_FA_mass_range_%i_%i_%i", shift, smear, cent));
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)mcfile_zoomin->Get(Form("template_FA_mass_range_%i_%i_%i", shift, smear, cent));
					}
				}
				this->areanormalize(h_mc_signal[shift][smear][cent]);
				this->areanormalize(h_mc_signal_zoomin[shift][smear][cent]);
			}
		}
	}
}

chisquaretest::chisquaretest(TString s1, TString s2, TString s3, int type, bool iseta, TString s4)
{

	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acooff
	// type 5 = Nominal_no_bk
	// type 6 = rebinned
	// type 7 = massrange

	mcfilepath = s1;
	datafilepath = s2;
	bkfilepath = s3;
	TString mcfile_zoomin_path = s4;

	cout << "mcfilepath is " << mcfilepath << endl;

	mcfile = new TFile(mcfilepath, "READ");
	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath, "READ");
	mcfile_zoomin = new TFile(mcfile_zoomin_path, "READ");

	if (iseta)
	{
		if (type == 5)
		{
			this->lowbin_mass_shift = eta_pp_mass_shift_low_without_bk;
			this->highbin_mass_shift = eta_pp_mass_shift_high_without_bk;
			this->lowbin_smear = eta_pp_smear_low_without_bk;
			this->highbin_smear = eta_pp_smear_high_without_bk;
		}
		else if (type == 7)
		{
			this->lowbin_mass_shift = eta_pp_mass_shift_low_mass_range;
			this->highbin_mass_shift = eta_pp_mass_shift_high_mass_range;
			this->lowbin_smear = eta_pp_smear_low_mass_range;
			this->highbin_smear = eta_pp_smear_high_mass_range;
		}
		else
		{
			this->lowbin_mass_shift = eta_pp_mass_shift_low_with_bk;
			this->highbin_mass_shift = eta_pp_mass_shift_high_with_bk;
			this->lowbin_smear = eta_pp_smear_low_with_bk;
			this->highbin_smear = eta_pp_smear_high_with_bk;
		}
	}

	if (!iseta)
	{
		if (type == 5)
		{
			this->lowbin_mass_shift = raw_pp_mass_shift_low_without_bk;
			this->highbin_mass_shift = raw_pp_mass_shift_high_without_bk;
			this->lowbin_smear = raw_pp_smear_low_without_bk;
			this->highbin_smear = raw_pp_smear_high_without_bk;
		}
		else if (type == 7)
		{
			this->lowbin_mass_shift = raw_pp_mass_shift_low_mass_range;
			this->highbin_mass_shift = raw_pp_mass_shift_high_mass_range;
			this->lowbin_smear = raw_pp_smear_low_mass_range;
			this->highbin_smear = raw_pp_smear_high_mass_range;
		}
		else
		{
			this->lowbin_mass_shift = raw_pp_mass_shift_low_with_bk;
			this->highbin_mass_shift = raw_pp_mass_shift_high_with_bk;
			this->lowbin_smear = raw_pp_smear_low_with_bk;
			this->highbin_smear = raw_pp_smear_high_with_bk;
		}
	}

	this->h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
	this->h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
	this->h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
	this->h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;

	readlimitpp(type, iseta);

	double h_low_mass_shift_zoomin = placeholder_pp_mass_shift_array_low_zoomin - ((placeholder_pp_mass_shift_array_high_zoomin - placeholder_pp_mass_shift_array_low_zoomin) / (nbins_mass_shift - 1)) / 2;
	double h_high_mass_shift_zoomin = placeholder_pp_mass_shift_array_high_zoomin + ((placeholder_pp_mass_shift_array_high_zoomin - placeholder_pp_mass_shift_array_low_zoomin) / (nbins_mass_shift - 1)) / 2;
	double h_low_smear_zoomin = placeholder_pp_mass_smear_array_low_zoomin - ((placeholder_pp_mass_smear_array_high_zoomin - placeholder_pp_mass_smear_array_low_zoomin) / (nbins_smear - 1)) / 2;
	double h_high_smear_zoomin = placeholder_pp_mass_smear_array_high_zoomin + ((placeholder_pp_mass_smear_array_high_zoomin - placeholder_pp_mass_smear_array_low_zoomin) / (nbins_smear - 1)) / 2;

	cout << "x low is " << h_low_mass_shift_zoomin << "x high is " << h_high_mass_shift_zoomin << "y low is " << h_low_smear_zoomin << "y high is " << h_high_smear_zoomin << endl;

	h_mc_bk_pp = (TH1D *)bkfile->Get("Normalized_mc_bk_10");

	for (int i = 0; i < nbins_smear; i++)
	{
		double bincenter = ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		// FIX ME !!!
		contour_y_HI[10][i] = 2 * (i)*bincenter + lowbin_smear;
	}

	for (int runperiod = 0; runperiod < 23; runperiod++)
	{
		if (runperiod == 22)
		{
			c_2d_chisquare_ndf_pp_zoomin = new TCanvas(Form("c_2d_chisquare_ndf_pp_zoomin_%i", 22), "", 3200, 2400);
		}
		c_2d_chisquare_ndf_pp[runperiod] = new TCanvas(Form("c_2d_chisquare_ndf_pp_%i", runperiod), "", 3200, 2400);
		c_data_mc_raw_pp[runperiod] = new TCanvas(Form("c_data_mc_raw_pp_%i", runperiod), "", 800, 600);
		c_data_data_bk_pp[runperiod] = new TCanvas(Form("c_data_data_bk_pp_%i", runperiod), "", 800, 600);

		if (type == 1 || type == 5 || type == 6)
		{
			if (iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("Eta_nominal_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("Eta_nominal_%i", runperiod));
			}
			if (!iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("FA_nominal_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("FA_nominal_%i", runperiod));
			}
		}

		if (type == 2)
		{
			if (iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("Eta_tnpU_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("Eta_tnpU_%i", runperiod));
			}
			if (!iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("FA_tnpU_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("FA_tnpU_%i", runperiod));
			}
		}

		if (type == 3)
		{
			if (iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("Eta_tnpD_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("Eta_tnpD_%i", runperiod));
			}
			if (!iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("FA_tnpD_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("FA_tnpD_%i", runperiod));
			}
		}

		if (type == 4)
		{
			if (iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("Eta_AcoOff_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("Eta_AcoOff_%i", runperiod));
			}
			if (!iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("FA_AcoOff_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("FA_AcoOff_%i", runperiod));
			}
		}

		if (type == 7)
		{
			if (iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("Eta_mass_range_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("Eta_mass_range_%i", runperiod));
			}
			if (!iseta)
			{
				if (runperiod == 22)
					h_data_pp[runperiod] = (TH1D *)datafile->Get("FA_mass_range_inclusive");
				else
					h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("FA_mass_range_%i", runperiod));
			}
		}

		h_chisquare_pp[runperiod] = new TH2D(Form("h_chisquare_pp_%i", runperiod), Form("h_chisquare_pp_%i", runperiod), nbins_mass_shift, h_low_mass_shift, h_high_mass_shift, nbins_smear, h_low_smear, h_high_smear);
		if (runperiod == 22)
			h_chisquare_pp_zoomin = new TH2D(Form("h_chisquare_pp_zoomin_%i", 22), "", nbins_mass_shift, h_low_mass_shift_zoomin, h_high_mass_shift_zoomin, nbins_smear, h_low_smear_zoomin, h_high_smear_zoomin);
		this->areanormalize(h_data_pp[runperiod]);

		h_data_bksub_pp[runperiod] = (TH1D *)h_data_pp[runperiod]->Clone();
		h_data_bksub_pp[runperiod]->Add(h_mc_bk_pp, -1);
		h_data_bksub_pp[runperiod]->SetName(Form("%s_bksub", h_data_pp[runperiod]->GetName()));

		xposition[runperiod] = runperiod + 1;
		xposition_err[runperiod] = 0;

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				if (type == 1 || type == 2 || type == 3 || type == 5 || type == 6)
				{
					if (iseta)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("template_Eta_nominal_%i_%i_%i", shift, smear, 10));
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)mcfile_zoomin->Get(Form("template_Eta_nominal_%i_%i_%i", shift, smear, 10));
						}
					}
					if (!iseta)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("template_FA_nominal_%i_%i_%i", shift, smear, 10));
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)mcfile_zoomin->Get(Form("template_FA_nominal_%i_%i_%i", shift, smear, 10));
						}
					}
				}

				if (type == 4)
				{
					if (iseta)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("template_Eta_acooff_%i_%i_%i", shift, smear, 10));
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)mcfile_zoomin->Get(Form("template_Eta_acooff_%i_%i_%i", shift, smear, 10));
						}
					}
					if (!iseta)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("template_FA_acooff_%i_%i_%i", shift, smear, 10));
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)mcfile_zoomin->Get(Form("template_FA_acooff_%i_%i_%i", shift, smear, 10));
						}
					}
				}
				if (type == 7)
				{
					if (iseta)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("template_Eta_mass_range_%i_%i_%i", shift, smear, 10));
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)mcfile_zoomin->Get(Form("template_Eta_mass_range_%i_%i_%i", shift, smear, 10));
						}
					}
					if (!iseta)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("template_FA_mass_range_%i_%i_%i", shift, smear, 10));
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)mcfile_zoomin->Get(Form("template_FA_mass_range_%i_%i_%i", shift, smear, 10));
						}
					}
				}

				this->areanormalize(h_mc_signal_pp[shift][smear][runperiod]);
				if (runperiod == 22)
				{
					this->areanormalize(h_mc_signal_pp_zoomin[shift][smear]);
				}
			}
		}
	}
}

void chisquaretest::calculatechisq(bool isbk, bool iszoomin = 0)
{

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		if (!((cent < 4) || (cent == 10)))
			continue;

		cout << "Now running File " << cenlowlimit[cent] << " " << cenhighlimit[cent] << endl;
		for (int i = 1; i <= h_data[cent]->GetNbinsX(); i++)
		{
			if (TMath::Power((h_data[cent]->GetBinError(i)), 2) == 0)
				cout << "Zero Error for bin " << i << endl;
		}

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				double chisquarevalue = 0;
				if (!iszoomin)
				{
					if (isbk)
						chisquarevalue = myownfunctionchi2(h_data_bksub[cent], h_mc_signal[shift][smear][cent]);
					if (!isbk)
						chisquarevalue = myownfunctionchi2(h_data[cent], h_mc_signal[shift][smear][cent]);
				}
				if (iszoomin)
				{
					if (isbk)
						chisquarevalue = myownfunctionchi2(h_data_bksub[cent], h_mc_signal_zoomin[shift][smear][cent]);
					if (!isbk)
						chisquarevalue = myownfunctionchi2(h_data[cent], h_mc_signal_zoomin[shift][smear][cent]);
				}
				std::ostringstream stream;
				stream << std::fixed << std::setprecision(4) << chisquarevalue;
				double formattedBinContent = std::stod(stream.str());
				if (!iszoomin)
					h_chisquare[cent]->SetBinContent(shift + 1, smear + 1, formattedBinContent);
				if (iszoomin)
					h_chisquare_zoomin[cent]->SetBinContent(shift + 1, smear + 1, formattedBinContent);
			}
		}
	}
}

void chisquaretest::calculatechisqpp(bool isbk)
{

	for (int runperiod = 0; runperiod < 23; runperiod++)
	{
		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				double chisquarevalue = 0;
				double chisquarevalue_zoomin = 0;
				if (isbk)
				{
					chisquarevalue = myownfunctionchi2(h_data_bksub_pp[runperiod], h_mc_signal_pp[shift][smear][runperiod]);
					if (runperiod == 22)
					{
						chisquarevalue_zoomin = myownfunctionchi2(h_data_bksub_pp[22], h_mc_signal_pp_zoomin[shift][smear]);
					}
				}

				if (!isbk)
				{
					chisquarevalue = myownfunctionchi2(h_data_pp[runperiod], h_mc_signal_pp[shift][smear][runperiod]);

					if (runperiod == 22)
					{
						chisquarevalue_zoomin = myownfunctionchi2(h_data_pp[22], h_mc_signal_pp_zoomin[shift][smear]);
					}
				}

				std::ostringstream stream;
				std::ostringstream stream_zoomin;

				stream << std::fixed << std::setprecision(4) << chisquarevalue;
				stream_zoomin << std::fixed << std::setprecision(4) << chisquarevalue_zoomin;

				double formattedBinContent = std::stod(stream.str());
				double formattedBinContent_zoomin = std::stod(stream_zoomin.str());
				h_chisquare_pp[runperiod]->SetBinContent(shift + 1, smear + 1, formattedBinContent);

				if (runperiod == 22)
				{
					h_chisquare_pp_zoomin->SetBinContent(shift + 1, smear + 1, formattedBinContent_zoomin);
				}
			}
		}
	}
}

void chisquaretest::plottingandformatting(int type, bool iseta)
{
	TString chi2_title;
	TString chi2_saving_path;
	TString chi2_saving_path_zoomin;
	TString data_mc_title;
	TString data_data_title;
	TString data_mc_saving_path;
	TString data_data_saving_path;
	TString contour_saving_path;

	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acooff
	// type 5 = nominal without bk subtraction

	// Int_t my_case = mapthreecases(type, isbk, iseff);

	if (type == 1)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, Nominal, centrality: (%i-%i)";
			data_mc_title = "Eta cut, Nominal, Cent:(%i-%i)";
			data_data_title = "Eta_Nominal_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/Nominal/Eta_Nominal_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/Nominal/Eta_Nominal_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/Nominal/Eta_Nominal_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_Nominal_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/Eta_Nominal_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, Nominal, centrality: (%i-%i)";
			data_mc_title = "WholeAcceptance, Nominal, Cent:(%i-%i)";
			data_data_title = "Raw_Nominal_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/Nominal/Raw_Nominal_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/Nominal/Raw_Nominal_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/Nominal/Raw_Nominal_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_Nominal_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_Nominal_%i_%i.png";
		}
	}
	if (type == 2)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, tnpU, centrality: (%i-%i)";
			data_mc_title = "Eta, tnpU, Cent:(%i-%i)";
			data_data_title = "Eta_tnpU_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/tnpU/Raw_tnpU_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/tnpU/Raw_tnpU_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/tnpU/Raw_tnpU_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_tnpU_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/Eta_tnpU_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, tnpU, centrality: (%i-%i)";
			data_mc_title = "WholeAcceptance, tnpU, Cent:(%i-%i)";
			data_data_title = "Raw_tnpU_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/tnpU/Raw_tnpU_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/tnpU/Raw_tnpU_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/tnpU/Raw_tnpU_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_tnpU_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_tnpU_%i_%i.png";
		}
	}
	if (type == 3)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, tnpD, centrality: (%i-%i)";
			data_mc_title = "Eta cut, tnpD, Cent:(%i-%i)";
			data_data_title = "Eta_tnpD_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/tnpD/Eta_tnpD_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/tnpD/Eta_tnpD_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/tnpD/Eta_tnpD_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_tnpD_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/eta_tnpD_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, tnpD, centrality: (%i-%i)";
			data_mc_title = "WholeAcceptance, tnpD, Cent:(%i-%i)";
			data_data_title = "Raw_tnpD_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/tnpD/Raw_tnpD_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/tnpD/Raw_tnpD_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/tnpD/Raw_tnpD_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_tnpD_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_tnpD_%i_%i.png";
		}
	}
	if (type == 4)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, Acooff, centrality: (%i-%i)";
			data_mc_title = "Eta cut, Acooff, Cent:(%i-%i)";
			data_data_title = "Eta_Acooff_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/Acooff/Eta_Acooff_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/Acooff/Eta_Acooff_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/Acooff/Eta_Acooff_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_Acooff_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/Eta_Acooff_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, Acooff, centrality: (%i-%i)";
			data_mc_title = "WholeAcceptance, Acooff, Cent:(%i-%i)";
			data_data_title = "Raw_Acooff_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/Acooff/Raw_Acooff_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/Acooff/Raw_Acooff_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/Acooff/Raw_Acooff_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_Acooff_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_Acooff_%i_%i.png";
		}
	}
	if (type == 5)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, Nominal no bk, centrality: (%i-%i)";
			data_mc_title = "Eta cut, Nominal no bk, centrality: (%i-%i)";
			data_data_title = "Eta_Nominal_no_bk_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/Nominal_no_bk/Eta_nominal_no_bk_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/Nominal_no_bk/Eta_nominal_no_bk_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/Nominal_no_bk/Eta_nominal_no_bk_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_nominal_no_bk_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/Eta_nominal_no_bk_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, Nominal no bk, centrality: (%i-%i)";
			data_mc_title = "|#eta| < 2.4, Nominal no bk, centrality: (%i-%i)";
			data_data_title = "Raw_Nominal_no_bk_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/Nominal_no_bk/Raw_nominal_no_bk_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/Nominal_no_bk/Raw_nominal_no_bk_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/Nominal_no_bk/Raw_nominal_no_bk_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_nominal_no_bk_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_nominal_no_bk_%i_%i.png";
		}
	}
	if (type == 6)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, Nominal uniform rebin, centrality: (%i-%i)";
			data_mc_title = "Eta cut, Nominal uniform rebin, centrality: (%i-%i)";
			data_data_title = "Eta_Nominal_uniform_rebin_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/Nominal_uniform_rebin/Eta_nominal_uniform_rebin_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/Nominal_uniform_rebin/Eta_nominal_uniform_rebin_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/Nominal_uniform_rebin/Eta_nominal_uniform_rebin_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_nominal_uniform_rebin_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/Eta_nominal_uniform_rebin_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, Nominal uniform rebin, centrality: (%i-%i)";
			data_mc_title = "|#eta| < 2.4, Nominal uniform rebin, centrality: (%i-%i)";
			data_data_title = "Raw_Nominal_uniform_rebin_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/Nominal_uniform_rebin/Raw_nominal_uniform_rebin_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/Nominal_uniform_rebin/Raw_nominal_uniform_rebin_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/Nominal_uniform_rebin/Raw_nominal_uniform_rebin_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_nominal_uniform_rebin_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_nominal_uniform_rebin_%i_%i.png";
		}
	}
	if (type == 7)
	{
		if (iseta)
		{
			chi2_title = "PbPb, |#eta| < 1.0, Nominal mass range, centrality: (%i-%i)";
			data_mc_title = "Eta cut, Nominal mass range, centrality: (%i-%i)";
			data_data_title = "Eta_Nominal_mass_range_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/eta/Nominal_mass_range/Eta_nominal_mass_range_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/eta/Nominal_mass_range/Eta_nominal_mass_range_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/eta/Nominal_mass_range/Eta_nominal_mass_range_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_nominal_nominal_mass_range_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/Eta_nominal_nominal_mass_range_%i_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "PbPb, |#eta| < 2.4, Nominal mass range, centrality: (%i-%i)";
			data_mc_title = "|#eta| < 2.4, Nominal mass range, centrality: (%i-%i)";
			data_data_title = "Raw_Nominal_mass_range_%i_%i";

			chi2_saving_path = "./newchi2/chi2plots/raw/Nominal_mass_range/Raw_nominal_mass_range_%i_%i.png";
			chi2_saving_path_zoomin = "./newchi2/chi2plots/raw/Nominal_mass_range/Raw_nominal_mass_range_%i_%i_zoomin.png";
			data_mc_saving_path = "./newchi2/datamc/raw/Nominal_mass_range/Raw_nominal_mass_range_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_nominal_mass_range_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/raw_nominal_mass_range_%i_%i.png";
		}
	}

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		if (!((cent < 4) || (cent == 10)))
			continue;

		// Chi2/ndf plot
		c_2d_chisquare_ndf[cent]->cd();
		gStyle->SetPalette(kRainBow);
		h_chisquare[cent]->SetTitle(Form(chi2_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		h_chisquare[cent]->SetTitleFont(42);
		h_chisquare[cent]->Draw("COLZ");
		// h_chisquare[cent]->Draw("TEXTSAME");
		h_chisquare[cent]->GetXaxis()->SetNdivisions(21, 0, 0);
		h_chisquare[cent]->GetYaxis()->SetNdivisions(21, 0, 0);
		h_chisquare[cent]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare[cent]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare[cent]->GetXaxis()->SetTitle("Mass Shifted Amount (GeV)");
		h_chisquare[cent]->GetYaxis()->SetTitle("Width Shifted Amount (GeV)");

		for (int j = 1; j <= nbins_mass_shift; j++)
		{
			for (int k = 1; k <= nbins_smear; k++)
			{
				double xlow = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(k);
				double xup = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(k);
				double ylow = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(j);
				double yup = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(j);

				TBox *box = new TBox(xlow, ylow, xup, yup);
				box->SetFillStyle(0);	   // No fill
				box->SetLineColor(kBlack); // Black border
				box->SetLineWidth(1);	   // Border width

				box->Draw("same");
			}
		}

		Int_t minBinX = -1, minBinY = -1;
		Double_t minContent = h_chisquare[cent]->GetMaximum();
		for (Int_t binX = 1; binX <= h_chisquare[cent]->GetNbinsX(); ++binX)
		{
			for (Int_t binY = 1; binY <= h_chisquare[cent]->GetNbinsY(); ++binY)
			{
				Double_t content = h_chisquare[cent]->GetBinContent(binX, binY);
				if (content < minContent)
				{
					minContent = content;
					minBinX = binX;
					minBinY = binY;
				}
			}
		}

		// This is to export minimum region
		this->saveChi2Region(h_chisquare[cent], minBinX, minBinY, type, iseta, cent, 3, false);
		this->getcontour(h_chisquare[cent], 1, minBinX, minBinY, minContent, contour_x_left_onesig_HI, contour_x_right_onesig_HI);

		g_HI_contour_1sig_left = new TGraph(nbins_smear, contour_x_left_weighted_onesig, contour_y_HI[cent]);
		g_HI_contour_1sig_right = new TGraph(nbins_smear, contour_x_right_weighted_onesig, contour_y_HI[cent]);
		g_HI_contour_2sig_left = new TGraph(nbins_smear, contour_x_left_weighted_twosig, contour_y_HI[cent]);
		g_HI_contour_2sig_right = new TGraph(nbins_smear, contour_x_right_weighted_twosig, contour_y_HI[cent]);

		bool iseta = false;
		if (type == 1)
			iseta = true;

		TString T_contour_saving_path = Form(contour_saving_path, cenlowlimit[cent], cenhighlimit[cent]);
		TString T_chi2_title = Form(chi2_title, cenlowlimit[cent], cenhighlimit[cent]);

		this->drawcontour(g_HI_contour_1sig_left, g_HI_contour_1sig_right, g_HI_contour_2sig_left, g_HI_contour_2sig_right, cent, iseta, false, T_chi2_title, T_contour_saving_path);

		Double_t xMin = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(minBinX);
		Double_t xMax = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(minBinX);
		Double_t yMin = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(minBinY);
		Double_t yMax = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(minBinY);

		Double_t xCenter = h_chisquare[cent]->GetXaxis()->GetBinCenter(minBinX);
		Double_t yCenter = h_chisquare[cent]->GetYaxis()->GetBinCenter(minBinY);

		// Here's drawing all the left right TBoxes.

		for (int nbinssmearing = 1; nbinssmearing <= nbins_smear; ++nbinssmearing)
		{

			c_2d_chisquare_ndf[cent]->cd();
			// cout << "contour_x_left_before_onesig_HI " << nbinssmearing - 1 << " is " << contour_x_left_before_onesig_HI[nbinssmearing-1] << endl;
			Double_t xMin_left_1 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_left_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMin_left_2 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_left_after_onesig_HI[nbinssmearing - 1]);

			Double_t xMax_left_1 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_left_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMax_left_2 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_left_after_onesig_HI[nbinssmearing - 1]);

			Double_t yMin_left_1 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_left_2 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_left_1 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_left_2 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			Double_t xMin_right_1 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_right_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMin_right_2 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_right_after_onesig_HI[nbinssmearing - 1]);

			Double_t xMax_right_1 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_right_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMax_right_2 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_right_after_onesig_HI[nbinssmearing - 1]);

			Double_t yMin_right_1 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_right_2 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_right_1 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_right_2 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			Double_t xMin_twosig_left_1 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_left_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMin_twosig_left_2 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_left_after_twosig_HI[nbinssmearing - 1]);

			Double_t xMax_twosig_left_1 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_left_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMax_twosig_left_2 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_left_after_twosig_HI[nbinssmearing - 1]);

			Double_t yMin_twosig_left_1 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_twosig_left_2 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_twosig_left_1 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_twosig_left_2 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			Double_t xMin_twosig_right_1 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_right_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMin_twosig_right_2 = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(contour_x_right_after_twosig_HI[nbinssmearing - 1]);

			Double_t xMax_twosig_right_1 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_right_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMax_twosig_right_2 = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(contour_x_right_after_twosig_HI[nbinssmearing - 1]);

			Double_t yMin_twosig_right_1 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_twosig_right_2 = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_twosig_right_1 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_twosig_right_2 = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			TBox *box_left_1 = new TBox(xMin_left_1, yMin_left_1, xMax_left_1, yMax_left_1);
			box_left_1->SetLineColor(kGreen);
			box_left_1->SetLineWidth(5);
			box_left_1->SetFillStyle(0);
			box_left_1->Draw("same");

			TBox *box_left_2 = new TBox(xMin_left_2, yMin_left_2, xMax_left_2, yMax_left_2);
			box_left_2->SetLineColor(kGreen);
			box_left_2->SetLineWidth(5);
			box_left_2->SetFillStyle(0);
			box_left_2->Draw("same");

			TBox *box_right_1 = new TBox(xMin_right_1, yMin_right_1, xMax_right_1, yMax_right_1);
			box_right_1->SetLineColor(kOrange);
			box_right_1->SetLineWidth(5);
			box_right_1->SetFillStyle(0);
			box_right_1->Draw("same");

			TBox *box_right_2 = new TBox(xMin_right_2, yMin_right_2, xMax_right_2, yMax_right_2);
			box_right_2->SetLineColor(kOrange);
			box_right_2->SetLineWidth(5);
			box_right_2->SetFillStyle(0);
			box_right_2->Draw("same");

			TBox *box_twosig_left_1 = new TBox(xMin_twosig_left_1, yMin_twosig_left_1, xMax_twosig_left_1, yMax_twosig_left_1);
			box_twosig_left_1->SetLineColor(kMagenta);
			box_twosig_left_1->SetLineWidth(5);
			box_twosig_left_1->SetFillStyle(0);
			box_twosig_left_1->Draw("same");

			TBox *box_twosig_left_2 = new TBox(xMin_twosig_left_2, yMin_twosig_left_2, xMax_twosig_left_2, yMax_twosig_left_2);
			box_twosig_left_2->SetLineColor(kMagenta);
			box_twosig_left_2->SetLineWidth(5);
			box_twosig_left_2->SetFillStyle(0);
			box_twosig_left_2->Draw("same");

			TBox *box_twosig_right_1 = new TBox(xMin_twosig_right_1, yMin_twosig_right_1, xMax_twosig_right_1, yMax_twosig_right_1);
			box_twosig_right_1->SetLineColor(kCyan);
			box_twosig_right_1->SetLineWidth(5);
			box_twosig_right_1->SetFillStyle(0);
			box_twosig_right_1->Draw("same");

			TBox *box_twosig_right_2 = new TBox(xMin_twosig_right_2, yMin_twosig_right_2, xMax_twosig_right_2, yMax_twosig_right_2);
			box_twosig_right_2->SetLineColor(kCyan);
			box_twosig_right_2->SetLineWidth(5);
			box_twosig_right_2->SetFillStyle(0);
			box_twosig_right_2->Draw("same");
		}

		double x_cross = 0;
		double y_cross = 0;

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
		box1->SetLineColor(kRed);
		box1->SetLineWidth(4);
		box1->SetFillStyle(0);
		box1->Draw("same");

		// This is Zoomed in version

		c_2d_chisquare_ndf_zoomin[cent]->cd();
		gStyle->SetPalette(kRainBow);
		h_chisquare_zoomin[cent]->SetTitle(Form(chi2_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		h_chisquare_zoomin[cent]->SetTitleFont(42);
		h_chisquare_zoomin[cent]->Draw("COLZ");
		// h_chisquare_zoomin[cent]->Draw("TEXTSAME");
		h_chisquare_zoomin[cent]->GetXaxis()->SetNdivisions(21, 0, 0);
		h_chisquare_zoomin[cent]->GetYaxis()->SetNdivisions(21, 0, 0);
		h_chisquare_zoomin[cent]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare_zoomin[cent]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare_zoomin[cent]->GetXaxis()->SetTitle("Mass Shifted Amount (GeV)");
		h_chisquare_zoomin[cent]->GetYaxis()->SetTitle("Width Shifted Amount (GeV)");

		for (int j = 1; j <= nbins_mass_shift; j++)
		{
			for (int k = 1; k <= nbins_smear; k++)
			{
				double xlow = h_chisquare_zoomin[cent]->GetXaxis()->GetBinLowEdge(k);
				double xup = h_chisquare_zoomin[cent]->GetXaxis()->GetBinUpEdge(k);
				double ylow = h_chisquare_zoomin[cent]->GetYaxis()->GetBinLowEdge(j);
				double yup = h_chisquare_zoomin[cent]->GetYaxis()->GetBinUpEdge(j);

				TBox *box = new TBox(xlow, ylow, xup, yup);
				box->SetFillStyle(0);	   // No fill
				box->SetLineColor(kBlack); // Black border
				box->SetLineWidth(1);	   // Border width

				box->Draw("same");
			}
		}

		Int_t minBinX_zoomin = -1, minBinY_zoomin = -1;
		Double_t minContent_zoomin = h_chisquare_zoomin[cent]->GetMaximum();

		for (Int_t binX = 1; binX <= h_chisquare_zoomin[cent]->GetNbinsX(); ++binX)
		{
			for (Int_t binY = 1; binY <= h_chisquare_zoomin[cent]->GetNbinsY(); ++binY)
			{
				Double_t content = h_chisquare_zoomin[cent]->GetBinContent(binX, binY);
				if (content < minContent_zoomin)
				{
					minContent_zoomin = content;
					minBinX_zoomin = binX;
					minBinY_zoomin = binY;
				}
			}
		}

		Double_t xMin_zoomin = h_chisquare_zoomin[cent]->GetXaxis()->GetBinLowEdge(minBinX_zoomin);
		Double_t xMax_zoomin = h_chisquare_zoomin[cent]->GetXaxis()->GetBinUpEdge(minBinX_zoomin);
		Double_t yMin_zoomin = h_chisquare_zoomin[cent]->GetYaxis()->GetBinLowEdge(minBinY_zoomin);
		Double_t yMax_zoomin = h_chisquare_zoomin[cent]->GetYaxis()->GetBinUpEdge(minBinY_zoomin);

		Double_t xCenter_zoomin = h_chisquare_zoomin[cent]->GetXaxis()->GetBinCenter(minBinX_zoomin);
		Double_t yCenter_zoomin = h_chisquare_zoomin[cent]->GetYaxis()->GetBinCenter(minBinY_zoomin);

		TBox *box1_zoomin = new TBox(xMin_zoomin, yMin_zoomin, xMax_zoomin, yMax_zoomin);
		box1_zoomin->SetLineColor(kRed);
		box1_zoomin->SetLineWidth(4);
		box1_zoomin->SetFillStyle(0);
		box1_zoomin->Draw("same");

		if (cent == 0 || cent == 1 || cent == 2 || cent == 3)
		{
			dMass_HI[cent + 1] = xCenter_zoomin;
			dWidth_HI[cent + 1] = yCenter_zoomin;

			// cout << "dWidth is " << dWidth_HI[cent + 1] << endl;

			dMass_Err_HI[cent + 1] = getuncertainty(h_chisquare[cent], 1, minBinX, minBinY);
			dWidth_Err_HI[cent + 1] = getuncertainty(h_chisquare[cent], 2, minBinX, minBinY);

			double diffindM = xCenter_zoomin - xCenter;
			dMass_Err_HI[cent + 1] = dMass_Err_HI[cent + 1] + diffindM;

			double diffindW = yCenter_zoomin - yCenter;
			dWidth_Err_HI[cent + 1] = dWidth_Err_HI[cent + 1] + diffindW;

			x_cross = dMass_Err_HI[cent + 1];
			y_cross = dWidth_Err_HI[cent + 1];
		}

		if (cent == 10)
		{
			dMass_HI[0] = xCenter_zoomin;
			dWidth_HI[0] = yCenter_zoomin;

			dMass_Err_HI[0] = getuncertainty(h_chisquare[cent], 1, minBinX, minBinY);
			dWidth_Err_HI[0] = getuncertainty(h_chisquare[cent], 2, minBinX, minBinY);

			double diffindM = xCenter_zoomin - xCenter;
			dMass_Err_HI[0] = dMass_Err_HI[0] + diffindM;

			double diffindW = yCenter_zoomin - yCenter;
			dWidth_Err_HI[0] = dWidth_Err_HI[0] + diffindW;

			x_cross = dMass_Err_HI[0];
			y_cross = dWidth_Err_HI[0];
		}

		c_2d_chisquare_ndf[cent]->cd();

		TMarker *marker_L = new TMarker(xCenter_zoomin - x_cross, yCenter_zoomin, 20); // Marker type 29 (big star)
		marker_L->SetMarkerColor(kRed);
		marker_L->SetMarkerSize(3.5);
		marker_L->Draw("same");

		TMarker *marker_R = new TMarker(xCenter_zoomin + x_cross, yCenter_zoomin, 20); // Marker type 29 (big star)
		marker_R->SetMarkerColor(kRed);
		marker_R->SetMarkerSize(3.5);
		marker_R->Draw("same");

		TMarker *marker_T = new TMarker(xCenter_zoomin, yCenter_zoomin - y_cross, 20); // Marker type 29 (big star)
		marker_T->SetMarkerColor(kRed);
		marker_T->SetMarkerSize(3.5);
		marker_T->Draw("same");

		TMarker *marker_B = new TMarker(xCenter_zoomin, yCenter_zoomin + y_cross, 20); // Marker type 29 (big star)
		marker_B->SetMarkerColor(kRed);
		marker_B->SetMarkerSize(3.5);
		marker_B->Draw("same");

		TMarker *marker_C = new TMarker(xCenter_zoomin, yCenter_zoomin, 20); // Marker type 29 (big star)
		marker_C->SetMarkerColor(kRed);
		marker_C->SetMarkerSize(3.5);
		marker_C->Draw("same");

		c_2d_chisquare_ndf[cent]->SaveAs(Form(chi2_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		c_2d_chisquare_ndf_zoomin[cent]->SaveAs(Form(chi2_saving_path_zoomin, this->cenlowlimit[cent], this->cenhighlimit[cent]));

		// This is Data and MC

		c_data_mc_raw[cent]->cd();
		c_data_mc_raw[cent]->SetLeftMargin(0.15);
		c_data_mc_raw[cent]->SetRightMargin(0.08);
		c_data_mc_raw[cent]->SetBottomMargin(0.13);
		c_data_mc_raw[cent]->SetTicks(1, 1);
		// c_data_mc_raw[cent]->SetLogy();

		if (type != 5)
		{
			h_data_bksub[cent]->SetTitle(Form(data_mc_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
			h_data_bksub[cent]->SetMarkerColor(kRed);
			h_data_bksub[cent]->SetMarkerSize(1.5);
			h_data_bksub[cent]->SetMarkerStyle(kFullCircle);
			h_data_bksub[cent]->GetYaxis()->SetTitle("Normalized counts");
			h_data_bksub[cent]->GetXaxis()->SetTitle("m_{u^{+}u^{-}} (GeV)");

			h_data_bksub[cent]->GetYaxis()->SetTitleFont(42);	// Times, bold
			h_data_bksub[cent]->GetYaxis()->SetLabelFont(42);	// Times, bold
			h_data_bksub[cent]->GetYaxis()->SetTitleSize(0.05); // Title size
			h_data_bksub[cent]->GetYaxis()->SetLabelSize(0.04); // Label size
			h_data_bksub[cent]->GetXaxis()->SetTitleFont(42);	// Times, bold
			h_data_bksub[cent]->GetXaxis()->SetLabelFont(42);	// Times, bold
			h_data_bksub[cent]->GetXaxis()->SetTitleSize(0.05); // Title size
			h_data_bksub[cent]->GetXaxis()->SetLabelSize(0.04); // Label size

			h_data_bksub[cent]->Draw("P");
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->SetMarkerColor(kGreen);
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->SetMarkerSize(1.5);
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->Draw("P SAME");
		}
		if (type == 5)
		{
			h_data[cent]->SetTitle(Form(data_mc_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
			h_data[cent]->SetMarkerColor(kRed);
			h_data[cent]->SetMarkerSize(1.5);
			h_data[cent]->SetMarkerStyle(kFullCircle);
			h_data[cent]->GetYaxis()->SetTitle("Normalized counts");
			h_data[cent]->GetXaxis()->SetTitle("m_{u^{+}u^{-}} (GeV)");

			h_data[cent]->GetYaxis()->SetTitleFont(42);	  // Times, bold
			h_data[cent]->GetYaxis()->SetLabelFont(42);	  // Times, bold
			h_data[cent]->GetYaxis()->SetTitleSize(0.05); // Title size
			h_data[cent]->GetYaxis()->SetLabelSize(0.04); // Label size
			h_data[cent]->GetXaxis()->SetTitleFont(42);	  // Times, bold
			h_data[cent]->GetXaxis()->SetLabelFont(42);	  // Times, bold
			h_data[cent]->GetXaxis()->SetTitleSize(0.05); // Title size
			h_data[cent]->GetXaxis()->SetLabelSize(0.04); // Label size

			h_data[cent]->Draw("P");
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->SetMarkerColor(kGreen);
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->SetMarkerSize(1.5);
			h_mc_signal[minBinX - 1][minBinY - 1][cent]->Draw("P SAME");
		}

		TPaveText *pt1 = new TPaveText(0.15, 0.7, 0.5, 0.8, "NDC");
		pt1->AddText("Red is data");
		pt1->AddText("Green is best template");
		pt1->SetTextSize(0.03);
		pt1->SetTextAlign(22); // Center alignment
		pt1->SetFillStyle(0);  // Make the background transparent
		pt1->SetBorderSize(0); // Remove the border
		pt1->SetLineColor(0);  // Remove the border line (optional)
		pt1->SetTextColor(1);  // Set text color (default: black)
		pt1->Draw();

		c_data_mc_raw[cent]->SaveAs(Form(data_mc_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));

		// This is data vs data_bksub
		c_data_data_bk[cent]->cd();
		c_data_data_bk[cent]->SetLogy();
		h_data[cent]->SetTitle(Form(data_data_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		h_data[cent]->SetMarkerColor(kRed);
		h_data[cent]->SetMarkerStyle(kFullCircle);
		h_data[cent]->Draw("P");
		h_data[cent]->Draw("TEXT SAME");
		h_data_bksub[cent]->SetMarkerColor(kGreen);
		h_data_bksub[cent]->SetMarkerStyle(kFullDotLarge);
		h_data_bksub[cent]->Draw("P SAME");
		h_mc_bk[cent]->SetMarkerColor(kBlue);
		// h_mc_bk[cent]->Draw("P SAME");

		TPaveText *pt = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
		pt->AddText("Red is data, green is data - bk, blue is bk");
		pt->SetTextSize(0.04);
		pt->SetTextAlign(22); // Center alignment
		pt->Draw();

		c_data_data_bk[cent]->SaveAs(Form(data_data_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));
	}

	g_HI_dmass = new TGraphErrors(5, xposition_HI, dMass_HI, xposition_err_HI, dMass_Err_HI);
	g_HI_dwidth = new TGraphErrors(5, xposition_HI, dWidth_HI, xposition_err_HI, dWidth_Err_HI);

	TFile *temp = new TFile("All_plots.root", "UPDATE");
	temp->cd();

	if (type == 1)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_nominal", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_nominal", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_nominal", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_nominal", 2);
		}
	}
	if (type == 2)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_tnpU", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_tnpU", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_tnpU", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_tnpU", 2);
		}
	}
	if (type == 3)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_tnpD", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_tnpD", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_tnpD", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_tnpD", 2);
		}
	}
	if (type == 4)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_acooff", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_acooff", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_acooff", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_acooff", 2);
		}
	}
	if (type == 5)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_nominal_no_bk", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_nominal_no_bk", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_nominal_no_bk", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_nominal_no_bk", 2);
		}
	}
	if (type == 6)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_nominal_uniform_rebin", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_nominal_uniform_rebin", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_nominal_uniform_rebin", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_nominal_uniform_rebin", 2);
		}
	}
	if (type == 7)
	{
		if (iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_eta_nominal_mass_range", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_eta_nominal_mass_range", 2);
		}
		if (!iseta)
		{
			g_HI_dmass->Write("HI_dM_chi2_raw_nominal_mass_range", 2);
			g_HI_dwidth->Write("HI_dWidth_chi2_raw_nominal_mass_range", 2);
		}
	}

	temp->Close();
}

void chisquaretest::plottingandformattingpp(int type, bool iseta)
{

	TString chi2_title;
	TString chi2_saving_path;
	TString chi2_saving_path_zoomin;
	TString data_mc_title;
	TString data_data_title;
	TString data_mc_saving_path;
	TString data_data_saving_path;
	TString contour_saving_path;

	// type 1 = nominal
	// type 2 = tnpU
	// type 3 = tnpD
	// type 4 = Acooff
	// type 5 = nominal without bk subtraction

	if (type == 1)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, Nominal, Period: (%i)";
			data_mc_title = "Eta cut, Nominal, Period: (%i)";
			data_data_title = "Eta_Nominal_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/Nominal/Eta_Nominal_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/Nominal/Eta_Nominal_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/Nominal/Eta_Nominal_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_Nominal_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/Eta_Nominal_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, Nominal, Period: (%i)";
			data_mc_title = "WholeAcceptance, Nominal, Period: (%i)";
			data_data_title = "Raw_Nominal_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/Nominal/Raw_Nominal_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/Nominal/Raw_Nominal_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/Nominal/Raw_Nominal_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_Nominal_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_Nominal_%i.png";
		}
	}
	if (type == 2)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, tnpU, Period: (%i)";
			data_mc_title = "Eta, tnpU, Period: (%i)";
			data_data_title = "Eta_tnpU_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/tnpU/Eta_tnpU_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/tnpU/Eta_tnpU_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/tnpU/Eta_tnpU_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_tnpU_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/Eta_tnpU_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, tnpU, Period: (%i)";
			data_mc_title = "WholeAcceptance, tnpU, Period: (%i)";
			data_data_title = "Raw_tnpU_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/tnpU/Raw_tnpU_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/tnpU/Raw_tnpU_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/tnpU/Raw_tnpU_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_tnpU_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_tnpU_%i.png";
		}
	}
	if (type == 3)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, tnpD, Period: (%i)";
			data_mc_title = "Eta cut, tnpD, Period: (%i)";
			data_data_title = "Eta_tnpD_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/tnpD/Eta_tnpD_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/tnpD/Eta_tnpD_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/tnpD/Eta_tnpD_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_tnpD_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/eta_tnpD_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, tnpD, Period: (%i)";
			data_mc_title = "WholeAcceptance, tnpD, Period: (%i)";
			data_data_title = "Raw_tnpD_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/tnpD/Raw_tnpD_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/tnpD/Raw_tnpD_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/tnpD/Raw_tnpD_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_tnpD_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_tnpD_%i.png";
		}
	}
	if (type == 4)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, Acooff, Period: (%i)";
			data_mc_title = "Eta cut, Acooff, Period: (%i)";
			data_data_title = "Eta_Acooff_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/Acooff/Eta_Acooff_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/Acooff/Eta_Acooff_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/Acooff/Eta_Acooff_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_Acooff_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/Eta_Acooff_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, Acooff, Period: (%i)";
			data_mc_title = "WholeAcceptance, Acooff, Period: (%i)";
			data_data_title = "Raw_Acooff_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/Acooff/Raw_Acooff_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/Acooff/Raw_Acooff_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/Acooff/Raw_Acooff_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_Acooff_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_Acooff_%i.png";
		}
	}

	if (type == 5)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, Nominal no bk, Period: (%i)";
			data_mc_title = "Eta cut, Nominal no bk, Period: (%i)";
			data_data_title = "Eta_Nominal_no_bk_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/Nominal_no_bk/Eta_nominal_no_bk_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/Nominal_no_bk/Eta_nominal_no_bk_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/Nominal_no_bk/Eta_nominal_no_bk_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_nominal_no_bk_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/Eta_nominal_no_bk_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, Nominal no bk, Period: (%i)";
			data_mc_title = "|#eta| < 2.4, Nominal no bk, Period: (%i)";
			data_data_title = "Raw_Nominal_no_bk_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/Nominal_no_bk/Raw_nominal_no_bk_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/Nominal_no_bk/Raw_nominal_no_bk_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/Nominal_no_bk/Raw_nominal_no_bk_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_nominal_no_bk_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_nominal_no_bk_%i.png";
		}
	}

	if (type == 6)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, Nominal uniform rebin , Period: (%i)";
			data_mc_title = "Eta cut, Nominal uniform rebin, Period: (%i)";
			data_data_title = "Eta_Nominal_uniform_rebin_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/Nominal_uniform_rebin/Eta_nominal_uniform_rebin_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/Nominal_uniform_rebin/Eta_nominal_uniform_rebin_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/Nominal_uniform_rebin/Eta_nominal_uniform_rebin_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_nominal_uniform_rebin_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/Eta_nominal_uniform_rebin_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, Nominal uniform rebin, Period: (%i)";
			data_mc_title = "|#eta| < 2.4, Nominal uniform rebin, Period: (%i)";
			data_data_title = "Raw_Nominal_uniform_rebin_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/Nominal_uniform_rebin/Raw_nominal_uniform_rebin_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/Nominal_uniform_rebin/Raw_nominal_uniform_rebin_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/Nominal_uniform_rebin/Raw_nominal_uniform_rebin_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_nominal_uniform_rebin_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_nominal_uniform_rebin_%i.png";
		}
	}

	if (type == 7)
	{
		if (iseta)
		{
			chi2_title = "pp, |#eta| < 1.0, Nominal mass range, Period: (%i)";
			data_mc_title = "Eta cut, Nominal mass range, Period: (%i)";
			data_data_title = "Eta_Nominal_mass_range_%i";

			chi2_saving_path = "./chi2pp/chi2plots/eta/Nominal_mass_range/Eta_nominal_mass_range_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/eta/Nominal_mass_range/Eta_nominal_mass_range_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/eta/Nominal_mass_range/Eta_nominal_mass_range_%i.png";
			data_data_saving_path = "./chi2pp/datadata/eta/Eta_nominal_mass_range_%i.png";
			contour_saving_path = "./chi2pp/contour/eta/Eta_nominal_mass_range_%i.png";
		}
		if (!iseta)
		{
			chi2_title = "pp, |#eta| < 2.4, Nominal mass range, Period: (%i)";
			data_mc_title = "|#eta| < 2.4, Nominal mass range, Period: (%i)";
			data_data_title = "Raw_Nominal_mass_range_%i";

			chi2_saving_path = "./chi2pp/chi2plots/raw/Nominal_mass_range/Raw_nominal_mass_range_%i.png";
			chi2_saving_path_zoomin = "./chi2pp/chi2plots/raw/Nominal_mass_range/Raw_nominal_mass_range_zoomin.png";
			data_mc_saving_path = "./chi2pp/datamc/raw/Nominal_mass_range/Raw_nominal_mass_range_%i.png";
			data_data_saving_path = "./chi2pp/datadata/raw/Raw_nominal_mass_range_%i.png";
			contour_saving_path = "./chi2pp/contour/raw/raw_nominal_mass_range_%i.png";
		}
	}

	for (int runperiod = 0; runperiod < 23; runperiod++)
	{
		// Chi2/ndf plot
		c_2d_chisquare_ndf_pp[runperiod]->cd();
		gStyle->SetPalette(kRainBow);
		h_chisquare_pp[runperiod]->SetTitle(Form(chi2_title, runperiod));
		h_chisquare_pp[runperiod]->SetTitleFont(42);
		h_chisquare_pp[runperiod]->Draw("COLZ");
		// h_chisquare_pp[runperiod]->Draw("TEXTSAME");
		h_chisquare_pp[runperiod]->GetXaxis()->SetNdivisions(21, 0, 0);
		h_chisquare_pp[runperiod]->GetYaxis()->SetNdivisions(21, 0, 0);
		h_chisquare_pp[runperiod]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare_pp[runperiod]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare_pp[runperiod]->GetXaxis()->SetTitle("Mass Shifted Amount (GeV)");
		h_chisquare_pp[runperiod]->GetYaxis()->SetTitle("Width Smeared Amount (GeV)");

		for (int j = 1; j <= nbins_mass_shift; j++)
		{
			for (int k = 1; k <= nbins_smear; k++)
			{
				double xlow = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(k);
				double xup = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(k);
				double ylow = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(j);
				double yup = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(j);

				TBox *box = new TBox(xlow, ylow, xup, yup);
				box->SetFillStyle(0);	   // No fill
				box->SetLineColor(kBlack); // Black border
				box->SetLineWidth(1);	   // Border width

				box->Draw("same");
			}
		}

		Int_t minBinX = -1, minBinY = -1;
		Double_t minContent = h_chisquare_pp[runperiod]->GetMaximum();
		for (Int_t binX = 1; binX <= h_chisquare_pp[runperiod]->GetNbinsX(); ++binX)
		{
			for (Int_t binY = 1; binY <= h_chisquare_pp[runperiod]->GetNbinsY(); ++binY)
			{
				Double_t content = h_chisquare_pp[runperiod]->GetBinContent(binX, binY);
				if (content < minContent)
				{
					minContent = content;
					minBinX = binX;
					minBinY = binY;
				}
			}
		}

		// This is to export minimum region
		this->saveChi2Region(h_chisquare_pp[runperiod], minBinX, minBinY, type, iseta, runperiod, 3, true);

		this->getcontour(h_chisquare_pp[runperiod], 2, minBinX, minBinY, minContent, contour_x_left_onesig_HI, contour_x_right_onesig_HI);

		g_HI_contour_1sig_left = new TGraph(nbins_smear, contour_x_left_weighted_onesig, contour_y_HI[10]);
		g_HI_contour_1sig_right = new TGraph(nbins_smear, contour_x_right_weighted_onesig, contour_y_HI[10]);
		g_HI_contour_2sig_left = new TGraph(nbins_smear, contour_x_left_weighted_twosig, contour_y_HI[10]);
		g_HI_contour_2sig_right = new TGraph(nbins_smear, contour_x_right_weighted_twosig, contour_y_HI[10]);

		TString T_contour_saving_path = Form(contour_saving_path, runperiod);
		TString T_chi2_title = Form(chi2_title, runperiod);

		this->drawcontour(g_HI_contour_1sig_left, g_HI_contour_1sig_right, g_HI_contour_2sig_left, g_HI_contour_2sig_right, runperiod, iseta, true, T_chi2_title, T_contour_saving_path);

		Double_t xMin = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(minBinX);
		Double_t xMax = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(minBinX);
		Double_t yMin = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(minBinY);
		Double_t yMax = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(minBinY);

		Double_t xCenter = h_chisquare_pp[runperiod]->GetXaxis()->GetBinCenter(minBinX);
		Double_t yCenter = h_chisquare_pp[runperiod]->GetYaxis()->GetBinCenter(minBinY);

		for (int nbinssmearing = 1; nbinssmearing <= nbins_smear; ++nbinssmearing)
		{
			c_2d_chisquare_ndf_pp[runperiod]->cd();

			Double_t xMin_left_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_left_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMin_left_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_left_after_onesig_HI[nbinssmearing - 1]);

			Double_t xMax_left_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_left_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMax_left_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_left_after_onesig_HI[nbinssmearing - 1]);

			Double_t yMin_left_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_left_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_left_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_left_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			Double_t xMin_right_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_right_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMin_right_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_right_after_onesig_HI[nbinssmearing - 1]);

			Double_t xMax_right_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_right_before_onesig_HI[nbinssmearing - 1]);
			Double_t xMax_right_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_right_after_onesig_HI[nbinssmearing - 1]);

			Double_t yMin_right_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_right_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_right_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_right_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			Double_t xMin_twosig_left_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_left_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMin_twosig_left_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_left_after_twosig_HI[nbinssmearing - 1]);

			Double_t xMax_twosig_left_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_left_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMax_twosig_left_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_left_after_twosig_HI[nbinssmearing - 1]);

			Double_t yMin_twosig_left_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_twosig_left_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_twosig_left_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_twosig_left_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			Double_t xMin_twosig_right_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_right_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMin_twosig_right_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(contour_x_right_after_twosig_HI[nbinssmearing - 1]);

			Double_t xMax_twosig_right_1 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_right_before_twosig_HI[nbinssmearing - 1]);
			Double_t xMax_twosig_right_2 = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(contour_x_right_after_twosig_HI[nbinssmearing - 1]);

			Double_t yMin_twosig_right_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);
			Double_t yMin_twosig_right_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(nbinssmearing);

			Double_t yMax_twosig_right_1 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);
			Double_t yMax_twosig_right_2 = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(nbinssmearing);

			TBox *box_left_1 = new TBox(xMin_left_1, yMin_left_1, xMax_left_1, yMax_left_1);
			box_left_1->SetLineColor(kGreen);
			box_left_1->SetLineWidth(4);
			box_left_1->SetFillStyle(0);
			box_left_1->Draw("same");

			TBox *box_left_2 = new TBox(xMin_left_2, yMin_left_2, xMax_left_2, yMax_left_2);
			box_left_2->SetLineColor(kGreen);
			box_left_2->SetLineWidth(4);
			box_left_2->SetFillStyle(0);
			box_left_2->Draw("same");

			TBox *box_right_1 = new TBox(xMin_right_1, yMin_right_1, xMax_right_1, yMax_right_1);
			box_right_1->SetLineColor(kOrange);
			box_right_1->SetLineWidth(4);
			box_right_1->SetFillStyle(0);
			box_right_1->Draw("same");

			TBox *box_right_2 = new TBox(xMin_right_2, yMin_right_2, xMax_right_2, yMax_right_2);
			box_right_2->SetLineColor(kOrange);
			box_right_2->SetLineWidth(4);
			box_right_2->SetFillStyle(0);
			box_right_2->Draw("same");

			TBox *box_twosig_left_1 = new TBox(xMin_twosig_left_1, yMin_twosig_left_1, xMax_twosig_left_1, yMax_twosig_left_1);
			box_twosig_left_1->SetLineColor(kMagenta);
			box_twosig_left_1->SetLineWidth(4);
			box_twosig_left_1->SetFillStyle(0);
			box_twosig_left_1->Draw("same");

			TBox *box_twosig_left_2 = new TBox(xMin_twosig_left_2, yMin_twosig_left_2, xMax_twosig_left_2, yMax_twosig_left_2);
			box_twosig_left_2->SetLineColor(kMagenta);
			box_twosig_left_2->SetLineWidth(4);
			box_twosig_left_2->SetFillStyle(0);
			box_twosig_left_2->Draw("same");

			TBox *box_twosig_right_1 = new TBox(xMin_twosig_right_1, yMin_twosig_right_1, xMax_twosig_right_1, yMax_twosig_right_1);
			box_twosig_right_1->SetLineColor(kCyan);
			box_twosig_right_1->SetLineWidth(4);
			box_twosig_right_1->SetFillStyle(0);
			box_twosig_right_1->Draw("same");

			TBox *box_twosig_right_2 = new TBox(xMin_twosig_right_2, yMin_twosig_right_2, xMax_twosig_right_2, yMax_twosig_right_2);
			box_twosig_right_2->SetLineColor(kCyan);
			box_twosig_right_2->SetLineWidth(4);
			box_twosig_right_2->SetFillStyle(0);
			box_twosig_right_2->Draw("same");
		}

		double x_cross = 0;
		double y_cross = 0;

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
		box1->SetLineColor(kRed);
		box1->SetLineWidth(4);
		box1->SetFillStyle(0);
		box1->Draw("same");

		if (runperiod == 22)
		{
			c_2d_chisquare_ndf_pp_zoomin->cd();
			gStyle->SetPalette(kRainBow);
			h_chisquare_pp_zoomin->SetTitle(Form(chi2_title, runperiod));
			h_chisquare_pp_zoomin->SetTitleFont(42);
			h_chisquare_pp_zoomin->Draw("COLZ");
			// h_chisquare_pp[runperiod]->Draw("TEXTSAME");
			h_chisquare_pp_zoomin->GetXaxis()->SetNdivisions(21, 0, 0);
			h_chisquare_pp_zoomin->GetYaxis()->SetNdivisions(21, 0, 0);
			h_chisquare_pp_zoomin->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
			h_chisquare_pp_zoomin->GetYaxis()->SetLabelSize(0.02);
			h_chisquare_pp_zoomin->GetXaxis()->SetTitle("Mass Shifted Amount (GeV)");
			h_chisquare_pp_zoomin->GetYaxis()->SetTitle("Width Smeared Amount (GeV)");

			for (int j = 1; j <= nbins_mass_shift; j++)
			{
				for (int k = 1; k <= nbins_smear; k++)
				{
					double xlow = h_chisquare_pp_zoomin->GetXaxis()->GetBinLowEdge(k);
					double xup = h_chisquare_pp_zoomin->GetXaxis()->GetBinUpEdge(k);
					double ylow = h_chisquare_pp_zoomin->GetYaxis()->GetBinLowEdge(j);
					double yup = h_chisquare_pp_zoomin->GetYaxis()->GetBinUpEdge(j);

					TBox *box = new TBox(xlow, ylow, xup, yup);
					box->SetFillStyle(0);	   // No fill
					box->SetLineColor(kBlack); // Black border
					box->SetLineWidth(1);	   // Border width

					box->Draw("same");
				}
			}

			Int_t minBinX_zoomin = -1, minBinY_zoomin = -1;
			Double_t minContent_zoomin = h_chisquare_pp_zoomin->GetMaximum();

			for (Int_t binX = 1; binX <= h_chisquare_pp_zoomin->GetNbinsX(); ++binX)
			{
				for (Int_t binY = 1; binY <= h_chisquare_pp_zoomin->GetNbinsY(); ++binY)
				{
					Double_t content = h_chisquare_pp_zoomin->GetBinContent(binX, binY);
					if (content < minContent_zoomin)
					{
						minContent_zoomin = content;
						minBinX_zoomin = binX;
						minBinY_zoomin = binY;
					}
				}
			}

			Double_t xMin_zoomin = h_chisquare_pp_zoomin->GetXaxis()->GetBinLowEdge(minBinX_zoomin);
			Double_t xMax_zoomin = h_chisquare_pp_zoomin->GetXaxis()->GetBinUpEdge(minBinX_zoomin);
			Double_t yMin_zoomin = h_chisquare_pp_zoomin->GetYaxis()->GetBinLowEdge(minBinY_zoomin);
			Double_t yMax_zoomin = h_chisquare_pp_zoomin->GetYaxis()->GetBinUpEdge(minBinY_zoomin);

			Double_t xCenter_zoomin = h_chisquare_pp_zoomin->GetXaxis()->GetBinCenter(minBinX_zoomin);
			Double_t yCenter_zoomin = h_chisquare_pp_zoomin->GetYaxis()->GetBinCenter(minBinY_zoomin);

			TBox *box1_zoomin = new TBox(xMin_zoomin, yMin_zoomin, xMax_zoomin, yMax_zoomin);
			box1_zoomin->SetLineColor(kRed);
			box1_zoomin->SetLineWidth(4);
			box1_zoomin->SetFillStyle(0);
			box1_zoomin->Draw("same");

			dMass_pp[22] = xCenter_zoomin;
			dWidth_pp[22] = yCenter_zoomin;

			// cout << "dWidth is " << dWidth_HI[cent + 1] << endl;

			dMass_Err_pp[22] = getuncertainty(h_chisquare_pp[22], 1, minBinX, minBinY);
			dWidth_Err_pp[22] = getuncertainty(h_chisquare_pp[22], 2, minBinX, minBinY);

			double diffindM = xCenter_zoomin - xCenter;
			dMass_Err_pp[22] = dMass_Err_pp[22] + diffindM;

			double diffindW = yCenter_zoomin - yCenter;
			dWidth_Err_pp[22] = dWidth_Err_pp[22] + diffindW;

			x_cross = dMass_Err_pp[22];
			y_cross = dWidth_Err_pp[22];

			c_2d_chisquare_ndf_pp[22]->cd();

			TMarker *marker_L = new TMarker(xCenter_zoomin - x_cross, yCenter_zoomin, 20); // Marker type 29 (big star)
			marker_L->SetMarkerColor(kRed);
			marker_L->SetMarkerSize(3.5);
			marker_L->Draw("same");

			TMarker *marker_R = new TMarker(xCenter_zoomin + x_cross, yCenter_zoomin, 20); // Marker type 29 (big star)
			marker_R->SetMarkerColor(kRed);
			marker_R->SetMarkerSize(3.5);
			marker_R->Draw("same");

			TMarker *marker_T = new TMarker(xCenter_zoomin, yCenter_zoomin - y_cross, 20); // Marker type 29 (big star)
			marker_T->SetMarkerColor(kRed);
			marker_T->SetMarkerSize(3.5);
			marker_T->Draw("same");

			TMarker *marker_B = new TMarker(xCenter_zoomin, yCenter_zoomin + y_cross, 20); // Marker type 29 (big star)
			marker_B->SetMarkerColor(kRed);
			marker_B->SetMarkerSize(3.5);
			marker_B->Draw("same");

			TMarker *marker_C = new TMarker(xCenter_zoomin, yCenter_zoomin, 20); // Marker type 29 (big star)
			marker_C->SetMarkerColor(kRed);
			marker_C->SetMarkerSize(3.5);
			marker_C->Draw("same");
		}
		else
		{

			dMass_pp[runperiod] = xCenter;
			dWidth_pp[runperiod] = yCenter;

			dMass_Err_pp[runperiod] = getuncertainty(h_chisquare_pp[runperiod], 1, minBinX, minBinY);
			dWidth_Err_pp[runperiod] = getuncertainty(h_chisquare_pp[runperiod], 2, minBinX, minBinY);

			x_cross = dMass_Err_pp[runperiod];
			y_cross = dWidth_Err_pp[runperiod];

			c_2d_chisquare_ndf_pp[runperiod]->cd();

			TMarker *marker_L = new TMarker(h_chisquare_pp[runperiod]->GetXaxis()->GetBinCenter(minBinX) - x_cross, h_chisquare_pp[runperiod]->GetYaxis()->GetBinCenter(minBinY), 29); // Marker type 29 (big star)
			marker_L->SetMarkerColor(kRed);
			marker_L->SetMarkerSize(5.5);
			marker_L->Draw("same");

			TMarker *marker_R = new TMarker(h_chisquare_pp[runperiod]->GetXaxis()->GetBinCenter(minBinX) + x_cross, h_chisquare_pp[runperiod]->GetYaxis()->GetBinCenter(minBinY), 29); // Marker type 29 (big star)
			marker_R->SetMarkerColor(kRed);
			marker_R->SetMarkerSize(5.5);
			marker_R->Draw("same");

			TMarker *marker_T = new TMarker(h_chisquare_pp[runperiod]->GetXaxis()->GetBinCenter(minBinX), h_chisquare_pp[runperiod]->GetYaxis()->GetBinCenter(minBinY) - y_cross, 29); // Marker type 29 (big star)
			marker_T->SetMarkerColor(kRed);
			marker_T->SetMarkerSize(5.5);
			marker_T->Draw("same");

			TMarker *marker_B = new TMarker(h_chisquare_pp[runperiod]->GetXaxis()->GetBinCenter(minBinX), h_chisquare_pp[runperiod]->GetYaxis()->GetBinCenter(minBinY) + y_cross, 29); // Marker type 29 (big star)
			marker_B->SetMarkerColor(kRed);
			marker_B->SetMarkerSize(5.5);
			marker_B->Draw("same");
		}

		c_2d_chisquare_ndf_pp[runperiod]->SaveAs(Form(chi2_saving_path, runperiod));
		if (runperiod == 22)
		{
			c_2d_chisquare_ndf_pp_zoomin->SaveAs(chi2_saving_path_zoomin);
		}

		// This is Data and MC

		c_data_mc_raw_pp[runperiod]->cd();
		c_data_mc_raw_pp[runperiod]->SetLeftMargin(0.15);
		c_data_mc_raw_pp[runperiod]->SetRightMargin(0.08);
		c_data_mc_raw_pp[runperiod]->SetBottomMargin(0.13);
		c_data_mc_raw_pp[runperiod]->SetTicks(1, 1);
		// c_data_mc_raw_pp[runperiod]->SetLogy();
		if (type != 5)
		{
			h_data_bksub_pp[runperiod]->SetTitle(Form(data_mc_title, runperiod));
			h_data_bksub_pp[runperiod]->SetMarkerColor(kRed);
			h_data_bksub_pp[runperiod]->SetMarkerSize(1.5);
			h_data_bksub_pp[runperiod]->SetMarkerStyle(kFullCircle);
			h_data_bksub_pp[runperiod]->GetYaxis()->SetTitle("Normalized counts");
			h_data_bksub_pp[runperiod]->GetXaxis()->SetTitle("m_{u^{+}u^{-}} (GeV)");

			h_data_bksub_pp[runperiod]->GetYaxis()->SetTitleFont(42);	// Times, bold
			h_data_bksub_pp[runperiod]->GetYaxis()->SetLabelFont(42);	// Times, bold
			h_data_bksub_pp[runperiod]->GetYaxis()->SetTitleSize(0.05); // Title size
			h_data_bksub_pp[runperiod]->GetYaxis()->SetLabelSize(0.04); // Label size
			h_data_bksub_pp[runperiod]->GetXaxis()->SetTitleFont(42);	// Times, bold
			h_data_bksub_pp[runperiod]->GetXaxis()->SetLabelFont(42);	// Times, bold
			h_data_bksub_pp[runperiod]->GetXaxis()->SetTitleSize(0.05); // Title size
			h_data_bksub_pp[runperiod]->GetXaxis()->SetLabelSize(0.04); // Label size

			h_data_bksub_pp[runperiod]->Draw("P");

			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerColor(kGreen);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerSize(1.5);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->Draw("P SAME");
		}
		if (type == 5)
		{
			h_data_pp[runperiod]->SetTitle(Form(data_mc_title, runperiod));
			h_data_pp[runperiod]->SetMarkerColor(kRed);
			h_data_pp[runperiod]->SetMarkerSize(1.5);
			h_data_pp[runperiod]->SetMarkerStyle(kFullCircle);
			h_data_pp[runperiod]->GetYaxis()->SetTitle("Normalized counts");
			h_data_pp[runperiod]->GetXaxis()->SetTitle("m_{u^{+}u^{-}} (GeV)");

			h_data_pp[runperiod]->GetYaxis()->SetTitleFont(42);	  // Times, bold
			h_data_pp[runperiod]->GetYaxis()->SetLabelFont(42);	  // Times, bold
			h_data_pp[runperiod]->GetYaxis()->SetTitleSize(0.05); // Title size
			h_data_pp[runperiod]->GetYaxis()->SetLabelSize(0.04); // Label size
			h_data_pp[runperiod]->GetXaxis()->SetTitleFont(42);	  // Times, bold
			h_data_pp[runperiod]->GetXaxis()->SetLabelFont(42);	  // Times, bold
			h_data_pp[runperiod]->GetXaxis()->SetTitleSize(0.05); // Title size
			h_data_pp[runperiod]->GetXaxis()->SetLabelSize(0.04); // Label size

			h_data_pp[runperiod]->Draw("P");
			// h_data_bksub_pp[runperiod]->SetMarkerColor(kBlue);
			// h_data_bksub_pp[runperiod]->SetMarkerStyle(kFullCircle);
			// h_data_bksub_pp[runperiod]->Draw("PSAME");
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerColor(kGreen);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerSize(1.5);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->Draw("P SAME");
		}

		TPaveText *pt1 = new TPaveText(0.15, 0.7, 0.5, 0.8, "NDC");
		pt1->AddText("Red is data");
		pt1->AddText("Green is best template");
		pt1->SetTextSize(0.03);
		pt1->SetTextAlign(22); // Center alignment
		pt1->SetFillStyle(0);  // Make the background transparent
		pt1->SetBorderSize(0); // Remove the border
		pt1->SetLineColor(0);  // Remove the border line (optional)
		pt1->SetTextColor(1);  // Set text color (default: black)
		pt1->Draw();

		c_data_mc_raw_pp[runperiod]->SaveAs(Form(data_mc_saving_path, runperiod));

		// This is data vs data_bksub
		c_data_data_bk_pp[runperiod]->cd();
		c_data_data_bk_pp[runperiod]->SetLogy();
		h_data_pp[runperiod]->SetTitle(Form(data_data_title, runperiod));
		h_data_pp[runperiod]->SetMarkerColor(kRed);
		h_data_pp[runperiod]->SetMarkerStyle(kFullCircle);
		h_data_pp[runperiod]->Draw("P");
		h_data_pp[runperiod]->Draw("TEXT SAME");
		h_data_bksub_pp[runperiod]->SetMarkerStyle(kFullDotLarge);
		h_data_bksub_pp[runperiod]->SetMarkerColor(kGreen);
		h_data_bksub_pp[runperiod]->Draw("P SAME");

		TPaveText *pt = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
		pt->AddText("Red is data, green is data - bk, blue is bk");
		pt->SetTextSize(0.04);
		pt->SetTextAlign(22); // Center alignment
		pt->Draw();

		c_data_data_bk_pp[runperiod]->SaveAs(Form(data_data_saving_path, runperiod));
	}

	g_pp_dmass = new TGraphErrors(22, xposition, dMass_pp, xposition_err, dMass_Err_pp);
	g_pp_dwidth = new TGraphErrors(22, xposition, dWidth_pp, xposition_err, dWidth_Err_pp);

	TFile *temp = new TFile("All_plots.root", "UPDATE");
	temp->cd();

	if (type == 1)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_nominal", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_nominal", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_nominal", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_nominal", 2);
		}
	}
	if (type == 2)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_tnpU", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_tnpU", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_tnpU", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_tnpU", 2);
		}
	}
	if (type == 3)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_tnpD", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_tnpD", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_tnpD", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_tnpD", 2);
		}
	}
	if (type == 4)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_acooff", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_acooff", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_acooff", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_acooff", 2);
		}
	}
	if (type == 5)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_nominal_no_bk", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_nominal_no_bk", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_nominal_no_bk", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_nominal_no_bk", 2);
		}
	}
	if (type == 6)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_nominal_uniform_rebin", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_nominal_uniform_rebin", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_nominal_uniform_rebin", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_nominal_uniform_rebin", 2);
		}
	}
	if (type == 7)
	{
		if (iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_eta_nominal_mass_range", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_eta_nominal_mass_range", 2);
		}
		if (!iseta)
		{
			g_pp_dmass->Write("pp_dM_chi2_raw_nominal_mass_range", 2);
			g_pp_dwidth->Write("pp_dWidth_chi2_raw_nominal_mass_range", 2);
		}
	}

	temp->Close();
}

void chisquaretest::bincontentcheck(bool isbk)
{
	int nbins;

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		if (isbk)
			nbins = h_data[cent]->GetNbinsX();
		if (!isbk)
			nbins = h_data_bksub[cent]->GetNbinsX();

		for (int bin = 1; bin <= nbins; bin++)
		{
			if (!isbk)
			{
				if (h_data[cent]->GetBinContent(bin) <= 0)
					cout << "Warning: " << "Cent " << this->cenlowlimit[cent] << " " << this->cenhighlimit[cent] << " Bin: " << bin << " is <= 0 !" << endl;
			}
			if (isbk)
			{
				if (h_data_bksub[cent]->GetBinContent(bin) <= 0)
					cout << "Warning(bksub): " << "Cent " << this->cenlowlimit[cent] << " " << this->cenhighlimit[cent] << " Bin: " << bin << " is <= 0 !" << endl;
			}
		}
	}
}

void chisquaretest::RebinAll(bool iseta, int type)
{

	std::vector<double> newBinEdges_0_100 = createCustomBinning(h_data[0], 60, 80, 10, 80, 86, 4, 95, 100, 5, 100, 120, 40);
	std::vector<double> newBinEdges_mass_range = createCustomBinning(h_data[0], 70, 80, 5, 80, 86, 4, 95, 100, 5, 100, 110, 20);

	int nNewBins = newBinEdges_0_100.size() - 1;
	int nNewBins_mass_range = newBinEdges_mass_range.size() - 1;

	Double_t *binEdgesArray = &newBinEdges_0_100[0];
	Double_t *binEdgesArray_mass_range = &newBinEdges_mass_range[0];

	if (iseta)
	{
		if (type == 6)
		{
			for (int cent = 0; cent < nbins_cent; cent++)
			{
				if (!((cent < 4) || (cent == 10)))
					continue;

				h_data[cent]->Rebin(4);
				h_data_bksub[cent]->Rebin(4);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{

						h_mc_signal[shift][smear][cent]->Rebin(4);
						h_mc_signal_zoomin[shift][smear][cent]->Rebin(4);
					}
				}
			}
		}
		else if (type == 7)
		{
			for (int cent = 0; cent < nbins_cent; cent++)
			{
				if (!((cent < 4) || (cent == 10)))
					continue;

				h_data[cent] = (TH1D *)h_data[cent]->Rebin(nNewBins_mass_range, Form("mass_array_data_witheta_rebin_%i", cent), binEdgesArray_mass_range);
				h_data_bksub[cent] = (TH1D *)h_data_bksub[cent]->Rebin(nNewBins_mass_range, Form("normalized_mc_bk_rebinned_%i", cent), binEdgesArray_mass_range);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)h_mc_signal[shift][smear][cent]->Rebin(nNewBins_mass_range, Form("modifiedmass_eta_without_eff_%i_%i_%i_new", shift, smear, cent), binEdgesArray_mass_range);
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)h_mc_signal_zoomin[shift][smear][cent]->Rebin(nNewBins_mass_range, Form("modifiedmass_eta_without_eff_%i_%i_%i_zoomin", shift, smear, cent), binEdgesArray_mass_range);
					}
				}
			}
		}
		else
		{
			for (int cent = 0; cent < nbins_cent; cent++)
			{
				if (!((cent < 4) || (cent == 10)))
					continue;

				h_data[cent] = (TH1D *)h_data[cent]->Rebin(nNewBins, Form("mass_array_data_witheta_rebin_%i", cent), binEdgesArray);
				h_data_bksub[cent] = (TH1D *)h_data_bksub[cent]->Rebin(nNewBins, Form("normalized_mc_bk_rebinned_%i", cent), binEdgesArray);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)h_mc_signal[shift][smear][cent]->Rebin(nNewBins, Form("modifiedmass_eta_without_eff_%i_%i_%i_new", shift, smear, cent), binEdgesArray);
						h_mc_signal_zoomin[shift][smear][cent] = (TH1D *)h_mc_signal_zoomin[shift][smear][cent]->Rebin(nNewBins, Form("modifiedmass_eta_without_eff_%i_%i_%i_zoomin", shift, smear, cent), binEdgesArray);
					}
				}
			}
		}
	}

	if (!iseta)
	{
		if (type == 6)
		{
			// Here is empty since we want to have no rebin
		}
		else
		{
			for (int cent = 0; cent < nbins_cent; cent++)
			{
				if (!((cent < 4) || (cent == 10)))
					continue;

				h_data[cent]->Rebin(4);
				h_data_bksub[cent]->Rebin(4);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{
						h_mc_signal[shift][smear][cent]->Rebin(4);
						h_mc_signal_zoomin[shift][smear][cent]->Rebin(4);
					}
				}
			}
		}
	}
}

void chisquaretest::RebinAllpp(bool iseta, int type)
{
	std::vector<double> newBinEdges_0_100 = createCustomBinning(h_data_pp[0], 60, 80, 10, 80, 86, 4, 95, 100, 5, 100, 120, 40);
	std::vector<double> newBinEdges_mass_range = createCustomBinning(h_data_pp[0], 70, 80, 5, 80, 86, 4, 95, 100, 5, 100, 110, 20);

	int nNewBins = newBinEdges_0_100.size() - 1;
	int nNewBins_mass_range = newBinEdges_mass_range.size() - 1;

	Double_t *binEdgesArray = &newBinEdges_0_100[0];
	Double_t *binEdgesArray_mass_range = &newBinEdges_mass_range[0];

	if (iseta)
	{
		if (type == 6)
		{
			for (int runperiod = 0; runperiod < 23; runperiod++)
			{

				h_data_pp[runperiod]->Rebin(4);
				h_data_bksub_pp[runperiod]->Rebin(4);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{
						if (runperiod == 0)
						{
							h_mc_signal_pp_zoomin[shift][smear]->Rebin(4);
							h_mc_signal_pp[shift][smear][runperiod]->Rebin(4);
						}
					}
				}
			}
		}
		else if (type == 7)
		{
			for (int runperiod = 0; runperiod < 23; runperiod++)
			{

				h_data_pp[runperiod] = (TH1D *)h_data_pp[runperiod]->Rebin(nNewBins_mass_range, Form("mass_array_data_witheta_rebin_%i", runperiod), binEdgesArray_mass_range);
				h_data_bksub_pp[runperiod] = (TH1D *)h_data_bksub_pp[runperiod]->Rebin(nNewBins_mass_range, Form("normalized_mc_bk_rebinned_%i", runperiod), binEdgesArray_mass_range);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)h_mc_signal_pp[shift][smear][runperiod]->Rebin(nNewBins_mass_range, Form("modifiedmass_eta_without_eff_%i_%i_%i_new", shift, smear, runperiod), binEdgesArray_mass_range);
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)h_mc_signal_pp_zoomin[shift][smear]->Rebin(nNewBins_mass_range, Form("modifiedmass_eta_without_eff_%i_%i_%i_zoomin", shift, smear, runperiod), binEdgesArray_mass_range);
						}
					}
				}
			}
		}
		else
		{
			for (int runperiod = 0; runperiod < 23; runperiod++)
			{

				h_data_pp[runperiod] = (TH1D *)h_data_pp[runperiod]->Rebin(nNewBins, Form("mass_array_data_witheta_rebin_%i", runperiod), binEdgesArray);
				h_data_bksub_pp[runperiod] = (TH1D *)h_data_bksub_pp[runperiod]->Rebin(nNewBins, Form("normalized_mc_bk_rebinned_%i", runperiod), binEdgesArray);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{

						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)h_mc_signal_pp[shift][smear][runperiod]->Rebin(nNewBins, Form("modifiedmass_eta_without_eff_%i_%i_%i_new", shift, smear, runperiod), binEdgesArray);
						if (runperiod == 22)
						{
							h_mc_signal_pp_zoomin[shift][smear] = (TH1D *)h_mc_signal_pp_zoomin[shift][smear]->Rebin(nNewBins, Form("modifiedmass_eta_without_eff_%i_%i_%i_zoomin", shift, smear, runperiod), binEdgesArray);
						}
					}
				}
			}
		}
	}

	if (!iseta)
	{
		if (type == 6)
		{
			// Here is empty since we want to have no rebin
		}
		else
		{
			for (int runperiod = 0; runperiod < 23; runperiod++)
			{
				h_data_pp[runperiod]->Rebin(4);
				h_data_bksub_pp[runperiod]->Rebin(4);

				for (int shift = 0; shift < nbins_mass_shift; shift++)
				{
					for (int smear = 0; smear < nbins_smear; smear++)
					{

						if (runperiod == 0)
						{
							h_mc_signal_pp[shift][smear][runperiod]->Rebin(4);
							h_mc_signal_pp_zoomin[shift][smear]->Rebin(4);
						}
					}
				}
			}
		}
	}
}

Double_t chisquaretest::getuncertainty(TH2D *h_1, int type, Int_t minBinX, Int_t minBinY)
{
	double threshold = 2.3;

	if (type == 1)
	{
		// Get the chi² minimum value
		double chi2_min = h_1->GetBinContent(minBinX, minBinY);
		double chi2_target = chi2_min + threshold;

		int bin_before_L = -1, bin_after_L = -1;
		int bin_before_R = -1, bin_after_R = -1;
		double chi2_before_L = 0, chi2_after_L = 0;
		double chi2_before_R = 0, chi2_after_R = 0;

		// Scan to the RIGHT of minBinX (increasing X)
		for (int binx = minBinX + 1; binx <= h_1->GetNbinsX(); ++binx)
		{
			double chi2_bin = h_1->GetBinContent(binx, minBinY);

			if (chi2_bin > chi2_target)
			{
				bin_before_R = binx - 1;
				bin_after_R = binx;
				chi2_before_R = h_1->GetBinContent(bin_before_R, minBinY);
				chi2_after_R = chi2_bin;
				break;
			}
		}

		// Scan to the LEFT of minBinX (decreasing X)
		for (int binx = minBinX - 1; binx >= 1; --binx)
		{
			double chi2_bin = h_1->GetBinContent(binx, minBinY);

			if (chi2_bin > chi2_target)
			{
				bin_before_L = binx;
				bin_after_L = binx + 1;
				chi2_before_L = chi2_bin;
				chi2_after_L = h_1->GetBinContent(bin_after_L, minBinY);
				break;
			}
		}

		double x_cross_L = -999, x_cross_R = -999;
		bool found_L = (bin_before_L != -1 && bin_after_L != -1);
		bool found_R = (bin_before_R != -1 && bin_after_R != -1);

		// Linear interpolation if crossing exists on the LEFT
		if (found_L)
		{
			double x_before_L = h_1->GetXaxis()->GetBinCenter(bin_before_L);
			double x_after_L = h_1->GetXaxis()->GetBinCenter(bin_after_L);
			x_cross_L = (chi2_before_L * x_before_L + chi2_after_L * x_after_L) / (chi2_before_L + chi2_after_L);
		}

		// Linear interpolation if crossing exists on the RIGHT
		if (found_R)
		{
			double x_before_R = h_1->GetXaxis()->GetBinCenter(bin_before_R);
			double x_after_R = h_1->GetXaxis()->GetBinCenter(bin_after_R);
			x_cross_R = (chi2_before_R * x_before_R + chi2_after_R * x_after_R) / (chi2_before_R + chi2_after_R);
		}

		// Determine which crossing to use
		if (found_L && found_R)
		{
			double dist_L = fabs(x_cross_L - h_1->GetXaxis()->GetBinCenter(minBinX));
			double dist_R = fabs(x_cross_R - h_1->GetXaxis()->GetBinCenter(minBinX));

			if (dist_L > dist_R)
			{
				std::cout << "Using LEFT threshold crossing at x = " << x_cross_L << std::endl;
				return dist_L;
			}
			else
			{
				std::cout << "Using RIGHT threshold crossing at x = " << x_cross_R << std::endl;
				return dist_R;
			}
		}
		else if (found_L)
		{
			std::cout << "Using LEFT threshold crossing at x = " << x_cross_L << std::endl;
			return fabs(x_cross_L - h_1->GetXaxis()->GetBinCenter(minBinX));
		}
		else if (found_R)
		{
			std::cout << "Using RIGHT threshold crossing at x = " << x_cross_R << std::endl;
			return fabs(x_cross_R - h_1->GetXaxis()->GetBinCenter(minBinX));
		}
		else
		{
			std::cout << "No threshold crossing found!" << std::endl;
			return -99;
		}
	}

	else if (type == 2)
	{
		// This is for dW (Y-axis uncertainty)
		double threshold = 2.3; // Set the threshold above the minimum chi²

		// Get the chi² minimum value
		double chi2_min = h_1->GetBinContent(minBinX, minBinY);
		double chi2_target = chi2_min + threshold;

		int bin_before_B = -1, bin_after_B = -1;
		int bin_before_T = -1, bin_after_T = -1;
		double chi2_before_B = 0, chi2_after_B = 0;
		double chi2_before_T = 0, chi2_after_T = 0;

		// Scan DOWN (decreasing Y, towards bottom)
		for (int biny = minBinY - 1; biny >= 1; --biny)
		{
			double chi2_bin = h_1->GetBinContent(minBinX, biny);

			if (chi2_bin > chi2_target)
			{
				bin_before_B = biny + 1;
				bin_after_B = biny;
				chi2_before_B = h_1->GetBinContent(bin_before_B, minBinX);
				chi2_after_B = chi2_bin;
				break;
			}
		}

		// Scan UP (increasing Y, towards top)
		for (int biny = minBinY + 1; biny <= h_1->GetNbinsY(); ++biny)
		{
			double chi2_bin = h_1->GetBinContent(minBinX, biny);

			if (chi2_bin > chi2_target)
			{
				bin_before_T = biny - 1;
				bin_after_T = biny;
				chi2_before_T = h_1->GetBinContent(bin_before_T, minBinX);
				chi2_after_T = chi2_bin;
				break;
			}
		}

		double y_cross_B = -999, y_cross_T = -999;
		bool found_B = (bin_before_B != -1 && bin_after_B != -1);
		bool found_T = (bin_before_T != -1 && bin_after_T != -1);

		// Linear interpolation if crossing exists at the BOTTOM
		if (found_B)
		{
			double y_before_B = h_1->GetYaxis()->GetBinCenter(bin_before_B);
			double y_after_B = h_1->GetYaxis()->GetBinCenter(bin_after_B);
			y_cross_B = (chi2_before_B * y_before_B + chi2_after_B * y_after_B) / (chi2_before_B + chi2_after_B);
		}

		// Linear interpolation if crossing exists at the TOP
		if (found_T)
		{
			double y_before_T = h_1->GetYaxis()->GetBinCenter(bin_before_T);
			double y_after_T = h_1->GetYaxis()->GetBinCenter(bin_after_T);
			y_cross_T = (chi2_before_T * y_before_T + chi2_after_T * y_after_T) / (chi2_before_T + chi2_after_T);
		}

		// Determine which crossing to use
		if (found_B && found_T)
		{
			double dist_B = fabs(y_cross_B - h_1->GetYaxis()->GetBinCenter(minBinY));
			double dist_T = fabs(y_cross_T - h_1->GetYaxis()->GetBinCenter(minBinY));

			if (dist_B > dist_T)
			{
				std::cout << "Using BOTTOM threshold crossing at y = " << y_cross_B << std::endl;
				return dist_B;
			}
			else
			{
				std::cout << "Using TOP threshold crossing at y = " << y_cross_T << std::endl;
				return dist_T;
			}
		}
		else if (found_B)
		{
			std::cout << "Using BOTTOM threshold crossing at y = " << y_cross_B << std::endl;
			return fabs(y_cross_B - h_1->GetYaxis()->GetBinCenter(minBinY));
		}
		else if (found_T)
		{
			std::cout << "Using TOP threshold crossing at y = " << y_cross_T << std::endl;
			return fabs(y_cross_T - h_1->GetYaxis()->GetBinCenter(minBinY));
		}
		else
		{
			std::cout << "No threshold crossing found!" << std::endl;
			return -99;
		}
	}
	else
	{
		cout << "You should never see this" << endl;
		return -99;
	}
}

Double_t chisquaretest::myownfunctionchi2(TH1D *h1, TH1D *h2)
{
	int nbinsh1 = h1->GetNbinsX();
	int nbinsh2 = h2->GetNbinsX();

	if (nbinsh1 != nbinsh2)
	{
		std::cerr << "Warning: h1 and h2 have different bin numbers!" << nbinsh1 << " " << nbinsh2 << std::endl;
		return -99; // Optionally, return an error value
	}

	Double_t chi2 = 0;
	for (int i = 1; i <= nbinsh1; i++)
	{
		double num = TMath::Power((h1->GetBinContent(i) - h2->GetBinContent(i)), 2);
		double den = TMath::Power((h1->GetBinError(i)), 2) + TMath::Power((h2->GetBinError(i)), 2);
		// double den = TMath::Power((h1->GetBinError(i)), 2) + TMath::Power((h2->GetBinError(i)), 2);

		if (TMath::Power((h1->GetBinContent(i)), 1) <= 0)
		{
			std::cerr << "Warning: Bin " << i << " has zero errors/missing bins/ missing bins after bksub for data. Not Skipped" << std::endl;
			// continue; // Skip this bin to avoid division by zero
		}

		chi2 += num / den;
	}

	Double_t ndf = nbinsh1 - 2;

	return chi2;
}

void chisquaretest::getcontour(TH2D *h1, int type, Int_t minBinX, Int_t minBinY, Double_t minBinContent, Double_t *arrayleft, Double_t *arrayright)
{
	localmin_x = h1->GetXaxis()->GetBinCenter(minBinX);
	localmin_y = h1->GetYaxis()->GetBinCenter(minBinY);
	// cleanning up the array
	for (int i = 0; i < nbins_smear; i++)
	{
		contour_x_left_before_onesig_HI[i] = -99;
		contour_x_left_after_onesig_HI[i] = -99;
		contour_x_right_before_onesig_HI[i] = -99;
		contour_x_right_after_onesig_HI[i] = -99;

		contour_x_left_before_twosig_HI[i] = -99;
		contour_x_left_after_twosig_HI[i] = -99;
		contour_x_right_before_twosig_HI[i] = -99;
		contour_x_right_after_twosig_HI[i] = -99;

		contour_x_left_weighted_twosig[i] = -99;
		contour_x_right_weighted_twosig[i] = -99;
		contour_x_left_weighted_onesig[i] = -99;
		contour_x_right_weighted_onesig[i] = -99;
	}
	// Now I am only doing one dim scanning, i.e. fix smearing, shift in mass
	double nbinsx = h1->GetNbinsX();
	double nbinsy = h1->GetNbinsY();
	double threshold = minBinContent + 2.30;
	double threshold_2 = minBinContent + 5.99;

	for (int yindex = 1; yindex <= nbinsy; yindex++)
	{
		bool isleft = 1;
		for (int xindex = 1; xindex <= nbinsx; xindex++)
		{
			double bincontent = h1->GetBinContent(xindex, yindex);
			if (bincontent < threshold && isleft)
			{
				// cout << "The first bin crossing threshold on the left is " << xindex << endl;
				double binleft_1 = h1->GetBinContent(xindex - 1, yindex);
				double binleft_2 = bincontent;
				double binleft_1_center = h1->GetXaxis()->GetBinCenter(xindex - 1);
				double binleft_2_center = h1->GetXaxis()->GetBinCenter(xindex);

				contour_x_left_before_onesig_HI[yindex - 1] = xindex - 1;
				contour_x_left_after_onesig_HI[yindex - 1] = xindex;

				double weighted_left = (binleft_1_center * binleft_1 + binleft_2_center * binleft_2) / (binleft_1 + binleft_2);
				contour_x_left_weighted_onesig[yindex - 1] = weighted_left;
				// arrayleft[yindex - 1] = weighted_left;

				// cout << "Here's the weighted left" << weighted_left << " For smearing " << yindex << endl;

				isleft = 0;
			}
			if (bincontent > threshold && !isleft)
			{
				// cout << "The first bin crossing threshold on the right is " << xindex << endl;
				double binright_1 = h1->GetBinContent(xindex - 1, yindex);
				double binright_2 = bincontent;
				double binright_1_center = h1->GetXaxis()->GetBinCenter(xindex - 1);
				double binright_2_center = h1->GetXaxis()->GetBinCenter(xindex);

				contour_x_right_before_onesig_HI[yindex - 1] = xindex - 1;
				contour_x_right_after_onesig_HI[yindex - 1] = xindex;

				double weighted_right = (binright_1_center * binright_1 + binright_2_center * binright_2) / (binright_1 + binright_2);
				contour_x_right_weighted_onesig[yindex - 1] = weighted_right;
				// cout << "Here's the weighted right" << weighted_right << " For smearing " << yindex << endl;
				// arrayright[yindex - 1] = weighted_right;
				break;
			}
		}
	}

	for (int yindex = 1; yindex <= nbinsy; yindex++)
	{
		bool isleft = 1;
		for (int xindex = 1; xindex <= nbinsx; xindex++)
		{
			double bincontent = h1->GetBinContent(xindex, yindex);
			if (bincontent < threshold_2 && isleft)
			{
				// cout << "The first bin crossing threshold on the left is " << xindex << endl;
				double binleft_1 = h1->GetBinContent(xindex - 1, yindex);
				double binleft_2 = bincontent;
				double binleft_1_center = h1->GetXaxis()->GetBinCenter(xindex - 1);
				double binleft_2_center = h1->GetXaxis()->GetBinCenter(xindex);

				contour_x_left_before_twosig_HI[yindex - 1] = xindex - 1;
				contour_x_left_after_twosig_HI[yindex - 1] = xindex;

				double weighted_left = (binleft_1_center * binleft_1 + binleft_2_center * binleft_2) / (binleft_1 + binleft_2);
				contour_x_left_weighted_twosig[yindex - 1] = weighted_left;
				// arrayleft[yindex - 1] = weighted_left;

				// cout << "Here's the weighted left" << weighted_left << " For smearing " << yindex << endl;

				isleft = 0;
			}
			if (bincontent > threshold_2 && !isleft)
			{
				// cout << "The first bin crossing threshold on the right is " << xindex << endl;
				double binright_1 = h1->GetBinContent(xindex - 1, yindex);
				double binright_2 = bincontent;
				double binright_1_center = h1->GetXaxis()->GetBinCenter(xindex - 1);
				double binright_2_center = h1->GetXaxis()->GetBinCenter(xindex);

				contour_x_right_before_twosig_HI[yindex - 1] = xindex - 1;
				contour_x_right_after_twosig_HI[yindex - 1] = xindex;

				double weighted_right = (binright_1_center * binright_1 + binright_2_center * binright_2) / (binright_1 + binright_2);
				contour_x_right_weighted_twosig[yindex - 1] = weighted_right;
				//  cout << "Here's the weighted right" << weighted_right << " For smearing " << yindex << endl;
				// arrayright[yindex - 1] = weighted_right;
				break;
			}
		}
	}
}

std::vector<double> chisquaretest::createCustomBinning(
	TH1D *hist,
	double range1_min, double range1_max, int rebin1,
	double range2_min, double range2_max, int rebin2,
	double range3_min, double range3_max, int rebin3,
	double range4_min = -1, double range4_max = -1, int rebin4 = 1) // Defaults for 3-region rebinning
{
	// Get the original histogram range and binning
	double xMin = hist->GetXaxis()->GetXmin();
	double xMax = hist->GetXaxis()->GetXmax();
	int nBins = hist->GetNbinsX();
	double binWidth = (xMax - xMin) / nBins; // Original bin width

	std::vector<double> newBins;

	// Function to check if the rebinning factor is proper
	auto is_proper_rebin = [](double range_min, double range_max, double binWidth, int rebin) -> bool
	{
		if (range_min < 0 || range_max < 0)
			return true; // Ignore invalid ranges
		double nOriginalBins = (range_max - range_min) / binWidth;
		return std::fmod(nOriginalBins, rebin) == 0; // Proper if nOriginalBins is divisible by rebin
	};

	// Check rebinning validity
	if (!is_proper_rebin(range1_min, range1_max, binWidth, rebin1))
		std::cerr << "Warning: Rebinning factor for range [" << range1_min << ", " << range1_max << "] is improper.\n";
	if (!is_proper_rebin(range2_min, range2_max, binWidth, rebin2))
		std::cerr << "Warning: Rebinning factor for range [" << range2_min << ", " << range2_max << "] is improper.\n";
	if (!is_proper_rebin(range3_min, range3_max, binWidth, rebin3))
		std::cerr << "Warning: Rebinning factor for range [" << range3_min << ", " << range3_max << "] is improper.\n";
	if (!is_proper_rebin(range4_min, range4_max, binWidth, rebin4) && range4_min >= 0)
		std::cerr << "Warning: Rebinning factor for range [" << range4_min << ", " << range4_max << "] is improper.\n";

	// Step 1: Loop over the histogram bins and create bin edges
	double currentEdge = xMin;

	while (currentEdge < xMax)
	{
		// Apply rebinning only to valid ranges
		if (currentEdge >= range1_min && currentEdge < range1_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin1 * binWidth;
		}
		else if (currentEdge >= range2_min && currentEdge < range2_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin2 * binWidth;
		}
		else if (currentEdge >= range3_min && currentEdge < range3_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin3 * binWidth;
		}
		else if (range4_min >= 0 && currentEdge >= range4_min && currentEdge < range4_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin4 * binWidth;
		}
		else
		{
			newBins.push_back(currentEdge);
			currentEdge += binWidth; // Regular bin width for non-rebinning regions
		}
	}

	newBins.push_back(xMax); // Ensure the last bin edge is included

	return newBins;
}

void chisquaretest::drawcontour(TGraph *onesig_left, TGraph *onesig_right, TGraph *twosig_left, TGraph *twosig_right, int iteration, bool iseta, bool ispp, TString contourtitle, TString contoursaving)
{

	onesig_left = RemoveInvalidPoints(onesig_left);
	onesig_right = RemoveInvalidPoints(onesig_right);
	twosig_left = RemoveInvalidPoints(twosig_left);
	twosig_right = RemoveInvalidPoints(twosig_right);

	onesig_left = CombineGraphsToCircle(onesig_left, onesig_right);
	twosig_left = CombineGraphsToCircle(twosig_left, twosig_right);

	onesig_left->GetXaxis()->SetLabelSize(0.03);
	onesig_left->GetYaxis()->SetLabelSize(0.03);
	TCanvas *temp_c1 = new TCanvas("temp_c1", "", 1000, 1000);
	temp_c1->cd();
	temp_c1->SetRightMargin(0.05);
	temp_c1->SetLeftMargin(0.14);
	temp_c1->SetTickx(1);
	temp_c1->SetTicky(1);

	if (iseta && !ispp)
	{
		onesig_left->GetXaxis()->SetLimits(eta_mass_shift_array_low[iteration], eta_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetLimits(eta_mass_smear_array_low[iteration], eta_mass_smear_array_high[iteration]);
		onesig_left->GetXaxis()->SetRangeUser(eta_mass_shift_array_low[iteration], eta_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetRangeUser(eta_mass_smear_array_low[iteration], eta_mass_smear_array_high[iteration]);
	}
	if (!iseta && !ispp)
	{
		onesig_left->GetXaxis()->SetLimits(raw_mass_shift_array_low[iteration], raw_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetLimits(raw_mass_smear_array_low[iteration], raw_mass_smear_array_high[iteration]);
		onesig_left->GetXaxis()->SetRangeUser(raw_mass_shift_array_low[iteration], raw_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetRangeUser(raw_mass_smear_array_low[iteration], raw_mass_smear_array_high[iteration]);
	}
	/*if (ispp)
	{
		if (iseta)
		{
			if (isbk)
			{
				onesig_left->GetXaxis()->SetLimits(eta_pp_mass_shift_low_with_bk, eta_pp_mass_shift_high_with_bk);
				onesig_left->GetYaxis()->SetLimits(eta_pp_smear_low_with_bk, eta_pp_smear_high_with_bk);
				onesig_left->GetXaxis()->SetRangeUser(eta_pp_mass_shift_low_with_bk, eta_pp_mass_shift_high_with_bk);
				onesig_left->GetYaxis()->SetRangeUser(eta_pp_smear_low_with_bk, eta_pp_smear_high_with_bk);
			}
			if (!isbk)
			{
				onesig_left->GetXaxis()->SetLimits(eta_pp_mass_shift_low_without_bk, eta_pp_mass_shift_high_without_bk);
				onesig_left->GetYaxis()->SetLimits(eta_pp_smear_low_without_bk, eta_pp_smear_high_without_bk);
				onesig_left->GetXaxis()->SetRangeUser(eta_pp_mass_shift_low_without_bk, eta_pp_mass_shift_high_without_bk);
				onesig_left->GetYaxis()->SetRangeUser(eta_pp_smear_low_without_bk, eta_pp_smear_high_without_bk);
			}
		}
		if (!iseta)
		{
			if (isbk)
			{
				onesig_left->GetXaxis()->SetLimits(raw_pp_mass_shift_low_with_bk, raw_pp_mass_shift_high_with_bk);
				onesig_left->GetYaxis()->SetLimits(raw_pp_smear_low_with_bk, raw_pp_smear_high_with_bk);
				onesig_left->GetXaxis()->SetRangeUser(raw_pp_mass_shift_low_with_bk, raw_pp_mass_shift_high_with_bk);
				onesig_left->GetYaxis()->SetRangeUser(raw_pp_smear_low_with_bk, raw_pp_smear_high_with_bk);
			}
			if (!isbk)
			{
				onesig_left->GetXaxis()->SetLimits(raw_pp_mass_shift_low_without_bk, raw_pp_mass_shift_high_without_bk);
				onesig_left->GetYaxis()->SetLimits(raw_pp_smear_low_without_bk, raw_pp_smear_high_without_bk);
				onesig_left->GetXaxis()->SetRangeUser(raw_pp_mass_shift_low_without_bk, raw_pp_mass_shift_high_without_bk);
				onesig_left->GetYaxis()->SetRangeUser(raw_pp_smear_low_without_bk, raw_pp_smear_high_without_bk);
			}
		}
	}*/
	// not solved until I finished modify pp

	TH1 *frame = onesig_left->GetHistogram(); // Get the underlying histogram for customization

	onesig_left->SetTitle(contourtitle);

	onesig_left->GetXaxis()->SetTitle("Shifted Amount (GeV)");
	onesig_left->GetYaxis()->SetTitle("Smeared Amount (GeV)");

	onesig_left->SetMarkerSize(2);
	onesig_left->SetMarkerColor(kGreen - 3);
	onesig_left->SetLineWidth(3);		   // Set line width to 2
	onesig_left->SetLineColor(kGreen - 3); // Set line color to blue
	onesig_left->SetLineStyle(2);
	onesig_left->SetMarkerStyle(21);

	onesig_right->SetMarkerSize(2);
	onesig_right->SetMarkerColor(kGreen);
	onesig_right->SetLineWidth(3);		// Set line width to 2
	onesig_right->SetLineColor(kGreen); // Set line color to blue
	onesig_right->SetMarkerStyle(21);

	twosig_left->SetMarkerSize(2);
	twosig_left->SetMarkerColor(kRose);
	twosig_left->SetLineWidth(3); // Set line width to 2
	twosig_left->SetLineStyle(2);
	twosig_left->SetLineColor(kRose); // Set line color to blue
	twosig_left->SetMarkerStyle(21);

	twosig_right->SetMarkerSize(2);
	twosig_right->SetMarkerColor(kMagenta);
	twosig_right->SetLineWidth(3);		  // Set line width to 2
	twosig_right->SetLineColor(kMagenta); // Set line color to blue
	twosig_right->SetMarkerStyle(21);

	onesig_left->Draw("AP L");

	for (int i = 0; i < onesig_left->GetN(); ++i)
	{
		double x, y;
		onesig_left->GetPoint(i, x, y);
		// std::cout << "Point " << i << ": (" << x << ", " << y << ")" << std::endl;
	}
	// onesig_right->Draw("PL SAME");
	twosig_left->Draw("PL SAME");
	// twosig_right->Draw("PL SAME");

	TLegend *legend = new TLegend(0.8, 0.75, 0.9, 0.85); // x1, y1, x2, y2 in NDC (normalized device coordinates)

	// Add entries to the legend
	legend->AddEntry(onesig_left, "1#sigma CL", "PL");
	// legend->AddEntry(onesig_right, "1#sigma Right", "P");
	// legend->AddEntry(twosig_left, "2#sigma Left", "P");
	legend->AddEntry(twosig_left, "2#sigma CL", "PL");

	legend->SetTextSize(0.03);				// Set text size
	legend->SetTextFont(42);				// Use a modern, clean font
	legend->SetBorderSize(0);				// Set border size (0 for no border)
	legend->SetLineColor(kBlack);			// Border color (if any)
	legend->SetLineWidth(0);				// Border line width
	legend->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
	// legend->SetShadowColor(kGray);			// Add a subtle shadow effect

	// Draw the legend
	legend->Draw();
	temp_c1->SaveAs(contoursaving);
}

TGraph *chisquaretest::RemoveInvalidPoints(TGraph *originalGraph)
{
	std::vector<double> validX;
	std::vector<double> validY;
	std::vector<double> validEX;
	std::vector<double> validEY;

	int nPoints = originalGraph->GetN();
	for (int i = 0; i < nPoints; ++i)
	{
		double x, y;
		// double ex = originalGraph->GetErrorX(i);
		// double ey = originalGraph->GetErrorY(i);
		originalGraph->GetPoint(i, x, y);

		// Exclude points where x == -99
		if (x != -99)
		{
			validX.push_back(x);
			validY.push_back(y);
			// validEX.push_back(ex);
			// validEY.push_back(ey);
		}
	}

	// Create a new graph with valid points
	auto *cleanedGraph = new TGraph(validX.size(), validX.data(), validY.data());
	return cleanedGraph;
}

TGraph *chisquaretest::CombineGraphsToCircle(TGraph *graph1, TGraph *graph2)
{
	// Vectors to hold combined points
	std::vector<double> xCombined, yCombined;

	// Get points from the first graph (left semi-circle)
	int n1 = graph1->GetN();
	for (int i = 0; i < n1; ++i)
	{
		double x, y;
		graph1->GetPoint(i, x, y);
		xCombined.push_back(x);
		yCombined.push_back(y);
	}

	// Get points from the second graph (right semi-circle)
	int n2 = graph2->GetN();
	std::vector<double> xRight, yRight;
	for (int i = 0; i < n2; ++i)
	{
		double x, y;
		graph2->GetPoint(i, x, y);
		xRight.push_back(x);
		yRight.push_back(y);
	}

	// Reverse the order of the second graph points
	std::reverse(xRight.begin(), xRight.end());
	std::reverse(yRight.begin(), yRight.end());

	// Append the reversed points to the combined vectors
	xCombined.insert(xCombined.end(), xRight.begin(), xRight.end());
	yCombined.insert(yCombined.end(), yRight.begin(), yRight.end());

	// Add the first point from the left semi-circle to close the loop
	double xStart, yStart;
	graph1->GetPoint(0, xStart, yStart);
	xCombined.push_back(xStart);
	yCombined.push_back(yStart);

	// Create the new graph with the combined points
	TGraph *combinedGraph = new TGraph(xCombined.size(), xCombined.data(), yCombined.data());
	return combinedGraph;
}
Int_t chisquaretest::mapthreecases(int type, bool isbk, bool iseff)
{
	// Trying to return a Int such that I don't need to deal with all those if statement.
	if (type == 0)
	{
		if (isbk)
		{
			if (iseff)
			{
				// raw, bk, eff
				return 1;
			}
			if (!iseff)
			{ // raw, bk
				return 2;
			}
		}
		if (!isbk)
		{
			if (iseff)
			{
				// raw, eff
				return 3;
			}
			if (!iseff)
			{
				// raw
				return 4;
			}
		}
	}
	if (type == 1)
	{
		if (isbk)
		{
			if (iseff)
			{
				// eta, bk, eff
				return 5;
			}
			if (!iseff)
			{
				// eta, bk
				return 6;
			}
		}
		if (!isbk)
		{
			if (iseff)
			{
				// eta, eff
				return 7;
			}
			if (!iseff)
			{
				// eta
				return 8;
			}
		}
	}
	return -99;
}
void chisquaretest::saveChi2Region(TH2D *hist, int binX_min, int binY_min, int my_case, bool iseta, int iteration, int region_size = 3, bool ispp = true)
{
	if (region_size % 2 == 0)
	{
		std::cerr << "Warning: Region size (" << region_size
				  << ") is even. It is recommended to use an odd size for symmetrical regions.\n";
	}

	// Define the range (symmetric or asymmetric for even region_size)
	int half_size = region_size / 2;

	// Calculate the range in bins
	int binX_low = std::max(1, binX_min - half_size);
	int binX_high = std::min(hist->GetNbinsX(), binX_min + half_size - (region_size % 2 == 0 ? 1 : 0));
	int binY_low = std::max(1, binY_min - half_size);
	int binY_high = std::min(hist->GetNbinsY(), binY_min + half_size - (region_size % 2 == 0 ? 1 : 0));

	// Get the corresponding axis values
	double x_low = hist->GetXaxis()->GetBinLowEdge(binX_low);
	double x_high = hist->GetXaxis()->GetBinUpEdge(binX_high);
	double y_low = hist->GetYaxis()->GetBinLowEdge(binY_low);
	double y_high = hist->GetYaxis()->GetBinUpEdge(binY_high);

	double x_min = hist->GetXaxis()->GetBinCenter(binX_min);
	double y_min = hist->GetYaxis()->GetBinCenter(binY_min);

	TMarker *marker_L = new TMarker(x_low, y_min, 20); // Marker type 29 (big star)
	marker_L->SetMarkerColor(kGreen);
	marker_L->SetMarkerSize(3.5);
	marker_L->Draw("same");

	TMarker *marker_R = new TMarker(x_high, y_min, 20); // Marker type 29 (big star)
	marker_R->SetMarkerColor(kGreen);
	marker_R->SetMarkerSize(3.5);
	marker_R->Draw("same");

	TMarker *marker_T = new TMarker(x_min, y_high, 20); // Marker type 29 (big star)
	marker_T->SetMarkerColor(kGreen);
	marker_T->SetMarkerSize(3.5);
	marker_T->Draw("same");

	TMarker *marker_B = new TMarker(x_min, y_low, 20); // Marker type 29 (big star)
	marker_B->SetMarkerColor(kGreen);
	marker_B->SetMarkerSize(3.5);
	marker_B->Draw("same");

	// Save to a file
	TString txtname = "";
	TString txtname_prefix = "";
	TString cent = "";
	TString anotherprefix = "";

	if (my_case == 1)
		txtname = "nominal";
	if (my_case == 2)
		txtname = "tnpU";
	if (my_case == 3)
		txtname = "tnpD";
	if (my_case == 4)
		txtname = "acooff";
	if (my_case == 5)
		txtname = "nominal_no_bk_sub";
	if (my_case == 6)
		txtname = "nominal_binning";
	if (my_case == 7)
		txtname = "nominal_range";

	if (!ispp)
	{
		if (iteration == 0)
			cent = "_0-10";
		if (iteration == 1)
			cent = "_10-20";
		if (iteration == 2)
			cent = "_20-30";
		if (iteration == 3)
			cent = "_30-100";
		if (iteration == 10)
			cent = "_0-100";
	}
	if (ispp)
	{
		cent = Form("_%i", iteration);
	}

	if (ispp)
		txtname_prefix = "pp_";
	if (!ispp)
		txtname_prefix = "PbPb_";

	if (iseta)
		anotherprefix = "eta_";
	if (!iseta)
		anotherprefix = "FA_";

	std::ofstream outfile("./zoomin/" + txtname_prefix + anotherprefix + txtname + cent + ".txt");

	outfile << std::fixed << std::setprecision(4); // Set fixed-point notation with 2 decimal places

	outfile << "Local Minimum Bin: (" << binX_min << ", " << binY_min << ")\n";
	outfile << "Region Size: " << region_size << "x" << region_size << "\n";
	outfile << "X-axis Range: [" << x_low << ", " << x_high << "]\n";
	outfile << "Y-axis Range: [" << y_low << ", " << y_high << "]\n";
	outfile.close();
}
TH1D *chisquaretest::reducebin(TH1D *h_old_data, bool iseta, bool isleft)
{
	if (!h_old_data)
	{
		std::cerr << "Error: Null pointer passed to reducebin()" << std::endl;
		return nullptr;
	}

	int old_bins = h_old_data->GetNbinsX();
	int bins_to_remove = (iseta) ? 2 : 5; // Remove 2 bins if iseta=true, 5 otherwise
	int new_bins = old_bins - bins_to_remove;

	if (new_bins <= 0)
	{
		std::cerr << "Error: Too many bins removed. No bins left!" << std::endl;
		return nullptr;
	}

	// Retrieve bin edges
	std::vector<double> bin_edges(old_bins + 1);
	for (int i = 0; i <= old_bins; i++)
	{
		bin_edges[i] = h_old_data->GetXaxis()->GetBinLowEdge(i + 1);
	}

	// Define the new bin edges based on the removal direction
	std::vector<double> new_bin_edges;
	if (isleft)
	{
		new_bin_edges.assign(bin_edges.begin() + bins_to_remove, bin_edges.end());
	}
	else
	{
		new_bin_edges.assign(bin_edges.begin(), bin_edges.end() - bins_to_remove);
	}

	// Create the new histogram with variable binning
	TH1D *h_new = new TH1D(Form("%s_reduced", h_old_data->GetName()), h_old_data->GetTitle(), new_bins, new_bin_edges.data());

	// Copy bin contents while shifting appropriately
	if (isleft)
	{
		for (int i = bins_to_remove + 1; i <= old_bins; i++)
		{
			int new_index = i - bins_to_remove;
			h_new->SetBinContent(new_index, h_old_data->GetBinContent(i));
			h_new->SetBinError(new_index, h_old_data->GetBinError(i));
		}
	}
	else
	{
		for (int i = 1; i <= old_bins - bins_to_remove; i++)
		{
			h_new->SetBinContent(i, h_old_data->GetBinContent(i));
			h_new->SetBinError(i, h_old_data->GetBinError(i));
		}
	}

	return h_new;
}

TH1D *chisquaretest::reducebinpp(TH1D *h_old_data, bool isleft)
{
	if (!h_old_data)
	{
		std::cerr << "Error: Null pointer passed to reducebin()" << std::endl;
		return nullptr;
	}

	int old_bins = h_old_data->GetNbinsX();
	int bins_to_remove = 20;
	int new_bins = old_bins - bins_to_remove;

	if (new_bins <= 0)
	{
		std::cerr << "Error: Too many bins removed. No bins left!" << std::endl;
		return nullptr;
	}

	// Retrieve bin edges
	std::vector<double> bin_edges(old_bins + 1);
	for (int i = 0; i <= old_bins; i++)
	{
		bin_edges[i] = h_old_data->GetXaxis()->GetBinLowEdge(i + 1);
	}

	// Define the new bin edges based on the removal direction
	std::vector<double> new_bin_edges;
	if (isleft)
	{
		new_bin_edges.assign(bin_edges.begin() + bins_to_remove, bin_edges.end());
	}
	else
	{
		new_bin_edges.assign(bin_edges.begin(), bin_edges.end() - bins_to_remove);
	}

	// Create the new histogram with variable binning
	TH1D *h_new = new TH1D(Form("%s_pp_reduced", h_old_data->GetName()), h_old_data->GetTitle(), new_bins, new_bin_edges.data());

	// Copy bin contents while shifting appropriately
	if (isleft)
	{
		for (int i = bins_to_remove + 1; i <= old_bins; i++)
		{
			int new_index = i - bins_to_remove;
			h_new->SetBinContent(new_index, h_old_data->GetBinContent(i));
			h_new->SetBinError(new_index, h_old_data->GetBinError(i));
		}
	}
	else
	{
		for (int i = 1; i <= old_bins - bins_to_remove; i++)
		{
			h_new->SetBinContent(i, h_old_data->GetBinContent(i));
			h_new->SetBinError(i, h_old_data->GetBinError(i));
		}
	}

	return h_new;
}
void chisquaretest::readlimit(int type, bool iseta, int cent)
{
	TString variation[7] = {"nominal", "tnpU", "tnpD", "acooff", "nominal_no_bk_sub", "nominal_binning", "nominal_range"};
	TString filename = "";
	TString type_str = "";
	TString savedname = "";

	if (iseta == 1)
		type_str = "eta_";
	if (iseta == 0)
		type_str = "FA_";

	filename = TString("./zoomin/") + TString("PbPb_") + type_str + variation[type - 1] + "_" + Form("%i-%i.txt", this->cenlowlimit[cent], this->cenhighlimit[cent]);

	std::ifstream file(filename);
	if (!file)
	{
		std::cerr << "Error: Unable to open file " << filename << std::endl;
		return;
	}

	std::vector<double> ranges(4); // Stores xmin, xmax, ymin, ymax
	std::string line;

	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		std::string key;
		if (line.find("X-axis Range:") != std::string::npos)
		{
			char ch;
			iss >> key >> key >> ch >> ranges[0] >> ch >> ranges[1]; // Extracts xmin, xmax
		}
		else if (line.find("Y-axis Range:") != std::string::npos)
		{
			char ch;
			iss >> key >> key >> ch >> ranges[2] >> ch >> ranges[3]; // Extracts ymin, ymax
		}
	}

	file.close();

	// Debug print
	std::cout << "Extracted ranges: [" << ranges[0] << ", " << ranges[1] << ", " << ranges[2] << ", " << ranges[3] << "]" << std::endl;
	placeholder_mass_shift_array_low_zoomin = ranges[0];
	placeholder_mass_shift_array_high_zoomin = ranges[1];
	placeholder_mass_smear_array_low_zoomin = ranges[2];
	placeholder_mass_smear_array_high_zoomin = ranges[3];
}

void chisquaretest::readlimitpp(int type, bool iseta)
{
	TString variation[7] = {"nominal", "tnpU", "tnpD", "acooff", "nominal_no_bk_sub", "nominal_binning", "nominal_range"};
	TString filename = "";
	TString type_str = "";
	TString savedname = "";

	if (iseta == 1)
		type_str = "eta_";
	if (iseta == 0)
		type_str = "FA_";

	filename = TString("./zoomin/") + TString("pp_") + type_str + variation[type - 1] + "_" + Form("%i.txt", 22);

	std::ifstream file(filename);
	if (!file)
	{
		std::cerr << "Error: Unable to open file " << filename << std::endl;
		return;
	}

	std::vector<double> ranges(4); // Stores xmin, xmax, ymin, ymax
	std::string line;

	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		std::string key;
		if (line.find("X-axis Range:") != std::string::npos)
		{
			char ch;
			iss >> key >> key >> ch >> ranges[0] >> ch >> ranges[1]; // Extracts xmin, xmax
		}
		else if (line.find("Y-axis Range:") != std::string::npos)
		{
			char ch;
			iss >> key >> key >> ch >> ranges[2] >> ch >> ranges[3]; // Extracts ymin, ymax
		}
	}

	file.close();

	// Debug print
	std::cout << "Extracted ranges: [" << ranges[0] << ", " << ranges[1] << ", " << ranges[2] << ", " << ranges[3] << "]" << std::endl;
	placeholder_pp_mass_shift_array_low_zoomin = ranges[0];
	placeholder_pp_mass_shift_array_high_zoomin = ranges[1];
	placeholder_pp_mass_smear_array_low_zoomin = ranges[2];
	placeholder_pp_mass_smear_array_high_zoomin = ranges[3];
}