#include "plotting_helper.h"
#include "tdrStyle.C"
#include "CMS_lumi.C"

class chisquaretest : public plotting_helper
{
public:
	chisquaretest();
	chisquaretest(TString s2, TString s3, int type, bool isNew);
	chisquaretest(TString s1, TString s2, TString s3, bool iseta, bool isNew);
	~chisquaretest();
	Double_t myownfunctionchi2(TH1D *h1, TH1D *h2);
	void calculatechisq(bool isbk);
	void calculatechisqpp(bool isbk);
	void plottingandformatting(int type, bool isbk);
	void plottingandformattingpp(bool iseta, bool isbk);
	void bincontentcheck(bool isbk);
	void RebinAll(int type);
	TGraph *RemoveInvalidPoints(TGraph *originalGraph);
	TGraph *CombineGraphsToCircle(TGraph *graph1, TGraph *graph2);
	void RebinAllpp(int x);
	Double_t getuncertainty(TH2D *h_1, int type, Int_t minBinX, Int_t minBinY, int ndf);
	void getcontour(TH2D *h1, int type, int ndf, Int_t minBinX, Int_t minBinY, Double_t minBinContent, Double_t *arrayleft, Double_t *arrayright);
	void drawcontour(TGraph *onesig_left, TGraph *onesig_right, TGraph *twosig_left, TGraph *twosig_right, int iteration, bool isbk, bool iseta, bool ispp);
	std::vector<double> createCustomBinning(
		TH1D *hist,
		double range1_min, double range1_max, int rebin1,
		double range2_min, double range2_max, int rebin2,
		double range3_min, double range3_max, int rebin3,
		double range4_min, double range4_max, int rebin4);

	static const int nbins_mass_shift = 21;
	static const int nbins_smear = 21;
	static const int nbins_cent = 11;

	// those default values now only affect pp

	double lowbin_mass_shift = -0.15;
	double highbin_mass_shift = -0.1;
	double lowbin_smear = 0.006;
	double highbin_smear = 0.012;

	double h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
	double h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

	double h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
	double h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;

	TH1D *h_mc_signal[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *h_data[nbins_cent];
	TH1D *h_mc_bk[nbins_cent];
	TH1D *h_data_bksub[nbins_cent];
	TH2D *h_chisquare[nbins_cent];

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

	TH1D *h_mc_signal_pp[nbins_mass_shift][nbins_smear][22];
	TH1D *h_data_pp[22];
	TH1D *h_mc_bk_pp;
	TH1D *h_data_bksub_pp[22];
	TH2D *h_chisquare_pp[22];

	Double_t dMass_pp[22];
	Double_t dMass_Err_pp[22];

	Double_t dWidth_pp[22];
	Double_t dWidth_Err_pp[22];

	Double_t xposition[22];
	Double_t xposition_err[22];

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
	TCanvas *c_contour_HI[nbins_cent];
	// What I need: invariant mass comparison data vs mc raw, data vs data bk sub.
	TCanvas *c_data_mc_raw[nbins_cent];
	TCanvas *c_data_data_bk[nbins_cent];

	TCanvas *c_2d_chisquare_ndf_pp[22];
	TCanvas *c_data_mc_raw_pp[22];
	TCanvas *c_data_data_bk_pp[22];

	TFile *mcfile;
	TFile *datafile;
	TFile *bkfile;

	double eta_mass_shift_array_low[nbins_cent] = {-0.45, -0.3, -0.45, -0.5, 0, 0, 0, 0, 0, 0, -0.35};
	double eta_mass_shift_array_high[nbins_cent] = {-0.05, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0.0};
	double eta_mass_smear_array_low[nbins_cent] = {-0.25, -0.05, -0.15, -0.3, 0, 0, 0, 0, 0, 0, 0};
	double eta_mass_smear_array_high[nbins_cent] = {0.3, 0.65, 0.65, 0.5, 0, 0, 0, 0, 0, 0, 0.35};

	double raw_mass_shift_array_low[nbins_cent] = {-0.3, -0.35, -0.4, -0.3, 0, 0, 0, 0, 0, 0, -0.25};
	double raw_mass_shift_array_high[nbins_cent] = {0.05, 0.0, 0.0, 0.1, 0, 0, 0, 0, 0, 0, 0.0};
	double raw_mass_smear_array_low[nbins_cent] = {0.05, 0.15, -0.1, 0.2, 0, 0, 0, 0, 0, 0, 0.2};
	double raw_mass_smear_array_high[nbins_cent] = {0.4, 0.6, 0.4, 0.7, 0, 0, 0, 0, 0, 0, 0.45};

	double eta_pp_mass_shift_low = -0.2;
	double eta_pp_mass_shift_high = 0.05;
	double eta_pp_smear_low = 0.2;
	double eta_pp_smear_high = 0.4;

	double raw_pp_mass_shift_low = -0.15;
	double raw_pp_mass_shift_high = 0.05;
	double raw_pp_smear_low = 0.26;
	double raw_pp_smear_high = 0.33;

};
chisquaretest::chisquaretest()
{
}

chisquaretest::chisquaretest(TString s2, TString s3, int type, bool isNew)
{

	cout << "We are running HI, with bin dim " << nbins_mass_shift << " * " << nbins_smear << " " << " Shift: " << lowbin_mass_shift << " " << highbin_mass_shift << " Smear: " << lowbin_smear << " " << highbin_smear << endl;
	// std::this_thread::sleep_for(std::chrono::seconds(3));

	datafilepath = s2;
	bkfilepath = s3;

	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath, "READ");
	mcfile = new TFile("../ZBoson_18/rootfile/new_template_reco_gen.root", "READ");

	for (int i = 0; i < 5; i++)
	{
		xposition_HI[i] = i + 1;
		xposition_err_HI[i] = 0;
	}

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		if (!((cent < 4) || (cent == 10)))
			continue;

		cout << "cent is " << cent << endl;

		/*if (type == 1)
		{
			if (cent == 0)
			{
				// mcfilepath = "../ZBoson_18/rootfile/eta_0_10_shift_-0.44_-0.33_smear_0_0.0085_modified_signal_21_21_100.root";
				mcfilepath = "../ZBoson_18/rootfile/eta_0_10_shift_-0.45_-0.05_smear_0_0.01_modified_signal_21_21_1000.root";
			}
			if (cent == 1)
			{
				mcfilepath = "../ZBoson_18/rootfile/eta_10_20_shift_-0.3_0.2_smear_0_0.01_modified_signal_21_21_1000.root";
			}
			if (cent == 2)
			{
				mcfilepath = "../ZBoson_18/rootfile/eta_20_30_shift_-0.45_0.1_smear_0_0.01_modified_signal_21_21_1000.root";
			}
			if (cent == 3)
			{
				mcfilepath = "../ZBoson_18/rootfile/eta_30_100_shift_-0.5_0.0_smear_0_0.01_modified_signal_21_21_1000.root";
			}
			if (cent == 10)
			{
				mcfilepath = "../ZBoson_18/rootfile/eta_0_100_shift_-0.34_-0.08_smear_0_0.007_modified_signal_21_21_1000.root";
			}
		}
		if (type == 0)
		{
			if (cent == 0)
			{
				mcfilepath = "../ZBoson_18/rootfile/raw_0_10_shift_-0.3_0.0_smear_0_0.008_modified_signal_21_21_1000.root";
			}
			if (cent == 1)
			{
				mcfilepath = "../ZBoson_18/rootfile/raw_10_20_shift_-0.4_-0.05_smear_0.001_0.012_modified_signal_21_21_1000.root";
			}
			if (cent == 2)
			{
				mcfilepath = "../ZBoson_18/rootfile/raw_20_30_shift_-0.4_0.0_smear_0.0_0.008_modified_signal_21_21_1000.root";
			}
			if (cent == 3)
			{
				mcfilepath = "../ZBoson_18/rootfile/raw_30_100_shift_-0.4_0.05_smear_0.002_0.013_modified_signal_21_21_1000.root";
			}
			if (cent == 10)
			{
				// This is the same as 0-10
				mcfilepath = "../ZBoson_18/rootfile/raw_0_10_shift_-0.3_0.0_smear_0_0.008_modified_signal_21_21_1000.root";
			}
		}*/

		// mcfile = new TFile(mcfilepath, "READ");
		if (type == 1)
		{
			this->lowbin_mass_shift = eta_mass_shift_array_low[cent];
			this->highbin_mass_shift = eta_mass_shift_array_high[cent];
			this->lowbin_smear = eta_mass_smear_array_low[cent];
			this->highbin_smear = eta_mass_smear_array_high[cent];

			this->h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
			this->h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

			this->h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
			this->h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		}
		if (type == 0)
		{
			this->lowbin_mass_shift = raw_mass_shift_array_low[cent];
			this->highbin_mass_shift = raw_mass_shift_array_high[cent];
			this->lowbin_smear = raw_mass_smear_array_low[cent];
			this->highbin_smear = raw_mass_smear_array_high[cent];

			this->h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
			this->h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

			this->h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
			this->h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		}

		for (int i = 0; i < nbins_smear; i++)
		{
			double bincenter = ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
			// FIX ME !!!
			contour_y_HI[cent][i] = 2 * (i)*bincenter + lowbin_smear;
		}

		c_2d_chisquare_ndf[cent] = new TCanvas(Form("c_2d_chisquare_ndf_%i", cent), "", 3200, 2400);
		c_data_mc_raw[cent] = new TCanvas(Form("c_data_mc_raw_%i", cent), "", 800, 800);
		c_data_data_bk[cent] = new TCanvas(Form("c_data_data_bk_%i", cent), "", 800, 600);
		c_contour_HI[cent] = new TCanvas(Form("c_contour_HI_%i", cent), "", 800, 600);

		if (type == 0)
			h_data[cent] = (TH1D *)datafile->Get(Form("mass_array_data_%i", cent));
		if (type == 1)
			h_data[cent] = (TH1D *)datafile->Get(Form("mass_array_data_witheta_%i", cent));

		h_mc_bk[cent] = (TH1D *)bkfile->Get(Form("Normalized_mc_bk_%i", cent));
		h_chisquare[cent] = new TH2D(Form("h_chisquare_%i", cent), Form("Cent_%i_%i", this->cenlowlimit[cent], this->cenhighlimit[cent]), nbins_mass_shift, h_low_mass_shift, h_high_mass_shift, nbins_smear, h_low_smear, h_high_smear); // Equation here: half bin to the left and half bin to the right, bin width = range / (21-1)

		this->areanormalize(h_data[cent]);

		h_data_bksub[cent] = (TH1D *)h_data[cent]->Clone();
		h_data_bksub[cent]->Add(h_mc_bk[cent], -1);

		// cout << "nbins is " << h_data_bksub[cent]->GetNbinsX() << endl;

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				if (type == 0)
				{
					if (isNew)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("mass_array_with_eff_template_%i_%i_%i", shift, smear, cent));
					}
					else
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_raw_without_eff_%i_%i_%i", shift, smear, cent));
					}
				}
				if (type == 1)
				{
					if (isNew)
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("mass_array_witheta_witheff_template_%i_%i_%i", shift, smear, cent));
					}
					else
					{
						h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_eta_without_eff_%i_%i_%i", shift, smear, cent));
					}
				}

				this->areanormalize(h_mc_signal[shift][smear][cent]);
			}
		}
	}
}

chisquaretest::chisquaretest(TString s1, TString s2, TString s3, bool iseta, bool isNew)
{
	cout << "We are running pp, with bin dim " << nbins_mass_shift << " * " << nbins_smear << " " << " Shift: " << lowbin_mass_shift << " " << highbin_mass_shift << " Smear: " << lowbin_smear << " " << highbin_smear << endl;
	// std::this_thread::sleep_for(std::chrono::seconds(3));

	// mcfilepath = s1;
	datafilepath = s2;
	bkfilepath = s3;

	if (iseta)
		mcfilepath = "../ZBoson_18/rootfile/new_template_pp_reco_gen.root";
		//mcfilepath = "../ZBoson_18/rootfile/new_template_pp_reco_gen_raw_test_-1_1_0_2.root";
	if (!iseta)
	{
		mcfilepath = "../ZBoson_18/rootfile/new_template_pp_reco_gen.root";
		//mcfilepath = "../ZBoson_18/rootfile/new_template_pp_reco_gen_raw_test_-1_1_0_2.root";
		//mcfilepath = "../ZBoson_18/rootfile/new_template_pp_reco_gen_raw_test_-1_1_-1_1.root";
	}

	// mcfilepath = "../ZBoson_18/rootfile/new_template_pp_reco_gen.root";
	// mcfilepath = "../ZBoson_18/rootfile/raw_pp_shift_-0.12_-0.08_smear_0.0085_0.012_modified_signal_21_21_1000.root";

	cout << "mcfilepath is " << mcfilepath << endl;

	mcfile = new TFile(mcfilepath, "READ");
	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath, "READ");

	if (iseta)
	{
		this->lowbin_mass_shift = eta_pp_mass_shift_low;
		this->highbin_mass_shift = eta_pp_mass_shift_high;
		this->lowbin_smear = eta_pp_smear_low;
		this->highbin_smear = eta_pp_smear_high;

		this->h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
		this->h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

		this->h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		this->h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
	}

	if (!iseta)
	{
		this->lowbin_mass_shift = raw_pp_mass_shift_low;
		this->highbin_mass_shift = raw_pp_mass_shift_high;
		this->lowbin_smear = raw_pp_smear_low;
		this->highbin_smear = raw_pp_smear_high;

		this->h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
		this->h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

		this->h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		this->h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
	}

	h_mc_bk_pp = (TH1D *)bkfile->Get("Normalized_mc_bk_10");

	for (int i = 0; i < nbins_smear; i++)
	{
		double bincenter = ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
		// FIX ME !!!
		contour_y_HI[10][i] = 2 * (i)*bincenter + lowbin_smear;
	}

	for (int runperiod = 0; runperiod < 22; runperiod++)
	{
		c_2d_chisquare_ndf_pp[runperiod] = new TCanvas(Form("c_2d_chisquare_ndf_pp_%i", runperiod), "", 3200, 2400);
		c_data_mc_raw_pp[runperiod] = new TCanvas(Form("c_data_mc_raw_pp_%i", runperiod), "", 800, 600);
		c_data_data_bk_pp[runperiod] = new TCanvas(Form("c_data_data_bk_pp_%i", runperiod), "", 800, 600);

		if (iseta)
		{
			// cout << "WE ARE HERE" << endl;
			h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("h_mass_array_eta_%i", runperiod));
		}
		else
		{
			h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("h_mass_array_raw_%i", runperiod));
		}

		h_chisquare_pp[runperiod] = new TH2D(Form("h_chisquare_pp_%i", runperiod), Form("h_chisquare_pp_%i", runperiod), nbins_mass_shift, h_low_mass_shift, h_high_mass_shift, nbins_smear, h_low_smear, h_high_smear);
		this->areanormalize(h_data_pp[runperiod]);

		h_data_bksub_pp[runperiod] = (TH1D *)h_data_pp[runperiod]->Clone();
		h_data_bksub_pp[runperiod]->Add(h_mc_bk_pp, -1);

		xposition[runperiod] = runperiod + 1;
		xposition_err[runperiod] = 0;

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				if (iseta)
				{
					if (isNew)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("mass_array_witheta_witheff_template_%i_%i_%i", shift, smear, 10));
					}
					else
					{
						//h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_eta_without_eff_%i_%i_%i", shift, smear, 10));
					}
					// h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_%i_%i_%i", shift, smear, 10));
				}
				else
				{
					if (isNew)
					{
						h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("mass_array_with_eff_template_%i_%i_%i", shift, smear, 10));
					}
					else
					{
						//h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_raw_without_eff_%i_%i_%i", shift, smear, 10));
					}

					// h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_%i_%i_%i", shift, smear, 10));
				}

				this->areanormalize(h_mc_signal_pp[shift][smear][runperiod]);
			}
		}
	}
}

void chisquaretest::calculatechisq(bool isbk)
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
				if (isbk)
					chisquarevalue = myownfunctionchi2(h_data_bksub[cent], h_mc_signal[shift][smear][cent]);
				if (!isbk)
					chisquarevalue = myownfunctionchi2(h_data[cent], h_mc_signal[shift][smear][cent]);
				std::ostringstream stream;
				stream << std::fixed << std::setprecision(2) << chisquarevalue;
				double formattedBinContent = std::stod(stream.str());
				h_chisquare[cent]->SetBinContent(shift + 1, smear + 1, formattedBinContent);
			}
		}
	}
}

void chisquaretest::calculatechisqpp(bool isbk)
{

	for (int runperiod = 0; runperiod < 22; runperiod++)
	{
		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				double chisquarevalue = 0;
				if (isbk)
					chisquarevalue = myownfunctionchi2(h_data_bksub_pp[runperiod], h_mc_signal_pp[shift][smear][runperiod]);
				if (!isbk)
					chisquarevalue = myownfunctionchi2(h_data_pp[runperiod], h_mc_signal_pp[shift][smear][runperiod]);

				if (runperiod == 0)
				{
					if (shift == 10 && (smear == 10 || smear == 0))
					{
						cout << "For smear " << smear << " For shift " << shift << " The chi2 value is " << chisquarevalue << endl;
						cout << "data name is " << h_data_pp[runperiod]->GetName()<< endl;
						cout << "Name is " << h_mc_signal_pp[shift][smear][runperiod]->GetName()<<endl; 

					}
				}
				std::ostringstream stream;
				stream << std::fixed << std::setprecision(2) << chisquarevalue;
				double formattedBinContent = std::stod(stream.str());
				h_chisquare_pp[runperiod]->SetBinContent(shift + 1, smear + 1, formattedBinContent);
			}
		}
	}
}

void chisquaretest::plottingandformatting(int type, bool isbk)
{

	TString chi2_title;
	TString chi2_saving_path;
	TString data_mc_title;
	TString data_data_title;
	TString data_mc_saving_path;
	TString data_data_saving_path;
	TString contour_saving_path;

	if (type == 0)
	{
		if (isbk)
		{
			chi2_title = "PbPb, |#eta| < 2.4 with bksub, centrality: (%i-%i)";
			data_mc_title = "WholeAcceptance, bksub, Cent:(%i-%i)";
			data_data_title = "Raw_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/raw/bksub/Raw_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/raw/bksub/Raw_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/bksub/raw_with_bksub_%i_%i.png";
		}
		if (!isbk)
		{
			chi2_title = "PbPb, |#eta| < 2.4 without bksub, centrality: (%i-%i)";
			data_mc_title = "WholeAcceptance, no bksub, Cent:(%i-%i)";
			data_data_title = "Raw_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/raw/nobksub/Raw_without_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/raw/nobksub/Raw_without_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_%i_%i.png";
			contour_saving_path = "./newchi2/contour/raw/nobksub/raw_without_bksub_%i_%i.png";
		}
	}
	if (type == 1)
	{
		if (isbk)
		{
			chi2_title = "PbPb, |#eta| < 1.0 with bksub, centrality: (%i-%i)";
			data_mc_title = "|#eta| < 1.0, bksub, centrality: (%i-%i)";
			data_data_title = "Eta_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/eta/bksub/Eta_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/eta/bksub/Eta_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/bksub/eta_with_bksub_%i_%i.png";
		}
		if (!isbk)
		{
			chi2_title = "PbPb, |#eta| < |1.0| without bksub, centrality: (%i-%i)";
			data_mc_title = "|#eta| < 1.0, no bksub, centrality: (%i-%i)";
			data_data_title = "Eta_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/eta/nobksub/Eta_without_bksub_%i_%i.png";
			chi2_saving_path = "./newchi2/chi2plots/eta/nobksub/Eta_without_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/eta/nobksub/Eta_without_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_%i_%i.png";
			contour_saving_path = "./newchi2/contour/eta/nobksub/raw_without_bksub_%i_%i.png";
		}
	}

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		if (!((cent < 4) || (cent == 10)))
			continue;
		// Chi2/ndf plot
		c_2d_chisquare_ndf[cent]->cd();
		// c_2d_chisquare_ndf[cent]->SetLogz();
		gStyle->SetPalette(kRainBow);
		h_chisquare[cent]->SetTitle(Form(chi2_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		h_chisquare[cent]->SetTitleFont(42);
		h_chisquare[cent]->Draw("COLZ");
		h_chisquare[cent]->Draw("TEXTSAME");
		h_chisquare[cent]->GetXaxis()->SetNdivisions(nbins_mass_shift, 0, 0);
		h_chisquare[cent]->GetYaxis()->SetNdivisions(nbins_smear, 0, 0);
		h_chisquare[cent]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare[cent]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare[cent]->GetXaxis()->SetTitle("Shifted Amount (GeV)");
		h_chisquare[cent]->GetYaxis()->SetTitle("Smeared Amount (GeV)");

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

		// Here's the place to get the contour
		cout << "Now we are running cent " << cenlowlimit[cent] << " " << cenhighlimit[cent] << endl;

		int numberofDF = -99;

		if (cent == 0)
			numberofDF = 44;
		if (cent == 1)
			numberofDF = 39;
		if (cent == 2)
			numberofDF = 29;
		if (cent == 3)
			numberofDF = 25;
		if (cent == 4)
			numberofDF = 11;
		if (cent == 10)
			numberofDF = 71;

		if (type == 0)
		{
			numberofDF = 27;
		}

		this->getcontour(h_chisquare[cent], 1, numberofDF, minBinX, minBinY, minContent, contour_x_left_onesig_HI, contour_x_right_onesig_HI);

		g_HI_contour_1sig_left = new TGraph(nbins_smear, contour_x_left_weighted_onesig, contour_y_HI[cent]);
		g_HI_contour_1sig_right = new TGraph(nbins_smear, contour_x_right_weighted_onesig, contour_y_HI[cent]);
		g_HI_contour_2sig_left = new TGraph(nbins_smear, contour_x_left_weighted_twosig, contour_y_HI[cent]);
		g_HI_contour_2sig_right = new TGraph(nbins_smear, contour_x_right_weighted_twosig, contour_y_HI[cent]);

		bool iseta = false;
		if (type == 1)
			iseta = true;

		this->drawcontour(g_HI_contour_1sig_left, g_HI_contour_1sig_right, g_HI_contour_2sig_left, g_HI_contour_2sig_right, cent, isbk, iseta, false);

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

		if (cent == 0 || cent == 1 || cent == 2 || cent == 3)
		{
			dMass_HI[cent + 1] = xCenter;
			dWidth_HI[cent + 1] = yCenter;

			//cout << "dWidth is " << dWidth_HI[cent + 1] << endl;

			dMass_Err_HI[cent + 1] = getuncertainty(h_chisquare[cent], 1, minBinX, minBinY, numberofDF);
			dWidth_Err_HI[cent + 1] = getuncertainty(h_chisquare[cent], 2, minBinX, minBinY, numberofDF);
		}

		if (cent == 10)
		{
			dMass_HI[0] = xCenter;
			dWidth_HI[0] = yCenter;

			dMass_Err_HI[0] = getuncertainty(h_chisquare[cent], 1, minBinX, minBinY, numberofDF);
			cout << "The Uncertainty for 0 - 100 is " << dMass_Err_HI[0] << endl;
			dWidth_Err_HI[0] = getuncertainty(h_chisquare[cent], 2, minBinX, minBinY, numberofDF);
		}

		cout << "yCenter is " << yCenter << endl;

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
		box1->SetLineColor(kRed);
		box1->SetLineWidth(4);
		box1->SetFillStyle(0);
		box1->Draw("same");

		c_2d_chisquare_ndf[cent]->SaveAs(Form(chi2_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));

		// This is Data and MC

		c_data_mc_raw[cent]->cd();
		c_data_mc_raw[cent]->SetLeftMargin(0.15);
		c_data_mc_raw[cent]->SetRightMargin(0.08);
		c_data_mc_raw[cent]->SetBottomMargin(0.13);
		c_data_mc_raw[cent]->SetTicks(1, 1);
		// c_data_mc_raw[cent]->SetLogy();

		if (isbk)
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
		if (!isbk)
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
	if (isbk && type == 1)
	{
		g_HI_dmass->Write("HI_dM_chi2_eta_bksub", 2);
		g_HI_dwidth->Write("HI_dWidth_chi2_eta_bksub", 2);
	}
	if (!isbk && type == 1)
	{
		g_HI_dmass->Write("HI_dM_chi2_eta", 2);
		g_HI_dwidth->Write("HI_dWidth_chi2_eta", 2);
	}
	if (isbk && type == 0)
	{
		g_HI_dmass->Write("HI_dM_chi2_raw_bksub", 2);
		g_HI_dwidth->Write("HI_dWidth_chi2_raw_bksub", 2);
	}
	if (!isbk && type == 0)
	{
		g_HI_dmass->Write("HI_dM_chi2_raw", 2);
		g_HI_dwidth->Write("HI_dWidth_chi2_raw", 2);
	}
	temp->Close();
}

void chisquaretest::plottingandformattingpp(bool iseta, bool isbk)
{

	TString chi2_title;
	TString chi2_saving_path;
	TString data_mc_title;
	TString data_data_title;
	TString data_mc_saving_path;
	TString data_data_saving_path;

	if (isbk && iseta)
	{
		chi2_title = "pp, |#eta| < 1.0 with bksub, Period: (%i)";
		data_mc_title = "eta_with_bksub_%i";
		data_data_title = "eta_%i";
		chi2_saving_path = "./chi2pp/chi2plots/eta/bksub/eta_with_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/eta/bksub/eta_with_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/eta/eta_%i.png";
	}
	if (!isbk && iseta)
	{
		chi2_title = "pp, |#eta| < 1.0 without bksub, Period: (%i)";
		data_mc_title = "eta_without_bksub_%i";
		data_data_title = "eta_%i";
		chi2_saving_path = "./chi2pp/chi2plots/eta/nobksub/eta_without_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/eta/nobksub/eta_without_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/eta/eta_%i.png";
	}
	if (!isbk && !iseta)
	{
		chi2_title = "pp, |#eta| < 2.4 without bksub, Period: (%i)";
		data_mc_title = "Raw_without_bksub_%i";
		data_data_title = "Raw_%i";
		chi2_saving_path = "./chi2pp/chi2plots/raw/nobksub/raw_without_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/raw/nobksub/raw_without_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/raw/raw_%i.png";
	}
	if (isbk && !iseta)
	{
		chi2_title = "pp, |#eta| < 2.4 with bksub, Period: (%i)";
		data_mc_title = "Raw_with_bksub_%i";
		data_data_title = "Raw_%i";
		chi2_saving_path = "./chi2pp/chi2plots/raw/bksub/raw_with_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/raw/bksub/raw_with_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/raw/raw_%i.png";
	}

	for (int runperiod = 0; runperiod < 22; runperiod++)
	{
		// Chi2/ndf plot
		c_2d_chisquare_ndf_pp[runperiod]->cd();
		gStyle->SetPalette(kRainBow);
		// c_2d_chisquare_ndf[cent]->SetLogz();
		h_chisquare_pp[runperiod]->SetTitle(Form(chi2_title, runperiod));
		h_chisquare_pp[runperiod]->SetTitleFont(42);
		h_chisquare_pp[runperiod]->Draw("COLZ");
		h_chisquare_pp[runperiod]->Draw("TEXTSAME");
		h_chisquare_pp[runperiod]->GetXaxis()->SetNdivisions(nbins_mass_shift, 0, 0);
		h_chisquare_pp[runperiod]->GetYaxis()->SetNdivisions(nbins_smear, 0, 0);
		h_chisquare_pp[runperiod]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare_pp[runperiod]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare_pp[runperiod]->GetXaxis()->SetTitle("Shifted Amount (GeV)");
		h_chisquare_pp[runperiod]->GetYaxis()->SetTitle("Smeared Amount (GeV)");

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

		Double_t xMin = h_chisquare_pp[runperiod]->GetXaxis()->GetBinLowEdge(minBinX);
		Double_t xMax = h_chisquare_pp[runperiod]->GetXaxis()->GetBinUpEdge(minBinX);
		Double_t yMin = h_chisquare_pp[runperiod]->GetYaxis()->GetBinLowEdge(minBinY);
		Double_t yMax = h_chisquare_pp[runperiod]->GetYaxis()->GetBinUpEdge(minBinY);

		Double_t xCenter = h_chisquare_pp[runperiod]->GetXaxis()->GetBinCenter(minBinX);
		Double_t yCenter = h_chisquare_pp[runperiod]->GetYaxis()->GetBinCenter(minBinY);

		dMass_pp[runperiod] = xCenter;
		dWidth_pp[runperiod] = yCenter;

		this->getcontour(h_chisquare_pp[runperiod], 2, 117, minBinX, minBinY, minContent, contour_x_left_onesig_HI, contour_x_right_onesig_HI);

		g_HI_contour_1sig_left = new TGraph(nbins_smear, contour_x_left_weighted_onesig, contour_y_HI[10]);
		g_HI_contour_1sig_right = new TGraph(nbins_smear, contour_x_right_weighted_onesig, contour_y_HI[10]);
		g_HI_contour_2sig_left = new TGraph(nbins_smear, contour_x_left_weighted_twosig, contour_y_HI[10]);
		g_HI_contour_2sig_right = new TGraph(nbins_smear, contour_x_right_weighted_twosig, contour_y_HI[10]);

		this->drawcontour(g_HI_contour_1sig_left, g_HI_contour_1sig_right, g_HI_contour_2sig_left, g_HI_contour_2sig_right, runperiod, isbk, iseta, true);

		for (int nbinssmearing = 1; nbinssmearing <= nbins_smear; ++nbinssmearing)
		{
			c_2d_chisquare_ndf_pp[runperiod]->cd();

			// cout << "contour_x_left_before_onesig_HI " << nbinssmearing - 1 << " is " << contour_x_left_before_onesig_HI[nbinssmearing-1] << endl;
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

		dMass_Err_pp[runperiod] = getuncertainty(h_chisquare_pp[runperiod], 1, minBinX, minBinY, 117);
		dWidth_Err_pp[runperiod] = getuncertainty(h_chisquare_pp[runperiod], 2, minBinX, minBinY, 117);

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
		box1->SetLineColor(kRed);
		box1->SetLineWidth(4);
		box1->SetFillStyle(0);
		box1->Draw("same");

		c_2d_chisquare_ndf_pp[runperiod]->SaveAs(Form(chi2_saving_path, runperiod));

		// This is Data and MC

		c_data_mc_raw_pp[runperiod]->cd();
		c_data_mc_raw_pp[runperiod]->SetLogy();
		if (isbk)
		{
			h_data_bksub_pp[runperiod]->SetTitle(Form(data_mc_title, runperiod));
			h_data_bksub_pp[runperiod]->SetMarkerColor(kRed);
			h_data_bksub_pp[runperiod]->SetMarkerStyle(kFullCircle);
			h_data_bksub_pp[runperiod]->Draw("P");
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerColor(kGreen);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->Draw("P SAME");
		}
		if (!isbk)
		{
			h_data_pp[runperiod]->SetTitle(Form(data_mc_title, runperiod));
			h_data_pp[runperiod]->SetMarkerColor(kRed);
			h_data_pp[runperiod]->SetMarkerStyle(kFullCircle);
			h_data_pp[runperiod]->Draw("P");
			// h_data_bksub_pp[runperiod]->SetMarkerColor(kBlue);
			// h_data_bksub_pp[runperiod]->SetMarkerStyle(kFullCircle);
			// h_data_bksub_pp[runperiod]->Draw("PSAME");
			// cout << "Missing Bin content is " << h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->GetBinContent(120) << endl;
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerColor(kGreen);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal_pp[minBinX - 1][minBinY - 1][runperiod]->Draw("P SAME");
		}

		TPaveText *pt1 = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
		pt1->AddText("Red is data, green is best signal");
		pt1->SetTextSize(0.04);
		pt1->SetTextAlign(22); // Center alignment
		pt1->Draw();

		c_data_mc_raw_pp[runperiod]->SaveAs(Form(data_mc_saving_path, runperiod));

		// This is data vs data_bksub
		c_data_data_bk_pp[runperiod]->cd();
		c_data_data_bk_pp[runperiod]->SetLogy();
		h_data_pp[runperiod]->SetTitle(Form(data_data_title, runperiod));
		h_data_pp[runperiod]->SetMarkerColor(kRed);
		h_data_pp[runperiod]->SetMarkerStyle(kFullCircle);
		h_data_pp[runperiod]->Draw("P HIST");
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
	if (isbk && iseta)
	{
		g_pp_dmass->Write("pp_dM_chi2_eta_bksub", 2);
		g_pp_dwidth->Write("pp_dWidth_chi2_eta_bksub", 2);
	}
	if (!isbk && iseta)
	{
		g_pp_dmass->Write("pp_dM_chi2_eta_nobksub", 2);
		g_pp_dwidth->Write("pp_dWidth_chi2_eta_nobksub", 2);
	}
	if (isbk && !iseta)
	{
		g_pp_dmass->Write("pp_dM_chi2_raw_bksub", 2);
		g_pp_dwidth->Write("pp_dWidth_chi2_raw_bksub", 2);
	}
	if (!isbk && !iseta)
	{
		g_pp_dmass->Write("pp_dM_chi2_raw_nobksub", 2);
		g_pp_dwidth->Write("pp_dWidth_chi2_raw_nobksub", 2);
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

void chisquaretest::RebinAll(int type)
{

	std::vector<double> newBinEdges_0_100 = createCustomBinning(h_data[0], 60, 80, 10, 80, 86, 4, 95, 100, 5, 100, 120, 40);

	int nNewBins = newBinEdges_0_100.size() - 1;

	Double_t *binEdgesArray = &newBinEdges_0_100[0];

	if (type == 1) // Now only have the version for eta cut
	{
		for (int cent = 0; cent < nbins_cent; cent++)
		{
			if (!((cent < 4) || (cent == 10)))
				continue;

			cout << "We are here" << cent << endl;
			cout << "nbinsX is " << h_data_bksub[cent]->GetNbinsX() << endl;

			h_data[cent] = (TH1D *)h_data[cent]->Rebin(nNewBins, Form("mass_array_data_witheta_rebin_%i", cent), binEdgesArray);

			h_data_bksub[cent] = (TH1D *)h_data_bksub[cent]->Rebin(nNewBins, Form("normalized_mc_bk_rebinned_%i", cent), binEdgesArray);
			cout << "222222" << endl;
			// h_data[cent]->Rebin(4);
			// h_data_bksub[cent]->Rebin(4);

			for (int shift = 0; shift < nbins_mass_shift; shift++)
			{
				for (int smear = 0; smear < nbins_smear; smear++)
				{

					// h_mc_signal[shift][smear][cent]->Rebin(4);
					h_mc_signal[shift][smear][cent] = (TH1D *)h_mc_signal[shift][smear][cent]->Rebin(nNewBins, Form("modifiedmass_eta_without_eff_%i_%i_%i_new", shift, smear, cent), binEdgesArray);
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

			h_data[cent]->Rebin(4);
			h_data_bksub[cent]->Rebin(4);

			for (int shift = 0; shift < nbins_mass_shift; shift++)
			{
				for (int smear = 0; smear < nbins_smear; smear++)
				{
					h_mc_signal[shift][smear][cent]->Rebin(4);
				}
			}
		}
	}
}

void chisquaretest::RebinAllpp(int x)
{
	for (int runperiod = 0; runperiod < 22; runperiod++)
	{
		h_data_pp[runperiod]->Rebin(x);
		h_data_bksub_pp[runperiod]->Rebin(x);

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				h_mc_signal_pp[shift][smear][runperiod]->Rebin(x);
			}
		}
	}
}

Double_t chisquaretest::getuncertainty(TH2D *h_1, int type, Int_t minBinX, Int_t minBinY, int ndf)
{
	double threshold = 2.3;

	if (type == 1)
	{
		// this is the dM

		double closest_difference = 1e9;
		int bin_closest = -1;
		for (int binx = 1; binx <= h_1->GetNbinsX(); ++binx)
		{
			if (binx == minBinX)
				continue;

			double chi2_bin = h_1->GetBinContent(binx, minBinY);
			double difference = fabs(chi2_bin - h_1->GetBinContent(minBinX, minBinY));
			// cout << "Lowest bin content is " << h_1->GetBinContent(minBinX, minBinY) << endl;
			// cout << "For bin " << binx << "Difference is " << difference << endl;
			if (fabs(difference - threshold) < fabs(closest_difference - threshold))
			{
				closest_difference = difference;
				bin_closest = binx; // Store the bin number
			}
		}
		if (bin_closest != -1)
		{
			double chi2_closest = h_1->GetBinContent(bin_closest, minBinY);
			std::cout << "X Bin with chi2/ndf closest to 1 unit difference from minimum: " << bin_closest << std::endl;
			std::cout << "X Chi2/ndf value in that bin: " << chi2_closest << std::endl;
			std::cout << "X Difference with minimum chi2/ndf: " << closest_difference << std::endl;
			return (fabs(h_1->GetXaxis()->GetBinCenter(bin_closest) - h_1->GetXaxis()->GetBinCenter(minBinX)));
		}
		else
		{
			std::cout << "No bin found with a difference close to 1" << std::endl;
			return -99;
		}
	}

	else if (type == 2)
	{
		// this is the dW

		double closest_difference = 1e9;
		int bin_closest = -1;
		for (int biny = 1; biny <= h_1->GetNbinsY(); ++biny)
		{
			if (biny == minBinY)
				continue;
			double chi2_bin = h_1->GetBinContent(minBinX, biny);
			double difference = fabs(chi2_bin - h_1->GetBinContent(minBinX, minBinY));
			if (fabs(difference - threshold) < fabs(closest_difference - threshold))
			{
				closest_difference = difference;
				bin_closest = biny; // Store the bin number
			}
		}
		if (bin_closest != -1)
		{
			double chi2_closest = h_1->GetBinContent(minBinX, bin_closest);
			// std::cout << "Y Bin with chi2/ndf closest to 1 unit difference from minimum: " << bin_closest << std::endl;
			// std::cout << "Y Chi2/ndf value in that bin: " << chi2_closest << std::endl;
			// std::cout << "Y Difference with minimum chi2/ndf: " << closest_difference << std::endl;
			// double sigmasmear = 91.1876 * h_1->GetYaxis()->GetBinCenter(minBinY);
			// double errofsmear = fabs(h_1->GetYaxis()->GetBinCenter(bin_closest) - h_1->GetYaxis()->GetBinCenter(minBinY)) * 91.1876;
			// double deno = sqrt(TMath::Power(2.4955, 2) + TMath::Power(h_1->GetYaxis()->GetBinCenter(minBinY) * 91.1876, 2));
			// double errorofdwidth = (sigmasmear / deno) * errofsmear;
			return (fabs(h_1->GetYaxis()->GetBinCenter(bin_closest) - h_1->GetYaxis()->GetBinCenter(minBinY)));
		}
		else
		{
			std::cout << "No bin found with a difference close to 1" << std::endl;
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

	// cout << "NBin for data is " << nbinsh1 << endl;

	if (nbinsh1 != nbinsh2)
	{
		std::cerr << "Warning: h1 and h2 have different bin numbers!" << std::endl;
		return -99; // Optionally, return an error value
	}

	Double_t chi2 = 0;
	for (int i = 1; i <= nbinsh1; i++)
	{
		double num = TMath::Power((h1->GetBinContent(i) - h2->GetBinContent(i)), 2);
		double den = TMath::Power((h1->GetBinError(i)), 2) + TMath::Power((h2->GetBinError(i)), 2);

		if (TMath::Power((h1->GetBinContent(i)), 1) <= 0)
		{
			std::cerr << "Warning: Bin " << i << " has zero errors/missing bins/ missing bins after bksub for data. Skipping this bin." << std::endl;
			continue; // Skip this bin to avoid division by zero
		}

		chi2 += num / den;
	}

	Double_t ndf = nbinsh1 - 2;

	return chi2;
}

void chisquaretest::getcontour(TH2D *h1, int type, int ndf, Int_t minBinX, Int_t minBinY, Double_t minBinContent, Double_t *arrayleft, Double_t *arrayright)
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

	if (ndf == -99)
		cout << "You should never see this, meaning we are running with undefined ndf" << endl;

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
	double range4_min, double range4_max, int rebin4)
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
		double nOriginalBins = (range_max - range_min) / binWidth;
		// cout << std::fmod(nOriginalBins, rebin) << endl;
		return std::fmod(nOriginalBins, rebin) == 0; // Proper if nOriginalBins is divisible by rebin
	};

	// Check if the rebinning factors are proper for each range
	if (!is_proper_rebin(range1_min, range1_max, binWidth, rebin1))
	{
		std::cerr << "Warning: Rebinning factor for range [" << range1_min << ", " << range1_max << "] is improper." << std::endl;
	}
	if (!is_proper_rebin(range2_min, range2_max, binWidth, rebin2))
	{
		std::cerr << "Warning: Rebinning factor for range [" << range2_min << ", " << range2_max << "] is improper." << std::endl;
	}
	if (!is_proper_rebin(range3_min, range3_max, binWidth, rebin3))
	{
		std::cerr << "Warning: Rebinning factor for range [" << range3_min << ", " << range3_max << "] is improper." << std::endl;
	}
	if (!is_proper_rebin(range4_min, range4_max, binWidth, rebin4))
	{
		std::cerr << "Warning: Rebinning factor for range [" << range4_min << ", " << range4_max << "] is improper." << std::endl;
	}

	// Step 1: Loop over the histogram bins and create bin edges
	double currentEdge = xMin;

	while (currentEdge < xMax)
	{
		// If we are in the first rebinning range
		if (currentEdge >= range1_min && currentEdge < range1_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin1 * binWidth;
		}
		// If we are in the second rebinning range
		else if (currentEdge >= range2_min && currentEdge < range2_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin2 * binWidth;
		}
		// If we are in the third rebinning range
		else if (currentEdge >= range3_min && currentEdge < range3_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin3 * binWidth;
		}
		// If we are in the fourth rebinning range
		else if (currentEdge >= range4_min && currentEdge < range4_max)
		{
			newBins.push_back(currentEdge);
			currentEdge += rebin4 * binWidth;
		}
		// If we are outside the custom rebinning ranges
		else
		{
			newBins.push_back(currentEdge);
			currentEdge += binWidth; // Regular bin width
		}
	}

	newBins.push_back(xMax); // Make sure to include the upper edge of the last bin

	return newBins;
}

void chisquaretest::drawcontour(TGraph *onesig_left, TGraph *onesig_right, TGraph *twosig_left, TGraph *twosig_right, int iteration, bool isbk, bool iseta, bool ispp)
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

	TString title = "";

	if (iseta && !ispp)
	{
		onesig_left->GetXaxis()->SetLimits(eta_mass_shift_array_low[iteration], eta_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetLimits(eta_mass_smear_array_low[iteration], eta_mass_smear_array_high[iteration]);
		onesig_left->GetXaxis()->SetRangeUser(eta_mass_shift_array_low[iteration], eta_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetRangeUser(eta_mass_smear_array_low[iteration], eta_mass_smear_array_high[iteration]);
		if (isbk)
			title = Form("PbPb Contour plot with bksub, |#eta| < 1.0, Cent: (%i-%i)", cenlowlimit[iteration], cenhighlimit[iteration]);
		if (!isbk)
			title = Form("PbPb Contour plot without bksub, |#eta| < 1.0, Cent: (%i-%i)", cenlowlimit[iteration], cenhighlimit[iteration]);
	}
	if (!iseta && !ispp)
	{
		onesig_left->GetXaxis()->SetLimits(raw_mass_shift_array_low[iteration], raw_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetLimits(raw_mass_smear_array_low[iteration], raw_mass_smear_array_high[iteration]);
		onesig_left->GetXaxis()->SetRangeUser(raw_mass_shift_array_low[iteration], raw_mass_shift_array_high[iteration]);
		onesig_left->GetYaxis()->SetRangeUser(raw_mass_smear_array_low[iteration], raw_mass_smear_array_high[iteration]);
		if (isbk)
			title = Form("PbPb Contour plot with bksub, |#eta| < 2.4, Cent: (%i-%i)", cenlowlimit[iteration], cenhighlimit[iteration]);
		if (!isbk)
			title = Form("PbPb Contour plot without bksub, |#eta| < 2.4, Cent: (%i-%i)", cenlowlimit[iteration], cenhighlimit[iteration]);
	}
	if (ispp)
	{
		if (iseta)
		{
			onesig_left->GetXaxis()->SetLimits(eta_pp_mass_shift_low, eta_pp_mass_shift_high);
			onesig_left->GetYaxis()->SetLimits(eta_pp_smear_low, eta_pp_smear_high);
			onesig_left->GetXaxis()->SetRangeUser(eta_pp_mass_shift_low, eta_pp_mass_shift_high);
			onesig_left->GetYaxis()->SetRangeUser(eta_pp_smear_low, eta_pp_smear_high);
			if (isbk)
				title = Form("pp Contour plot with bksub, |#eta| < 1.0, Period: (%i)", iteration);
			if (!isbk)
				title = Form("pp Contour plot without bksub, |#eta| < 1.0, Period: (%i)", iteration);
		}
		if (!iseta)
		{
			onesig_left->GetXaxis()->SetLimits(raw_pp_mass_shift_low, raw_pp_mass_shift_high);
			onesig_left->GetYaxis()->SetLimits(raw_pp_smear_low, raw_pp_smear_high);
			onesig_left->GetXaxis()->SetRangeUser(raw_pp_mass_shift_low, raw_pp_mass_shift_high);
			onesig_left->GetYaxis()->SetRangeUser(raw_pp_smear_low, raw_pp_smear_high);
			if (isbk)
				title = Form("pp Contour plot with bksub, |#eta| < 2.4, Period: (%i)", iteration);
			if (!isbk)
				title = Form("pp Contour plot without bksub, |#eta| < 2.4, Period: (%i)", iteration);
		}
	}
	// std::cout << "Title: " << title << std::endl;

	TH1 *frame = onesig_left->GetHistogram(); // Get the underlying histogram for customization

	onesig_left->SetTitle(title);

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
		//std::cout << "Point " << i << ": (" << x << ", " << y << ")" << std::endl;
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

	if (isbk && iseta == 1)
	{
		if (ispp)
			temp_c1->SaveAs(Form("./chi2pp/contour/eta/contour_eta_bksub_%i.png", iteration));
		if (!ispp)
			temp_c1->SaveAs(Form("./newchi2/contour/eta/contour_eta_bksub_%i_%i.png", this->cenlowlimit[iteration], this->cenhighlimit[iteration]));
	}
	if (!isbk && iseta == 1)
	{
		if (ispp)
			temp_c1->SaveAs(Form("./chi2pp/contour/eta/contour_eta_%i.png", iteration));
		if (!ispp)
			temp_c1->SaveAs(Form("./newchi2/contour/eta/contour_eta_%i_%i.png", this->cenlowlimit[iteration], this->cenhighlimit[iteration]));
	}
	if (isbk && iseta == 0)
	{
		if (ispp)
			temp_c1->SaveAs(Form("./chi2pp/contour/raw/contour_raw_bksub_%i.png", iteration));
		if (!ispp)
			temp_c1->SaveAs(Form("./newchi2/contour/raw/contour_raw_bksub_%i_%i.png", this->cenlowlimit[iteration], this->cenhighlimit[iteration]));
	}
	if (!isbk && iseta == 0)
	{
		if (ispp)
			temp_c1->SaveAs(Form("./chi2pp/contour/raw/contour_raw_%i.png", iteration));
		if (!ispp)
			temp_c1->SaveAs(Form("./newchi2/contour/raw/contour_raw_%i_%i.png", this->cenlowlimit[iteration], this->cenhighlimit[iteration]));
	}
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