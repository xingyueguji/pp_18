#include "plotting_helper.h"
class chisquaretest : public plotting_helper
{
public:
	chisquaretest(TString s1, TString s2, TString s3, int type);
	chisquaretest(TString s1, TString s2, TString s3, bool iseta);
	~chisquaretest();
	void calculatechisq(bool isbk);
	void calculatechisqpp(bool isbk);
	void plottingandformatting(int type, bool isbk);
	void plottingandformattingpp(bool iseta, bool isbk);
	void bincontentcheck(bool isbk);
	void RebinAll(int x);
	void RebinAllpp(int x);
	Double_t getuncertainty(TH2D *h_1, int type, Int_t minBinX, Int_t minBinY);

	static const int nbins_mass_shift = 26;
	static const int nbins_smear = 26;
	static const int nbins_cent = 11;
	static constexpr double lowbin_mass_shift = -0.5;
	static constexpr double highbin_mass_shift = 0.2;
	static constexpr double lowbin_smear = 0;
	static constexpr double highbin_smear = 0.02;

	double h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;
	double h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (nbins_mass_shift - 1)) / 2;

	double h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;
	double h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (nbins_smear - 1)) / 2;

	TH1D *h_mc_signal[nbins_mass_shift][nbins_smear][nbins_cent];
	TH1D *h_data[nbins_cent];
	TH1D *h_mc_bk[nbins_cent];
	TH1D *h_data_bksub[nbins_cent];
	TH2D *h_chisquare[nbins_cent];

	Double_t dMass_HI[6];
	Double_t dMass_Err_HI[6];

	Double_t dWidth_HI[6];
	Double_t dWidth_Err_HI[6];

	Double_t xposition_HI[6];
	Double_t xposition_err_HI[6];

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

	double chisquarearray[nbins_mass_shift][nbins_smear][nbins_cent] = {};

	TString mcfilepath;
	TString datafilepath;
	TString bkfilepath;

	TCanvas *c_2d_chisquare_ndf[nbins_cent];
	// What I need: invariant mass comparison data vs mc raw, data vs data bk sub.
	TCanvas *c_data_mc_raw[nbins_cent];
	TCanvas *c_data_data_bk[nbins_cent];

	TCanvas *c_2d_chisquare_ndf_pp[22];
	TCanvas *c_data_mc_raw_pp[22];
	TCanvas *c_data_data_bk_pp[22];

	TFile *mcfile;
	TFile *datafile;
	TFile *bkfile;
};

chisquaretest::chisquaretest(TString s1, TString s2, TString s3, int type)
{

	mcfilepath = s1;
	datafilepath = s2;
	bkfilepath = s3;

	mcfile = new TFile(mcfilepath, "READ");
	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath, "READ");

	for (int i = 0; i < 6; i++)
	{
		xposition_HI[i] = i + 1;
		xposition_err_HI[i] = 0;
	}

	for (int cent = 0; cent < nbins_cent; cent++)
	{

		c_2d_chisquare_ndf[cent] = new TCanvas(Form("c_2d_chisquare_ndf_%i", cent), "", 1600, 1200);
		c_data_mc_raw[cent] = new TCanvas(Form("c_data_mc_raw_%i", cent), "", 800, 600);
		c_data_data_bk[cent] = new TCanvas(Form("c_data_data_bk_%i", cent), "", 800, 600);

		if (type == 0)
			h_data[cent] = (TH1D *)datafile->Get(Form("mass_array_data_%i", cent));
		if (type == 1)
			h_data[cent] = (TH1D *)datafile->Get(Form("mass_array_data_witheta_%i", cent));

		h_mc_bk[cent] = (TH1D *)bkfile->Get(Form("Normalized_mc_bk_%i", cent));
		h_chisquare[cent] = new TH2D(Form("h_chisquare_%i", cent), Form("Cent_%i_%i", this->cenlowlimit[cent], this->cenhighlimit[cent]), nbins_mass_shift,h_low_mass_shift, h_high_mass_shift, nbins_smear, h_low_smear, h_high_smear); // Equation here: half bin to the left and half bin to the right, bin width = range / (21-1)

		this->areanormalize(h_data[cent]);

		h_data_bksub[cent] = (TH1D *)h_data[cent]->Clone();
		h_data_bksub[cent]->Add(h_mc_bk[cent], -1);

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				if (type == 0)
					h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_raw_without_eff_%i_%i_%i", shift, smear, cent));
				if (type == 1)
					h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_eta_without_eff_%i_%i_%i", shift, smear, cent));
				this->areanormalize(h_mc_signal[shift][smear][cent]);
			}
		}
	}
}

chisquaretest::chisquaretest(TString s1, TString s2, TString s3, bool iseta)
{
	mcfilepath = s1;
	datafilepath = s2;
	bkfilepath = s3;

	mcfile = new TFile(mcfilepath, "READ");
	datafile = new TFile(datafilepath, "READ");
	bkfile = new TFile(bkfilepath, "READ");

	h_mc_bk_pp = (TH1D *)bkfile->Get("Normalized_mc_bk_10");

	for (int runperiod = 0; runperiod < 22; runperiod++)
	{
		c_2d_chisquare_ndf_pp[runperiod] = new TCanvas(Form("c_2d_chisquare_ndf_pp_%i", runperiod), "", 1600, 1200);
		c_data_mc_raw_pp[runperiod] = new TCanvas(Form("c_data_mc_raw_pp_%i", runperiod), "", 800, 600);
		c_data_data_bk_pp[runperiod] = new TCanvas(Form("c_data_data_bk_pp_%i", runperiod), "", 800, 600);

		if (iseta)
		{
			cout << "WE ARE HERE" << endl;
			h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("h_mass_array_eta_%i", runperiod));
		}
		else
		{
			h_data_pp[runperiod] = (TH1D *)datafile->Get(Form("h_mass_array_raw_%i", runperiod));
		}

		h_chisquare_pp[runperiod] = new TH2D(Form("h_chisquare_pp_%i", runperiod), Form("h_chisquare_pp_%i", runperiod), nbins_mass_shift, -0.5175, 0.2175, nbins_smear, -0.0005, 0.0205);
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
					h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_eta_without_eff_%i_%i_%i", shift, smear, 10));
					//h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_%i_%i_%i", shift, smear, 10));
				}
				else
				{
					h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_raw_without_eff_%i_%i_%i", shift, smear, 10));
					//h_mc_signal_pp[shift][smear][runperiod] = (TH1D *)mcfile->Get(Form("modifiedmass_%i_%i_%i", shift, smear, 10));
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
		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				double chisquarevalue = 0;
				cout << "Now running File " << cent << " mass shift " << shift << " smearing " << smear << endl;
				if (isbk)
					chisquarevalue = h_data_bksub[cent]->Chi2Test(h_mc_signal[shift][smear][cent], "WW CHI2/NDF");
				if (!isbk)
					chisquarevalue = h_data[cent]->Chi2Test(h_mc_signal[shift][smear][cent], "WW CHI2/NDF");
				// cout << chisquarevalue << endl;
				// chisquarearray[shift][smear][cent] = chisquarevalue;
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
				cout << "Now running File " << runperiod << " mass shift " << shift << " smearing " << smear << endl;
				if (isbk)
					chisquarevalue = h_data_bksub_pp[runperiod]->Chi2Test(h_mc_signal_pp[shift][smear][runperiod], "WW CHI2/NDF");
				if (!isbk)
					chisquarevalue = h_data_pp[runperiod]->Chi2Test(h_mc_signal_pp[shift][smear][runperiod], "WW CHI2/NDF");
				// cout << chisquarevalue << endl;
				// chisquarearray[shift][smear][cent] = chisquarevalue;
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

	if (type == 0)
	{
		if (isbk)
		{
			chi2_title = "Raw_with_bksub_%i_%i";
			data_mc_title = "Raw_with_bksub_%i_%i";
			data_data_title = "Raw_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/raw/bksub/Raw_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/raw/bksub/Raw_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_%i_%i.png";
		}
		if (!isbk)
		{
			chi2_title = "Raw_without_bksub_%i_%i";
			data_mc_title = "Raw_without_bksub_%i_%i";
			data_data_title = "Raw_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/raw/nobksub/Raw_without_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/raw/nobksub/Raw_without_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/raw/Raw_%i_%i.png";
		}
	}
	if (type == 1)
	{
		if (isbk)
		{
			chi2_title = "Eta_with_bksub_%i_%i";
			data_mc_title = "Eta_with_bksub_%i_%i";
			data_data_title = "Eta_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/eta/bksub/Eta_with_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/eta/bksub/Eta_with_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_%i_%i.png";
		}
		if (!isbk)
		{
			chi2_title = "Eta_without_bksub_%i_%i";
			data_mc_title = "Eta_without_bksub_%i_%i";
			data_data_title = "Eta_%i_%i";
			chi2_saving_path = "./newchi2/chi2plots/eta/nobksub/Eta_without_bksub_%i_%i.png";
			chi2_saving_path = "./newchi2/chi2plots/eta/nobksub/Eta_without_bksub_%i_%i.png";
			data_mc_saving_path = "./newchi2/datamc/eta/nobksub/Eta_without_bksub_%i_%i.png";
			data_data_saving_path = "./newchi2/datadata/eta/Eta_%i_%i.png";
		}
	}

	for (int cent = 0; cent < nbins_cent; cent++)
	{
		// Chi2/ndf plot
		c_2d_chisquare_ndf[cent]->cd();
		gStyle->SetPalette(kRainBow);
		// c_2d_chisquare_ndf[cent]->SetLogz();
		h_chisquare[cent]->SetTitle(Form(chi2_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		h_chisquare[cent]->Draw("COLZ");
		h_chisquare[cent]->Draw("TEXTSAME");
		h_chisquare[cent]->GetXaxis()->SetNdivisions(nbins_mass_shift, 0, 0);
		h_chisquare[cent]->GetYaxis()->SetNdivisions(nbins_smear, 0, 0);
		h_chisquare[cent]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare[cent]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare[cent]->GetXaxis()->SetTitle("Shifted Amount (GeV)");
		h_chisquare[cent]->GetYaxis()->SetTitle("Smeared Amount (Sig)");

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

		Double_t xMin = h_chisquare[cent]->GetXaxis()->GetBinLowEdge(minBinX);
		Double_t xMax = h_chisquare[cent]->GetXaxis()->GetBinUpEdge(minBinX);
		Double_t yMin = h_chisquare[cent]->GetYaxis()->GetBinLowEdge(minBinY);
		Double_t yMax = h_chisquare[cent]->GetYaxis()->GetBinUpEdge(minBinY);

		Double_t xCenter = h_chisquare[cent]->GetXaxis()->GetBinCenter(minBinX);
		Double_t yCenter = h_chisquare[cent]->GetYaxis()->GetBinCenter(minBinY);
		yCenter = 91.1876 * (yCenter);

		if (cent == 0 || cent == 1 || cent == 2 || cent == 3 || cent == 4)
		{
			dMass_HI[cent + 1] = xCenter;
			dWidth_HI[cent + 1] = yCenter;

			dMass_Err_HI[cent + 1] = getuncertainty(h_chisquare[cent], 1, minBinX, minBinY);
			dWidth_Err_HI[cent + 1] = getuncertainty(h_chisquare[cent], 2, minBinX, minBinY);
		}

		if (cent == 10)
		{
			dMass_HI[0] = xCenter;
			dWidth_HI[0] = yCenter;

			dMass_Err_HI[0] = getuncertainty(h_chisquare[cent], 1, minBinX, minBinY);
			dWidth_Err_HI[0] = getuncertainty(h_chisquare[cent], 2, minBinX, minBinY);
		}

		TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
		box1->SetLineColor(kRed);
		box1->SetLineWidth(4);
		box1->SetFillStyle(0);
		box1->Draw("same");

		c_2d_chisquare_ndf[cent]->SaveAs(Form(chi2_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));

		// This is Data and MC

		c_data_mc_raw[cent]->cd();
		c_data_mc_raw[cent]->SetLogy();
		if (isbk)
		{
			h_data_bksub[cent]->SetTitle(Form(data_mc_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
			h_data_bksub[cent]->SetMarkerColor(kRed);
			h_data_bksub[cent]->SetMarkerStyle(kFullCircle);
			h_data_bksub[cent]->Draw("P");
			h_mc_signal[minBinX-1][minBinY-1][cent]->SetMarkerColor(kGreen);
			h_mc_signal[minBinX-1][minBinY-1][cent]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal[minBinX-1][minBinY-1][cent]->Draw("P SAME");
		}
		if (!isbk)
		{
			h_data[cent]->SetTitle(Form(data_mc_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
			h_data[cent]->SetMarkerColor(kRed);
			h_data[cent]->SetMarkerStyle(kFullCircle);
			h_data[cent]->Draw("P");
			h_mc_signal[minBinX-1][minBinY-1][cent]->SetMarkerColor(kGreen);
			h_mc_signal[minBinX-1][minBinY-1][cent]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal[minBinX-1][minBinY-1][cent]->Draw("P SAME");
		}

		TPaveText *pt1 = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
		pt1->AddText("Red is data, green is best signal");
		pt1->SetTextSize(0.04);
		pt1->SetTextAlign(22); // Center alignment
		pt1->Draw();

		c_data_mc_raw[cent]->SaveAs(Form(data_mc_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));

		// This is data vs data_bksub
		c_data_data_bk[cent]->cd();
		c_data_data_bk[cent]->SetLogy();
		h_data[cent]->SetTitle(Form(data_data_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
		h_data[cent]->SetMarkerColor(kRed);
		h_data[cent]->SetMarkerStyle(kFullCircle);
		h_data[cent]->Draw("P");
		h_data_bksub[cent]->SetMarkerColor(kGreen);
		h_data_bksub[cent]->SetMarkerStyle(kFullDotLarge);
		h_data_bksub[cent]->Draw("P SAME");
		h_mc_bk[cent]->SetMarkerColor(kBlue);
		//h_mc_bk[cent]->Draw("P SAME");

		TPaveText *pt = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
		pt->AddText("Red is data, green is data - bk, blue is bk");
		pt->SetTextSize(0.04);
		pt->SetTextAlign(22); // Center alignment
		pt->Draw();

		c_data_data_bk[cent]->SaveAs(Form(data_data_saving_path, this->cenlowlimit[cent], this->cenhighlimit[cent]));
	}
	g_HI_dmass = new TGraphErrors(6, xposition_HI, dMass_HI, xposition_err_HI, dMass_Err_HI);
	g_HI_dwidth = new TGraphErrors(6, xposition_HI, dWidth_HI, xposition_err_HI, dWidth_Err_HI);

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
		chi2_title = "eta_with_bksub_%i";
		data_mc_title = "eta_with_bksub_%i";
		data_data_title = "eta_%i";
		chi2_saving_path = "./chi2pp/chi2plots/eta/bksub/eta_with_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/eta/bksub/eta_with_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/eta/eta_%i.png";
	}
	if (!isbk && iseta)
	{
		chi2_title = "eta_without_bksub_%i";
		data_mc_title = "eta_without_bksub_%i";
		data_data_title = "eta_%i";
		chi2_saving_path = "./chi2pp/chi2plots/eta/nobksub/eta_without_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/eta/nobksub/eta_without_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/eta/eta_%i.png";
	}
	if (!isbk && !iseta)
	{
		chi2_title = "Raw_without_bksub_%i";
		data_mc_title = "Raw_without_bksub_%i";
		data_data_title = "Raw_%i";
		chi2_saving_path = "./chi2pp/chi2plots/raw/nobksub/raw_without_bksub_%i.png";
		data_mc_saving_path = "./chi2pp/datamc/raw/nobksub/raw_without_bksub_%i.png";
		data_data_saving_path = "./chi2pp/datadata/raw/raw_%i.png";
	}
	if (isbk && !iseta)
	{
		chi2_title = "Raw_with_bksub_%i";
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
		h_chisquare_pp[runperiod]->Draw("COLZ");
		h_chisquare_pp[runperiod]->Draw("TEXTSAME");
		h_chisquare_pp[runperiod]->GetXaxis()->SetNdivisions(nbins_mass_shift, 0, 0);
		h_chisquare_pp[runperiod]->GetYaxis()->SetNdivisions(nbins_smear, 0, 0);
		h_chisquare_pp[runperiod]->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
		h_chisquare_pp[runperiod]->GetYaxis()->SetLabelSize(0.02);
		h_chisquare_pp[runperiod]->GetXaxis()->SetTitle("Shifted Amount (GeV)");
		h_chisquare_pp[runperiod]->GetYaxis()->SetTitle("Smeared Amount (Sig)");

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
		yCenter = 91.1876 * (yCenter);

		dMass_pp[runperiod] = xCenter;
		dWidth_pp[runperiod] = yCenter;

		dMass_Err_pp[runperiod] = getuncertainty(h_chisquare_pp[runperiod], 1, minBinX, minBinY);
		dWidth_Err_pp[runperiod] = getuncertainty(h_chisquare_pp[runperiod], 2, minBinX, minBinY);

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
			h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->SetMarkerColor(kGreen);
			h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->Draw("P SAME");
		}
		if (!isbk)
		{
			h_data_pp[runperiod]->SetTitle(Form(data_mc_title, runperiod));
			h_data_pp[runperiod]->SetMarkerColor(kRed);
			h_data_pp[runperiod]->SetMarkerStyle(kFullCircle);
			h_data_pp[runperiod]->Draw("P");
			//h_data_bksub_pp[runperiod]->SetMarkerColor(kBlue);
			//h_data_bksub_pp[runperiod]->SetMarkerStyle(kFullCircle);
			//h_data_bksub_pp[runperiod]->Draw("PSAME");
			//cout << "Missing Bin content is " << h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->GetBinContent(120) << endl;
			h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->SetMarkerColor(kGreen);
			h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->SetMarkerStyle(kFullDotLarge);
			h_mc_signal_pp[minBinX-1][minBinY-1][runperiod]->Draw("P SAME");
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

void chisquaretest::RebinAll(int x)
{
	for (int cent = 0; cent < nbins_cent; cent++)
	{
		h_data[cent]->Rebin(x);
		h_data_bksub[cent]->Rebin(x);

		for (int shift = 0; shift < nbins_mass_shift; shift++)
		{
			for (int smear = 0; smear < nbins_smear; smear++)
			{
				h_mc_signal[shift][smear][cent]->Rebin(x);
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

Double_t chisquaretest::getuncertainty(TH2D *h_1, int type, Int_t minBinX, Int_t minBinY)
{
	if (type == 1)
	{
		// this is the dM

		double closest_difference = 1e9;
		int bin_closest = -1;
		for (int binx = 1; binx <= h_1->GetNbinsX(); ++binx)
		{
			double chi2_bin = h_1->GetBinContent(binx, minBinY);
			double difference = fabs(chi2_bin - h_1->GetBinContent(minBinX, minBinY));
			cout << "Lowest bin content is " << h_1->GetBinContent(minBinX, minBinY) << endl;
			cout << "For bin " << binx << "Difference is " << difference << endl;
			if (fabs(difference - 1) < fabs(closest_difference - 1))
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
		// this is the dM

		double closest_difference = 1e9;
		int bin_closest = -1;
		for (int biny = 1; biny <= h_1->GetNbinsY(); ++biny)
		{
			double chi2_bin = h_1->GetBinContent(minBinX, biny);
			double difference = fabs(chi2_bin - h_1->GetBinContent(minBinX, minBinY));
			if (fabs(difference - 1) < fabs(closest_difference - 1))
			{
				closest_difference = difference;
				bin_closest = biny; // Store the bin number
			}
		}
		if (bin_closest != -1)
		{
			double chi2_closest = h_1->GetBinContent(minBinX, bin_closest);
			std::cout << "Y Bin with chi2/ndf closest to 1 unit difference from minimum: " << bin_closest << std::endl;
			std::cout << "Y Chi2/ndf value in that bin: " << chi2_closest << std::endl;
			std::cout << "Y Difference with minimum chi2/ndf: " << closest_difference << std::endl;
			return (fabs(h_1->GetYaxis()->GetBinCenter(bin_closest) - h_1->GetYaxis()->GetBinCenter(minBinY)) * 91.1876);
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
