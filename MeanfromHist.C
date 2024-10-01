#include "fit.h"
void GetHistpp(TFile *f1)
{
	for (int i = 0; i < 22; i++)
	{
		h_mass_array_raw[i] = (TH1D *)f1->Get(Form("h_mass_array_raw_%i", i));
		h_mass_array_eta[i] = (TH1D *)f1->Get(Form("h_mass_array_eta_%i", i));
	}
}
void GetHistHI(TFile *f1)
{
	h_HI_mass_array_raw[0] = (TH1D *)f1->Get("mass_array_data_10");
	h_HI_mass_array_eta[0] = (TH1D *)f1->Get("mass_array_data_witheta_10");

	for (int i = 1; i < 6; i++)
	{
		h_HI_mass_array_raw[i] = (TH1D *)f1->Get(Form("mass_array_data_%i", i - 1));
		h_HI_mass_array_eta[i] = (TH1D *)f1->Get(Form("mass_array_data_witheta_%i", i - 1));
	}
}
void fitpp(bool isbksub)
{
	for (int i = 0; i < 22; i++)
	{
		xposition[i] = i + 1;
		xposition_err[i] = 0;
		bksub_hist_pp_raw[i] = (TH1D *)h_mass_array_raw[i]->Clone();
		bksub_hist_pp_eta[i] = (TH1D *)h_mass_array_eta[i]->Clone();

		if (!isbksub)
		{
			TH1D *bksub_hist_pp_raw_aftercut = new TH1D("bksub_hist_pp_raw_aftercut", "", 40, 80, 100);
			TH1D *bksub_hist_pp_eta_aftercut = new TH1D("bksub_hist_pp_eta_aftercut", "", 40, 80, 100);

			int bin_low_raw = bksub_hist_pp_raw[i]->FindBin(80);
			int bin_high_raw = bksub_hist_pp_raw[i]->FindBin(100);

			int bin_low_eta = bksub_hist_pp_eta[i]->FindBin(80);
			int bin_high_eta = bksub_hist_pp_eta[i]->FindBin(100);

			for (int j = bin_low_raw; j <= bin_high_raw; j++)
			{
				double content_pp_raw = bksub_hist_pp_raw[i]->GetBinContent(j);
				double error_pp_raw = bksub_hist_pp_raw[i]->GetBinError(j);
				double content_pp_eta = bksub_hist_pp_eta[i]->GetBinContent(j);
				double error_pp_eta = bksub_hist_pp_eta[i]->GetBinError(j);

				int new_bin_raw = bksub_hist_pp_raw_aftercut->FindBin(bksub_hist_pp_raw[i]->GetBinCenter(j));
				int new_bin_eta = bksub_hist_pp_eta_aftercut->FindBin(bksub_hist_pp_eta[i]->GetBinCenter(j));

				bksub_hist_pp_raw_aftercut->SetBinContent(new_bin_raw, content_pp_raw);
				bksub_hist_pp_raw_aftercut->SetBinError(new_bin_raw, error_pp_raw);
				bksub_hist_pp_eta_aftercut->SetBinContent(new_bin_eta, content_pp_eta);
				bksub_hist_pp_eta_aftercut->SetBinError(new_bin_eta, error_pp_eta);
			}

			cout << "pp_raw_before_cut is " << bksub_hist_pp_raw[i]->GetMean() << " After cut is " << bksub_hist_pp_raw_aftercut->GetMean() << endl;

			dMass_raw[i] = (bksub_hist_pp_raw_aftercut->GetMean()) - 91.1876;
			dMass_Err_raw[i] = (bksub_hist_pp_raw_aftercut->GetMeanError());
			dWidth_raw[i] = (bksub_hist_pp_raw_aftercut->GetStdDev()) - 2.4955;
			dWidth_Err_raw[i] = (bksub_hist_pp_raw_aftercut->GetStdDevError());

			dMass_eta[i] = (bksub_hist_pp_eta_aftercut->GetMean()) - 91.1876;
			dMass_Err_eta[i] = (bksub_hist_pp_eta_aftercut->GetMeanError());
			dWidth_eta[i] = (bksub_hist_pp_eta_aftercut->GetStdDev()) - 2.4955;
			dWidth_Err_eta[i] = (bksub_hist_pp_eta_aftercut->GetStdDevError());

			delete bksub_hist_pp_raw_aftercut;
			delete bksub_hist_pp_eta_aftercut;
		}
		else
		{
			TCanvas *c1 = new TCanvas("c1", "", 800, 600);
			TCanvas *c2 = new TCanvas("c2", "", 800, 600);

			dataHistarray_raw[i] = new RooDataHist("dataHistarray_raw", "RooDataHist from TH1D", RooArgList(*x), bksub_hist_pp_raw[i]);
			dataHistarray_eta[i] = new RooDataHist("dataHistarray_eta", "RooDataHist from TH1D", RooArgList(*x), bksub_hist_pp_eta[i]);
			double Meanbefore_raw = bksub_hist_pp_raw[i]->GetMean();
			double Meanbefore_eta = bksub_hist_pp_eta[i]->GetMean();

			c1->cd();
			bksub_hist_pp_raw[i]->Draw("HIST");
			c2->cd();
			bksub_hist_pp_eta[i]->Draw("HIST");

			RooFitResult *fitResult_raw = expExtended->fitTo(*dataHistarray_raw[i], RooFit::Range("LeftTail,RightTail"), RooFit::Save());
			fitHist_pp_raw[i] = expExtended->createHistogram(Form("fitHist_pp_raw_%i", i), *x, RooFit::Binning(120, 60, 120));

			RooFitResult *fitResult_eta = expExtended->fitTo(*dataHistarray_eta[i], RooFit::Range("LeftTail,RightTail"), RooFit::Save());
			fitHist_pp_eta[i] = expExtended->createHistogram(Form("fitHist_pp_eta_%i", i), *x, RooFit::Binning(120, 60, 120));

			bksub_hist_pp_raw[i]->Add(fitHist_pp_raw[i], -1);
			bksub_hist_pp_eta[i]->Add(fitHist_pp_eta[i], -1);

			double Meanafter_raw = bksub_hist_pp_raw[i]->GetMean();
			double Meanafter_eta = bksub_hist_pp_eta[i]->GetMean();

			cout << "The Mean_raw before is " << Meanbefore_raw << " The Mean_raw after is " << Meanafter_raw << endl;
			cout << "The Mean_eta before is " << Meanbefore_eta << " The Mean_eta after is " << Meanafter_eta << endl;

			c1->cd();
			fitHist_pp_raw[i]->SetLineColor(kRed);
			fitHist_pp_raw[i]->Draw("HISTSAME");

			c2->cd();
			fitHist_pp_eta[i]->SetLineColor(kRed);
			fitHist_pp_eta[i]->Draw("HISTSAME");

			c1->SaveAs(Form("./expbktest/pp/raw/fitHist_pp_raw_%i.png", i));
			c2->SaveAs(Form("./expbktest/pp/eta/fitHist_pp_eta_%i.png", i));
			delete c1;
			delete c2;

			TH1D *bksub_hist_pp_raw_aftercut = new TH1D("bksub_hist_pp_raw_aftercut", "", 40, 80, 100);
			TH1D *bksub_hist_pp_eta_aftercut = new TH1D("bksub_hist_pp_eta_aftercut", "", 40, 80, 100);

			int bin_low_raw = bksub_hist_pp_raw[i]->FindBin(80);
			int bin_high_raw = bksub_hist_pp_raw[i]->FindBin(100);

			int bin_low_eta = bksub_hist_pp_eta[i]->FindBin(80);
			int bin_high_eta = bksub_hist_pp_eta[i]->FindBin(100);

			for (int j = bin_low_raw; j <= bin_high_raw; j++)
			{
				double content_pp_raw = bksub_hist_pp_raw[i]->GetBinContent(j);
				double error_pp_raw = bksub_hist_pp_raw[i]->GetBinError(j);
				double content_pp_eta = bksub_hist_pp_eta[i]->GetBinContent(j);
				double error_pp_eta = bksub_hist_pp_eta[i]->GetBinError(j);

				int new_bin_raw = bksub_hist_pp_raw_aftercut->FindBin(bksub_hist_pp_raw[i]->GetBinCenter(j));
				int new_bin_eta = bksub_hist_pp_eta_aftercut->FindBin(bksub_hist_pp_eta[i]->GetBinCenter(j));

				bksub_hist_pp_raw_aftercut->SetBinContent(new_bin_raw, content_pp_raw);
				bksub_hist_pp_raw_aftercut->SetBinError(new_bin_raw, error_pp_raw);
				bksub_hist_pp_eta_aftercut->SetBinContent(new_bin_eta, content_pp_eta);
				bksub_hist_pp_eta_aftercut->SetBinError(new_bin_eta, error_pp_eta);
			}

			c1 = new TCanvas("c1", "", 800, 600);
			c1->cd();
			bksub_hist_pp_raw_aftercut->Draw("HIST");
			c1->SaveAs(Form("./expbktest/pp/raw/histogram_after_cut_%i.png", i));
			c1->Clear();
			bksub_hist_pp_eta_aftercut->Draw("HIST");
			c1->SaveAs(Form("./expbktest/pp/eta/histogram_after_cut_%i.png", i));
			delete c1;

			dMass_raw[i] = (bksub_hist_pp_raw_aftercut->GetMean()) - 91.1876;
			dMass_Err_raw[i] = (bksub_hist_pp_raw_aftercut->GetMeanError());
			dWidth_raw[i] = (bksub_hist_pp_raw_aftercut->GetStdDev()) - 2.4955;
			dWidth_Err_raw[i] = (bksub_hist_pp_raw_aftercut->GetStdDevError());

			dMass_eta[i] = (bksub_hist_pp_eta_aftercut->GetMean()) - 91.1876;
			dMass_Err_eta[i] = (bksub_hist_pp_eta_aftercut->GetMeanError());
			dWidth_eta[i] = (bksub_hist_pp_eta_aftercut->GetStdDev()) - 2.4955;
			dWidth_Err_eta[i] = (bksub_hist_pp_eta_aftercut->GetStdDevError());

			delete bksub_hist_pp_raw_aftercut;
			delete bksub_hist_pp_eta_aftercut;
		}
	}
}
void fitHI(bool isbksub = 0)
{
	for (int i = 0; i < 6; i++)
	{
		HI_xposition[i] = i + 1;
		HI_xposition_err[i] = 0;
		bksub_hist_HI_raw[i] = (TH1D *)h_HI_mass_array_raw[i]->Clone();
		bksub_hist_HI_eta[i] = (TH1D *)h_HI_mass_array_eta[i]->Clone();

		if (!isbksub)
		{
			TH1D *bksub_hist_HI_raw_aftercut = new TH1D("bksub_hist_HI_raw_aftercut", "", 40, 80, 100);
			TH1D *bksub_hist_HI_eta_aftercut = new TH1D("bksub_hist_HI_eta_aftercut", "", 40, 80, 100);

			int bin_low_raw = bksub_hist_HI_raw[i]->FindBin(80);
			int bin_high_raw = bksub_hist_HI_raw[i]->FindBin(100);

			int bin_low_eta = bksub_hist_HI_eta[i]->FindBin(80);
			int bin_high_eta = bksub_hist_HI_eta[i]->FindBin(100);

			for (int j = bin_low_raw; j <= bin_high_raw; j++)
			{
				double content_HI_raw = bksub_hist_HI_raw[i]->GetBinContent(j);
				double error_HI_raw = bksub_hist_HI_raw[i]->GetBinError(j);
				double content_HI_eta = bksub_hist_HI_eta[i]->GetBinContent(j);
				double error_HI_eta = bksub_hist_HI_eta[i]->GetBinError(j);

				int new_bin_raw = bksub_hist_HI_raw_aftercut->FindBin(bksub_hist_HI_raw[i]->GetBinCenter(j));
				int new_bin_eta = bksub_hist_HI_eta_aftercut->FindBin(bksub_hist_HI_eta[i]->GetBinCenter(j));

				bksub_hist_HI_raw_aftercut->SetBinContent(new_bin_raw, content_HI_raw);
				bksub_hist_HI_raw_aftercut->SetBinError(new_bin_raw, error_HI_raw);
				bksub_hist_HI_eta_aftercut->SetBinContent(new_bin_eta, content_HI_eta);
				bksub_hist_HI_eta_aftercut->SetBinError(new_bin_eta, error_HI_eta);
			}

			cout << "HI_raw_before_cut is " << bksub_hist_HI_raw[i]->GetMean() << " After cut is " << bksub_hist_HI_raw_aftercut->GetMean() << endl;

			HI_dMass_raw[i] = (bksub_hist_HI_raw_aftercut->GetMean()) - 91.1876;
			HI_dMass_Err_raw[i] = (bksub_hist_HI_raw_aftercut->GetMeanError());
			HI_dWidth_raw[i] = (bksub_hist_HI_raw_aftercut->GetStdDev()) - 2.4955;
			HI_dWidth_Err_raw[i] = (bksub_hist_HI_raw_aftercut->GetStdDevError());

			HI_dMass_eta[i] = (bksub_hist_HI_eta_aftercut->GetMean()) - 91.1876;
			HI_dMass_Err_eta[i] = (bksub_hist_HI_eta_aftercut->GetMeanError());
			HI_dWidth_eta[i] = (bksub_hist_HI_eta_aftercut->GetStdDev()) - 2.4955;
			HI_dWidth_Err_eta[i] = (bksub_hist_HI_eta_aftercut->GetStdDevError());

			delete bksub_hist_HI_raw_aftercut;
			delete bksub_hist_HI_eta_aftercut;
		}
		else
		{
			TCanvas *c1 = new TCanvas("c1", "", 800, 600);
			TCanvas *c2 = new TCanvas("c2", "", 800, 600);

			HI_dataHistarray_raw[i] = new RooDataHist("HI_dataHistarray_raw", "RooDataHist from TH1D", RooArgList(*x), bksub_hist_HI_raw[i]);
			HI_dataHistarray_eta[i] = new RooDataHist("HI_dataHistarray_eta", "RooDataHist from TH1D", RooArgList(*x), bksub_hist_HI_eta[i]);
			double Meanbefore_raw = bksub_hist_HI_raw[i]->GetMean();
			double Meanbefore_eta = bksub_hist_HI_eta[i]->GetMean();

			c1->cd();
			bksub_hist_HI_raw[i]->Draw("HIST");
			c2->cd();
			bksub_hist_HI_eta[i]->Draw("HIST");

			RooFitResult *fitResult_raw = expExtended->fitTo(*HI_dataHistarray_raw[i], RooFit::Range("LeftTail,RightTail"), RooFit::Save());
			fitHist_HI_raw[i] = expExtended->createHistogram(Form("fitHist_HI_raw_%i", i), *x, RooFit::Binning(120, 60, 120));

			RooFitResult *fitResult_eta = expExtended->fitTo(*HI_dataHistarray_eta[i], RooFit::Range("LeftTail,RightTail"), RooFit::Save());
			fitHist_HI_eta[i] = expExtended->createHistogram(Form("fitHist_HI_eta_%i", i), *x, RooFit::Binning(120, 60, 120));

			bksub_hist_HI_raw[i]->Add(fitHist_HI_raw[i], -1);
			bksub_hist_HI_eta[i]->Add(fitHist_HI_eta[i], -1);

			double Meanafter_raw = bksub_hist_HI_raw[i]->GetMean();
			double Meanafter_eta = bksub_hist_HI_eta[i]->GetMean();

			cout << "The Mean_raw before is " << Meanbefore_raw << " The Mean_raw after is " << Meanafter_raw << endl;
			cout << "The Mean_eta before is " << Meanbefore_eta << " The Mean_eta after is " << Meanafter_eta << endl;

			c1->cd();
			fitHist_HI_raw[i]->SetLineColor(kRed);
			fitHist_HI_raw[i]->Draw("HISTSAME");

			c2->cd();
			fitHist_HI_eta[i]->SetLineColor(kRed);
			fitHist_HI_eta[i]->Draw("HISTSAME");

			c1->SaveAs(Form("./expbktest/HI/raw/fitHist_HI_raw_%i.png", i));
			c2->SaveAs(Form("./expbktest/HI/eta/fitHist_HI_eta_%i.png", i));
			delete c1;
			delete c2;

			TH1D *bksub_hist_HI_raw_aftercut = new TH1D("bksub_hist_HI_raw_aftercut", "", 40, 80, 100);
			TH1D *bksub_hist_HI_eta_aftercut = new TH1D("bksub_hist_HI_eta_aftercut", "", 40, 80, 100);

			int bin_low_raw = bksub_hist_HI_raw[i]->FindBin(80);
			int bin_high_raw = bksub_hist_HI_raw[i]->FindBin(100);

			int bin_low_eta = bksub_hist_HI_eta[i]->FindBin(80);
			int bin_high_eta = bksub_hist_HI_eta[i]->FindBin(100);

			for (int j = bin_low_raw; j <= bin_high_raw; j++)
			{
				double content_HI_raw = bksub_hist_HI_raw[i]->GetBinContent(j);
				double error_HI_raw = bksub_hist_HI_raw[i]->GetBinError(j);
				double content_HI_eta = bksub_hist_HI_eta[i]->GetBinContent(j);
				double error_HI_eta = bksub_hist_HI_eta[i]->GetBinError(j);

				int new_bin_raw = bksub_hist_HI_raw_aftercut->FindBin(bksub_hist_HI_raw[i]->GetBinCenter(j));
				int new_bin_eta = bksub_hist_HI_eta_aftercut->FindBin(bksub_hist_HI_eta[i]->GetBinCenter(j));

				bksub_hist_HI_raw_aftercut->SetBinContent(new_bin_raw, content_HI_raw);
				bksub_hist_HI_raw_aftercut->SetBinError(new_bin_raw, error_HI_raw);
				bksub_hist_HI_eta_aftercut->SetBinContent(new_bin_eta, content_HI_eta);
				bksub_hist_HI_eta_aftercut->SetBinError(new_bin_eta, error_HI_eta);
			}

			c1 = new TCanvas("c1", "", 800, 600);
			c1->cd();
			bksub_hist_HI_raw_aftercut->Draw("HIST");
			c1->SaveAs(Form("./expbktest/HI/raw/histogram_after_cut_%i.png", i));
			c1->Clear();
			bksub_hist_HI_eta_aftercut->Draw("HIST");
			c1->SaveAs(Form("./expbktest/HI/eta/histogram_after_cut_%i.png", i));
			delete c1;

			HI_dMass_raw[i] = (bksub_hist_HI_raw_aftercut->GetMean()) - 91.1876;
			HI_dMass_Err_raw[i] = (bksub_hist_HI_raw_aftercut->GetMeanError());
			HI_dWidth_raw[i] = (bksub_hist_HI_raw_aftercut->GetStdDev()) - 2.4955;
			HI_dWidth_Err_raw[i] = (bksub_hist_HI_raw_aftercut->GetStdDevError());

			HI_dMass_eta[i] = (bksub_hist_HI_eta_aftercut->GetMean()) - 91.1876;
			HI_dMass_Err_eta[i] = (bksub_hist_HI_eta_aftercut->GetMeanError());
			HI_dWidth_eta[i] = (bksub_hist_HI_eta_aftercut->GetStdDev()) - 2.4955;
			HI_dWidth_Err_eta[i] = (bksub_hist_HI_eta_aftercut->GetStdDevError());

			delete bksub_hist_HI_raw_aftercut;
			delete bksub_hist_HI_eta_aftercut;
		}
	}
}
void MeanfromHist(bool isbksub = 1)
{
	TFile *f1 = new TFile("new_pp_data_file_stability_readonly.root", "READ");
	TFile *f2 = new TFile("data_file.root", "READ");

	x = new RooRealVar("x", "x", 60, 120);
	x->setRange("LeftTail", 60, 65);
	x->setRange("RightTail", 110, 120);
	decayVar = new RooRealVar("decayVar", "slope", -0.1, -5, 5);
	expo = new RooExponential("exp", "Exponential PDF", *x, *decayVar);
	nEvents = new RooRealVar("nEvents", "number of events", 1000, 0, 100000);
	expExtended = new RooExtendPdf("expExtended", "Extended Exponential PDF", *expo, *nEvents);

	GetHistHI(f2);
	GetHistpp(f1);

	fitpp(isbksub);
	fitHI(isbksub);

	TGraphErrors *g_pp_dM_raw = new TGraphErrors(22, xposition, dMass_raw, xposition_err, dMass_Err_raw);
	TGraphErrors *g_pp_dM_eta = new TGraphErrors(22, xposition, dMass_eta, xposition_err, dMass_Err_eta);
	TGraphErrors *g_pp_dWidth_raw = new TGraphErrors(22, xposition, dWidth_raw, xposition_err, dWidth_Err_raw);
	TGraphErrors *g_pp_dWidth_eta = new TGraphErrors(22, xposition, dWidth_eta, xposition_err, dWidth_Err_eta);

	TGraphErrors *g_HI_dM_raw = new TGraphErrors(6, HI_xposition, HI_dMass_raw, HI_xposition_err, HI_dMass_Err_raw);
	TGraphErrors *g_HI_dM_eta = new TGraphErrors(6, HI_xposition, HI_dMass_eta, HI_xposition_err, HI_dMass_Err_eta);
	TGraphErrors *g_HI_dWidth_raw = new TGraphErrors(6, HI_xposition, HI_dWidth_raw, HI_xposition_err, HI_dWidth_Err_raw);
	TGraphErrors *g_HI_dWidth_eta = new TGraphErrors(6, HI_xposition, HI_dWidth_eta, HI_xposition_err, HI_dWidth_Err_eta);

	TFile *f3 = new TFile("All_plots.root", "UPDATE");
	f3->cd();
	if (isbksub)
	{
		g_pp_dM_raw->Write("hist_pp_dM_raw_bksub", 2);
		g_pp_dM_eta->Write("hist_pp_dM_eta_bksub", 2);
		g_pp_dWidth_raw->Write("hist_pp_dWidth_raw_bksub", 2);
		g_pp_dWidth_eta->Write("hist_pp_dWidth_eta_bksub", 2);

		g_HI_dM_raw->Write("hist_HI_dM_raw_bksub", 2);
		g_HI_dM_eta->Write("hist_HI_dM_eta_bksub", 2);
		g_HI_dWidth_raw->Write("hist_HI_dWidth_raw_bksub", 2);
		g_HI_dWidth_eta->Write("hist_HI_dWidth_eta_bksub", 2);
	}
	else
	{
		g_pp_dM_raw->Write("hist_pp_dM_raw_nobksub", 2);
		g_pp_dM_eta->Write("hist_pp_dM_eta_nobksub", 2);
		g_pp_dWidth_raw->Write("hist_pp_dWidth_raw_nobksub", 2);
		g_pp_dWidth_eta->Write("hist_pp_dWidth_eta_nobksub", 2);

		g_HI_dM_raw->Write("hist_HI_dM_raw_nobksub", 2);
		g_HI_dM_eta->Write("hist_HI_dM_eta_nobksub", 2);
		g_HI_dWidth_raw->Write("hist_HI_dWidth_raw_nobksub", 2);
		g_HI_dWidth_eta->Write("hist_HI_dWidth_eta_nobksub", 2);
	}
}