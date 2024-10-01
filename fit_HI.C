#include "fit.h"

void setupfitvariable()
{
	x = (RooRealVar *)HI_mass_array_raw[0]->get()->find("roomass");
	x->setBinning(RooBinning(10000, 60, 120), "cache");
	bwmean = new RooRealVar("bwmean", "bwmean", 90.8513, 80, 100);
	width = new RooRealVar("width", "width", 2.911, 0, 10);
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width);
	cbmean = new RooRealVar("mean", "mean", 0);
	decayVar = new RooRealVar("decayVar", "decayVar", -0.11, -5, +5);
	expo = new RooExponential("expo", "Exponential PDF", *x, *decayVar);
}
void setupfitvariableHI(bool isfix, bool isbksub, bool isetacut)
{
	if (isfix)
	{
		if (isbksub)
		{
			if (isetacut)
			{
				cbalpha = new RooRealVar("alpha", "alpha", 1.8012); // Combined_alpha_exp_eta
				cbn = new RooRealVar("n", "n", 0.9523);				// Combined_n_exp_eta
																	// cbsigma = new RooRealVar("sigma", "sigma", 0.8281); // Combined_STD_exp_eta
			}
			else
			{
				cbalpha = new RooRealVar("alpha", "alpha", 1.8202); // Combined_alpha_exp_raw
				cbn = new RooRealVar("n", "n", 1.0352);				// Combined_n_exp_raw
																	// cbsigma = new RooRealVar("sigma", "sigma", 1.0787); // Combined_STD_exp_raw
			}
		}

		else
		{
			if (isetacut)
			{
				cbalpha = new RooRealVar("alpha", "alpha", 1.9404); // Combined_alpha_eta
				cbn = new RooRealVar("n", "n", 0.8015);				// Combined_n_eta
																	// cbsigma = new RooRealVar("sigma", "sigma", 0.8010); // Combined_STD_eta
			}
			else
			{
				cbalpha = new RooRealVar("alpha", "alpha", 1.8927); // Combined_alpha_raw
				cbn = new RooRealVar("n", "n", 0.9353);				// Combined_n_raw
																	// cbsigma = new RooRealVar("sigma", "sigma", 1.0805); // Combined_STD_raw
			}
		}
	}
	else
	{
		cbalpha = new RooRealVar("alpha", "alpha", 1.8, -5, 5);
		cbn = new RooRealVar("n", "n", 0.972, -5, 5);
	}
	cbsigma = new RooRealVar("sigma", "sigma", 1.1104, 0, 10);
	cb = new RooCBShape("cb", "cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);
	newconvpdf = new RooFFTConvPdf("newconvpdf", "newconvpdf", *x, *bw, *cb);
	fsig = new RooRealVar("fsig", "signal fraction", 0.99, 0, 1);
	purepdf = new RooAddPdf("purepdf", "cbconbw+exp", RooArgList(*newconvpdf, *expo), *fsig);
}

void dofit(RooDataSet *dataset, int iteration, string type)
{
	RooFitResult *fitresult;
	RooFitResult *fitresult1;
	// RooFit::SumW2Error(true)

	// RooFit::AsymptoticError(false)
	if (fittype == 1)
		fitresult = newconvpdf->fitTo(*dataset, RooFit::Save(true), RooFit::NumCPU(8), RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migrad"));
	if (fittype == 2)
		fitresult = purepdf->fitTo(*dataset, RooFit::Save(true), RooFit::NumCPU(8), RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migrad"));

	frame = x->frame(RooFit::Title("Z mass fit"));
	framecheck = x->frame(RooFit::Title("Background"));

	dataset->plotOn(frame, RooFit::Name("roodata"), RooFit::Binning(60), RooFit::MarkerColor(kBlack), RooFit::MarkerSize(0.1));
	if (fittype == 1)
	{
		newconvpdf->plotOn(frame, RooFit::Name("fit"), RooFit::LineWidth(1));
		newconvpdf->paramOn(frame, RooFit::Format("NEU", RooFit::AutoPrecision(3)), RooFit::Layout(0.6, 1, 0.9), RooFit::ShowConstants(kTRUE));
	}
	if (fittype == 2)
	{
		purepdf->plotOn(frame, RooFit::Name("fit"), RooFit::LineWidth(1));
		purepdf->paramOn(frame, RooFit::Format("NEU", RooFit::AutoPrecision(3)), RooFit::Layout(0.6, 1, 0.9), RooFit::ShowConstants(kTRUE));
		purepdf->plotOn(framecheck, RooFit::Components(*expo), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
	}
	residuals = frame->residHist("roodata", "fit", true, false); // true = pull, false = center of the bin
	int nParams = fitresult->floatParsFinal().getSize();
	double chi2ndf = frame->chiSquare(nParams);

	pullFrame = x->frame();
	pullFrame->addPlotable(residuals, "P");
	pullFrame->SetTitle("");
	pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
	frame->SetTitle("");
	frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");

	c1 = new TCanvas("c1", "", 1600, 1200);
	c1->Divide(2, 2);
	c1->cd(1);
	frame->Draw();

	double meanofBW = bwmean->getVal() - 91.1876;
	double meanerrorofBW = bwmean->getError();
	double widthofBW = width->getVal() - 2.4955;
	double widtherrorofBW = width->getError();
	Int_t numEntries = dataset->numEntries();

	Double_t alpha = cbalpha->getVal();
	Double_t alpha_err = cbalpha->getError();
	Double_t n = cbn->getVal();
	Double_t n_err = cbn->getError();
	Double_t sigma = cbsigma->getVal();
	Double_t sigma_err = cbsigma->getError();

	textBox = new TPaveText(0.1, 0.5, 0.3, 0.8, "NDC");
	textBox->SetFillColor(0);
	textBox->SetTextSize(0.03);
	textBox->AddText(Form("Mean value: %.2f", meanofBW));
	textBox->AddText(Form("Width value: %.2f", widthofBW));
	textBox->AddText(Form("# of Entries: %i", numEntries));
	textBox->AddText(Form("chi2/ndf: %f", chi2ndf));
	textBox->Draw();
	c1->cd(2);
	pullFrame->Draw();
	c1->cd(3);
	framecheck->Draw();

	if (type == "raw")
	{
		HI_dMass_raw[iteration] = meanofBW;
		HI_dMass_Err_raw[iteration] = meanerrorofBW;
		HI_dWidth_raw[iteration] = widthofBW;
		HI_dWidth_Err_raw[iteration] = widtherrorofBW;
		HI_Alpha_raw[iteration] = alpha;
		HI_Alpha_Err_raw[iteration] = alpha_err;
		HI_N_raw[iteration] = n;
		HI_N_Err_raw[iteration] = n_err;
		HI_Std_raw[iteration] = sigma;
		HI_Std_Err_raw[iteration] = sigma_err;
	}
	if (type == "eta")
	{
		HI_dMass_eta[iteration] = meanofBW;
		HI_dMass_Err_eta[iteration] = meanerrorofBW;
		HI_dWidth_eta[iteration] = widthofBW;
		HI_dWidth_Err_eta[iteration] = widtherrorofBW;
		HI_Alpha_eta[iteration] = alpha;
		HI_Alpha_Err_eta[iteration] = alpha_err;
		HI_N_eta[iteration] = n;
		HI_N_Err_eta[iteration] = n_err;
		HI_Std_eta[iteration] = sigma;
		HI_Std_Err_eta[iteration] = sigma_err;
	}

	if (fittype == 1)
	{
		if (type == "raw")
		{
			c1->SaveAs(Form("./fitHI/raw/HI_18_%i.png", iteration));
		}
		else
		{
			c1->SaveAs(Form("./fitHI/eta/HI_18_%i.png", iteration));
		}
	}

	if (fittype == 2)
	{
		if (type == "raw")
		{
			c1->SaveAs(Form("./fitHI/raw/HI_18_exp_%i.png", iteration));
		}
		else
		{
			c1->SaveAs(Form("./fitHI/eta/HI_18_exp_%i.png", iteration));
		}
	}
}

void fit_HI(int nobkorexp = 2, bool isfix = 1)
{
	fittype = nobkorexp;
	TFile *f1 = new TFile("data_file.root", "READ");

	HI_xposition[0] = 1;
	HI_xposition_err[0] = 0;
	HI_mass_array_raw[0] = (RooDataSet *)f1->Get("roodata_10");
	HI_mass_array_eta[0] = (RooDataSet *)f1->Get("roodata_eta_10");

	for (int i = 1; i <= 5; i++)
	{
		// Reserve 0th position for inclusive
		HI_mass_array_raw[i] = (RooDataSet *)f1->Get(Form("roodata_%i", i - 1));
		HI_mass_array_eta[i] = (RooDataSet *)f1->Get(Form("roodata_eta_%i", i - 1));
		HI_xposition[i] = i + 1;
		HI_xposition_err[i] = 0;
	}

	bool isbksub;

	if (nobkorexp == 1)
		isbksub = false;
	if (nobkorexp == 2)
		isbksub = true;

	setupfitvariable();

	if (!isfix)
	{
		setupfitvariableHI(isfix, isbksub, 0); // Here the iseta does not matter
	}

	for (int i = 0; i < 6; i++)
	{
		if (isfix)
		{
			setupfitvariableHI(isfix, isbksub, 0);
		}
		dofit(HI_mass_array_raw[i], i, "raw");
		if (isfix)
		{
			setupfitvariableHI(isfix, isbksub, 1);
		}
		dofit(HI_mass_array_eta[i], i, "eta");
	}

	cout << "HI_xposition_is" << HI_xposition[0] << endl;
	cout << "Mass is " << HI_dMass_raw[0] << endl;

	TGraphErrors *g_1 = new TGraphErrors(6, HI_xposition, HI_dMass_raw, HI_xposition_err, HI_dMass_Err_raw);
	TGraphErrors *g_2 = new TGraphErrors(6, HI_xposition, HI_dMass_eta, HI_xposition_err, HI_dMass_Err_eta);

	TGraphErrors *g_3 = new TGraphErrors(6, HI_xposition, HI_dWidth_raw, HI_xposition_err, HI_dWidth_Err_raw);
	TGraphErrors *g_4 = new TGraphErrors(6, HI_xposition, HI_dWidth_eta, HI_xposition_err, HI_dWidth_Err_eta);

	TGraphErrors *g_5 = new TGraphErrors(6, HI_xposition, HI_Alpha_raw, HI_xposition_err, HI_Alpha_Err_raw);
	TGraphErrors *g_6 = new TGraphErrors(6, HI_xposition, HI_Alpha_eta, HI_xposition_err, HI_Alpha_Err_raw);

	TGraphErrors *g_7 = new TGraphErrors(6, HI_xposition, HI_N_raw, HI_xposition_err, HI_N_Err_raw);
	TGraphErrors *g_8 = new TGraphErrors(6, HI_xposition, HI_N_eta, HI_xposition_err, HI_N_Err_eta);

	TGraphErrors *g_9 = new TGraphErrors(6, HI_xposition, HI_Std_raw, HI_xposition_err, HI_Std_Err_raw);
	TGraphErrors *g_10 = new TGraphErrors(6, HI_xposition, HI_Std_eta, HI_xposition_err, HI_Std_Err_eta);

	TFile *f2 = new TFile("All_plots.root", "UPDATE");
	f2->cd();
	if (fittype == 1)
	{
		g_1->Write("HI_raw_dM", 2);
		g_2->Write("HI_eta_dM", 2);

		g_3->Write("HI_raw_dW", 2);
		g_4->Write("HI_eta_dW", 2);

		g_5->Write("HI_raw_alpha", 2);
		g_6->Write("HI_eta_alpha", 2);

		g_7->Write("HI_raw_N", 2);
		g_8->Write("HI_eta_N", 2);

		g_9->Write("HI_raw_sigma", 2);
		g_10->Write("HI_eta_sigma", 2);
	}
	if (fittype == 2)
	{
		g_1->Write("HI_raw_dM_exp", 2);
		g_2->Write("HI_eta_dM_exp", 2);

		g_3->Write("HI_raw_dW_exp", 2);
		g_4->Write("HI_eta_dW_exp", 2);

		g_5->Write("HI_raw_alpha_exp", 2);
		g_6->Write("HI_eta_alpha_exp", 2);

		g_7->Write("HI_raw_N_exp", 2);
		g_8->Write("HI_eta_N_exp", 2);

		g_9->Write("HI_raw_sigma_exp", 2);
		g_10->Write("HI_eta_sigma_exp", 2);
	}
	f1->Close();
	f2->Close();
}