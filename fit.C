#include "fit.h"

void setupfitvariable()
{
	x = (RooRealVar *)roodataarray[0]->get()->find("x");
	x->setBinning(RooBinning(5000, 60, 120), "cache");
	bwmean = new RooRealVar("bwmean", "bwmean", 90.8513, 80, 100);
	width = new RooRealVar("width", "width", 2.911, 0, 10);
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width);
	cbmean = new RooRealVar("mean", "mean", 0);
	cbsigma = new RooRealVar("sigma", "sigma", 1.1104, 0, 10);
	cbalpha = new RooRealVar("alpha", "alpha", 1.8, -5, 5);
	cbn = new RooRealVar("n", "n", 0.972, -5, 5);
	cb = new RooCBShape("cb", "cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);

	decayVar = new RooRealVar("decayVar", "decayVar", -0.11, -5, +5);
	expo = new RooExponential("expo", "Exponential PDF", *x, *decayVar);

	newconvpdf = new RooFFTConvPdf("newconvpdf", "newconvpdf", *x, *bw, *cb);

	fsig = new RooRealVar("fsig", "signal fraction", 0.99, 0, 1);
	purepdf = new RooAddPdf("purepdf", "cbconbw+exp", RooArgList(*newconvpdf, *expo), *fsig);
}

void dofit(RooDataSet *dataset, int iteration, const char *type)
{
	// if iteration == -99 means fit the inclusive and without dM dWidth TGraphError filling.

	RooFitResult *fitresult;
	RooFitResult *fitresult1;

	if (fittype == 1)
		fitresult = newconvpdf->fitTo(*dataset, RooFit::Save(true), RooFit::NumCPU(8), RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migrad"));
	if (fittype == 2)
		fitresult = purepdf->fitTo(*dataset, RooFit::Save(true), RooFit::NumCPU(8), RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migrad"));

	frame = x->frame(RooFit::Title("Z mass fit"));
	framecheck = x->frame(RooFit::Title("Background"));

	dataset->plotOn(frame, RooFit::Name("datahistogram"), RooFit::Binning(60), RooFit::MarkerColor(kBlack), RooFit::MarkerSize(0.1));
	if (fittype == 1)
	{
		newconvpdf->plotOn(frame, RooFit::Name("fitcurve"), RooFit::LineWidth(1));
		newconvpdf->paramOn(frame, RooFit::Format("NEU", RooFit::AutoPrecision(3)), RooFit::Layout(0.6, 1, 0.9), RooFit::ShowConstants(kTRUE));
	}
	if (fittype == 2)
	{
		purepdf->plotOn(frame, RooFit::Name("fitcurve"), RooFit::LineWidth(1));
		purepdf->paramOn(frame, RooFit::Format("NEU", RooFit::AutoPrecision(3)), RooFit::Layout(0.6, 1, 0.9), RooFit::ShowConstants(kTRUE));
		purepdf->plotOn(framecheck, RooFit::Components(*expo), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
	}
	residuals = frame->residHist("datahistogram", "fitcurve", true, false); // true = pull, false = center of the bin

	// Calculate Chi2/ndf
	int nParams = fitresult->floatParsFinal().getSize();
	double chi2ndf = frame->chiSquare("datahistogram", "fitcurve", nParams);

	cout << "The chi2/ndf is " << chi2ndf << endl;

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
	Int_t numEntries = roodataarray[i]->numEntries();

	textBox = new TPaveText(0.1, 0.5, 0.3, 0.8, "NDC");
	textBox->SetFillColor(0);
	textBox->SetTextSize(0.03);
	textBox->AddText(Form("Mean value: %.2f", meanofBW));
	textBox->AddText(Form("Width value: %.2f", widthofBW));
	textBox->AddText(Form("# of Entries: %i", numEntries));
	textBox->AddText(Form("chi2/ndf: %i", chi2ndf));
	textBox->AddText("Raw,Asym");
	textBox->Draw();
	c1->cd(2);
	pullFrame->Draw();
	c1->cd(3);
	framecheck->Draw();

	if (iteration != -99)
	{
		if (type == "raw")
		{
			dMass_raw[iteration] = meanofBW;
			dMass_Err_raw[iteration] = meanerrorofBW;
			dWidth_raw[iteration] = widthofBW;
			dWidth_Err_raw[iteration] = widtherrorofBW;
		}
		if (type == "eta")
		{
			dMass_eta[iteration] = meanofBW;
			dMass_Err_eta[iteration] = meanerrorofBW;
			dWidth_eta[iteration] = widthofBW;
			dWidth_Err_eta[iteration] = widtherrorofBW;
		}

		if (fittype == 1)
		{
			if (type == "raw")
			{
				c1->SaveAs(Form("./fit/raw/run_period_%i.png", iteration));
			}
			else
			{
				c1->SaveAs(Form("./fit/eta/run_period_%i.png", iteration));
			}
		}

		if (fittype == 2)
			if (type == "raw")
			{
				c1->SaveAs(Form("./fit/raw/run_period_exp_%i.png", iteration));
			}
			else
			{
				c1->SaveAs(Form("./fit/eta/run_period_exp_%i.png", iteration));
			}
	}
}

void fit(int fittype = 2)
{

	TFile *f1 = new TFile("pp_data_file_stability.root", "READ");

	RooDataSet *mass_array_eta[22];
	RooDataSet *mass_array_raw[22];
	TH1D *h_mass_array_raw[22];
	TH1D *h_mass_array_eta[22];

	TH1D *h_mass_array_raw_inclusive;
	RooDataSet *mass_array_raw_inclusive;
	TH1D *h_mass_array_eta_inclusive;
	RooDataSet *mass_array_eta_inclusive;

	Double_t dMass_raw[22];
	Double_t dMass_Err_raw[22];
	Double_t dWidth_raw[22];
	Double_t dWidth_Err_raw[22];

	Double_t dMass_eta[22];
	Double_t dMass_Err_eta[22];
	Double_t dWidth_eta[22];
	Double_t dWidth_Err_eta[22];

	Double_t xposition[22];
	Double_t xposition_err[22];

	h_mass_array_raw_inclusive = (TH1D *)f1->Get("h_mass_array_raw_inclusive");
	h_mass_array_eta_inclusive = (TH1D *)f1->Get("h_mass_array_eta_inclusive");
	mass_array_raw_inclusive = (RooDataSet *)f1->Get("mass_array_raw_inclusive");
	mass_array_eta_inclusive = (RooDataSet *)f1->Get("mass_array_eta_inclusive");

	for (int i = 0; i < 22; i++)
	{
		mass_array_raw[i] = (RooDataSet *)f1->Get(Form("mass_array_raw_%i", i));
		mass_array_eta[i] = (RooDataSet *)f1->Get(Form("mass_array_eta_%i", i));
		h_mass_array_raw[i] = (TH1D *)f1->Get(Form("h_mass_array_raw_%i", i));
		h_mass_array_eta[i] = (TH1D *)f1->Get(Form("h_mass_array_eta_%i", i));
		xposition[i] = i + 1;
		xposition_err[i] = 0;
	}

	setupfitvariable();

	dofit(mass_array_raw_inclusive, -99, "raw");
	dofit(mass_array_eta_inclusive, -99, "eta");

	for (int i = 0; i < 22; i++)
	{
		dofit(mass_array_raw[i], i, "raw");
		dofit(mass_array_eta[i], i, "eta");
	}

	TGraphErrors *g_1 = new TGraphErrors(22, xposition, dMass_raw, xposition_err_raw, dMass_Err_raw);
	TGraphErrors *g_2 = new TGraphErrors(22, xposition, dMass_eta, xposition_err_eta, dMass_Err_eta);
	TGraphErrors *g_3 = new TGraphErrors(22, xposition, dWidth_raw, xposition_err_raw, dWidth_Err_raw);
	TGraphErrors *g_4 = new TGraphErrors(22, xposition, dWidth_eta, xposition_err_eta, dWidth_Err_eta);

	TFile *f2 = new TFile("./All_plots.root", "UPDATE");

	f2->cd();
	if (fittype == 1)
	{
		g_1->Write("pp_raw_dM", 2);
		g_2->Write("pp_eta_dM", 2);
		g_3->Write("pp_raw_dW", 2);
		g_4->Write("pp_eta_dW", 2);
	}
	if (fittype == 2)
	{
		g_1->Write("pp_raw_dM_exp", 2);
		g_2->Write("pp_eta_dM_exp", 2);
		g_3->Write("pp_raw_dW_exp", 2);
		g_4->Write("pp_eta_dW_exp", 2);
	}
	f2->Close();
	f1->Close();
}