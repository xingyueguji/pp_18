void MeanfromHist(){
	TFile *f1 = new TFile("pp_data_file_stability.root","UPDATE");
	TH1D* mass_array_data[22];

	Double_t dMass[22];
	Double_t dMass_Err[22];

	Double_t dWidth[22];
	Double_t dWidth_Err[22];

	Double_t dMass_before[22];
	Double_t dMass_Err_before[22];

	Double_t dWidth_before[22];
	Double_t dWidth_Err_before[22];

	Double_t xposition[22];
	Double_t xposition_err[22];

	RooDataHist *dataHistarray[22];
	RooRealVar x("x", "x", 60, 120);
	RooRealVar tau("tau", "slope", -0.1, -5, +5);
	RooExponential exp("exp", "Exponential PDF", x, tau);

	RooRealVar nEvents("nEvents", "number of events", 1000, 0, 100000);
	RooExtendPdf expExtended("expExtended", "Extended Exponential PDF", exp, nEvents);

	x.setRange("LeftTail",60,65);
	x.setRange("RightTail",110,120);

	TH1D* bksub_hist[22];
	TH1* fitHist[22];

	TCanvas *c1 = new TCanvas("c1","c1",800,600);

	for (int i=0; i<22; i++){
		mass_array_data[i] = (TH1D*)f1->Get(Form("mass_array_data_%i",i));
		xposition[i] = i+1; 
		xposition_err[i] = 0;

		bksub_hist[i] = (TH1D*) mass_array_data[i]->Clone();
		dMass_before[i] = (bksub_hist[i]->GetMean()) - 91.1876;
		dMass_Err_before[i] = (bksub_hist[i]->GetMeanError());

		dWidth_before[i] = (bksub_hist[i]->GetStdDev()) - 5.4955;
		dWidth_Err_before[i] = (bksub_hist[i]->GetStdDevError());

		cout << "STD is " << bksub_hist[i]->GetStdDev() << endl;

		dataHistarray[i] = new RooDataHist("dataHist", "RooDataHist from TH1D", RooArgList(x), bksub_hist[i]);

		Double_t meanbefore = bksub_hist[i]->GetMean();
		bksub_hist[i]->Draw("HIST");
		RooFitResult* fitResult = expExtended.fitTo(*dataHistarray[i], RooFit::Range("LeftTail,RightTail"),RooFit::Save());
		//RooFitResult* fitResult = exp.fitTo(*dataHistarray[i],RooFit::Save());
		//fitResult->Print("v");
		fitHist[i] = expExtended.createHistogram(Form("fitHist_%i",i), x, RooFit::Binning(120, 60, 120));

		bksub_hist[i]->Add(fitHist[i],-1);
		Double_t meanafter = bksub_hist[i]->GetMean();

		cout << "The Mean Before is " << meanbefore << " The Mean after is " << meanafter << endl;
		fitHist[i]->SetLineColor(kRed);
		fitHist[i]->Draw("HISTSAME");
		c1->SaveAs(Form("./expbktest/pp_%i.png",i));

		dMass[i] = (bksub_hist[i]->GetMean()) - 91.1876;
		dMass_Err[i] = (bksub_hist[i]->GetMeanError());

		dWidth[i] = (bksub_hist[i]->GetStdDev()) - 5.4955;
		dWidth_Err[i] = (bksub_hist[i]->GetStdDevError());

	}




	TGraphErrors *g_1 = new TGraphErrors(22,xposition,dMass,xposition_err,dMass_Err);
	TGraphErrors *g_2 = new TGraphErrors(22,xposition,dMass_before,xposition_err,dMass_Err_before);

	TGraphErrors *g_w_1 = new TGraphErrors(22,xposition,dWidth,xposition_err,dWidth_Err);
	TGraphErrors *g_w_2 = new TGraphErrors(22,xposition,dWidth_before,xposition_err,dWidth_Err_before);

	g_1->Write("run_period_dM_hist_bksub",2);
	g_2->Write("run_period_dM_hist",2);

	g_w_1->Write("run_period_dWidth_hist_bksub",2);
	g_w_2->Write("run_period_dWidth_hist",2);

}