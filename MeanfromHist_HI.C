void MeanfromHist_HI(){
	TFile *f1 = new TFile("data_file.root","UPDATE");
	TH1D* mass_array_data[6];

	Double_t dMass[6];
	Double_t dMass_Err[6];

	Double_t dWidth[6];
	Double_t dWidth_Err[6];

	Double_t xposition[6];
	Double_t xposition_err[6];

	Double_t dMass_before[6];
	Double_t dMass_Err_before[6];

	Double_t dWidth_before[6];
	Double_t dWidth_Err_before[6];

	RooDataHist *dataHistarray[6];
	RooRealVar x("x", "x", 60, 120);
	RooRealVar tau("tau", "slope", -0.1, -5, +5);
	RooExponential exp("exp", "Exponential PDF", x, tau);

	RooRealVar nEvents("nEvents", "number of events", 1000, 0, 100000);
	RooExtendPdf expExtended("expExtended", "Extended Exponential PDF", exp, nEvents);

	x.setRange("LeftTail",60,65);
	x.setRange("RightTail",110,120);

	TH1D* bksub_hist[6];
	TH1* fitHist[6];

	TCanvas *c1 = new TCanvas("c1","c1",800,600);

	for (int i=1; i<=5; i++){
		mass_array_data[i] = (TH1D*)f1->Get(Form("mass_array_data_witheta_%i",i-1));
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
		c1->SaveAs(Form("./expbktest/HI_%i.png",i));

		dMass[i] = (bksub_hist[i]->GetMean()) - 91.1876;
		dMass_Err[i] = (bksub_hist[i]->GetMeanError());

		dWidth[i] = (bksub_hist[i]->GetStdDev()) - 5.4955;
		cout << "The Width is " << dWidth[i] << endl;
		dWidth_Err[i] = (bksub_hist[i]->GetStdDevError());
	}

	xposition[0] = 1;
	xposition_err[0] = 0;
	mass_array_data[0] = (TH1D*)f1->Get("mass_array_data_witheta_10");

	bksub_hist[0] = (TH1D*) mass_array_data[0]->Clone();
	dMass_before[0] = (bksub_hist[0]->GetMean()) - 91.1876;
	dMass_Err_before[0] = (bksub_hist[0]->GetMeanError());

	dWidth_before[0] = (bksub_hist[0]->GetStdDev()) - 5.4955;
	dWidth_Err_before[0] = (bksub_hist[0]->GetStdDevError());

	cout << "STD is " << bksub_hist[0]->GetStdDev() << endl;

	dataHistarray[0] = new RooDataHist("dataHist", "RooDataHist from TH1D", RooArgList(x), bksub_hist[0]);

	Double_t meanbefore = bksub_hist[0]->GetMean();
	bksub_hist[0]->Draw("HIST");
	RooFitResult* fitResult = expExtended.fitTo(*dataHistarray[0], RooFit::Range("LeftTail,RightTail"),RooFit::Save());
	//RooFitResult* fitResult = exp.fitTo(*dataHistarray[i],RooFit::Save());
	//fitResult->Print("v");
	fitHist[0] = expExtended.createHistogram(Form("fitHist_%i",0), x, RooFit::Binning(120, 60, 120));
	bksub_hist[0]->Add(fitHist[0],-1);
	Double_t meanafter = bksub_hist[0]->GetMean();

	cout << "The Mean Before is " << meanbefore << " The Mean after is " << meanafter << endl;
	fitHist[0]->SetLineColor(kRed);
	fitHist[0]->Draw("HISTSAME");
	c1->SaveAs(Form("./expbktest/HI_%i.png",0));

	dMass[0] = (bksub_hist[0]->GetMean()) - 91.1876;
	dMass_Err[0] = (bksub_hist[0]->GetMeanError());

	dWidth[0] = (bksub_hist[0]->GetStdDev()) - 5.4955;
	dWidth_Err[0] = (bksub_hist[0]->GetStdDevError());


	TGraphErrors *g_1 = new TGraphErrors(6,xposition,dMass,xposition_err,dMass_Err);
	TGraphErrors *g_2 = new TGraphErrors(6,xposition,dMass_before,xposition_err,dMass_Err_before);

	TGraphErrors *g_w_1 = new TGraphErrors(6,xposition,dWidth,xposition_err,dWidth_Err);
	TGraphErrors *g_w_2 = new TGraphErrors(6,xposition,dWidth_before,xposition_err,dWidth_Err_before);

	TFile *f2 = new TFile("pp_data_file_stability.root","UPDATE");
	f2->cd();
	g_1->Write("HI_dM_bksub_hist",2);
	g_2->Write("HI_dM_hist",2);
	g_w_1->Write("HI_dWidth_bksub_hist",2);
	g_w_2->Write("HI_dWidth_hist",2);

	f2->Close();

}