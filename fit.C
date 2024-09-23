#include "fit.h"

void fit(int fittype = 2){

	TFile *f1 = new TFile("pp_data_file_stability.root","UPDATE");

	RooDataSet* roodataarray[22];
	TH1D* mass_array_data[22];

	Double_t dMass[22];
	Double_t dMass_Err[22];
	Double_t dWidth[22];
	Double_t dWidth_Err[22];
	Double_t xposition[22];
	Double_t xposition_err[22];

	for (int i=0; i<22; i++){
		roodataarray[i] = (RooDataSet*)f1->Get(Form("roodata_%i",i));
		mass_array_data[i] = (TH1D*)f1->Get(Form("mass_array_data_%i",i));
		xposition[i] = i+1; 
		xposition_err[i] = 0;
	}

	x = (RooRealVar*)roodataarray[0]->get()->find("x");
	x->setBinning(RooBinning(5000,60,120),"cache");
	bwmean = new RooRealVar("bwmean", "bwmean", 90.8513 ,80,100);
	width = new RooRealVar("width", "width",2.911,0,10 );
	bw = new RooBreitWigner("bw", "bw", *x, *bwmean, *width );
	cbmean = new RooRealVar("mean","mean",0);
	cbsigma = new RooRealVar("sigma","sigma",1.1104,0,10);
	cbalpha = new RooRealVar("alpha","alpha",1.8,-5,5);
	cbn = new RooRealVar("n","n", 0.972 ,-5,5);
	cb = new RooCBShape("cb","cb", *x, *cbmean, *cbsigma, *cbalpha, *cbn);

	decayVar = new RooRealVar("decayVar","decayVar",-0.11,-5,+5);
	expo = new RooExponential("expo", "Exponential PDF", *x, *decayVar);

	newconvpdf = new RooFFTConvPdf("newconvpdf","newconvpdf", *x, *bw, *cb);

	fsig = new RooRealVar("fsig","signal fraction",0.99,0,1);
    purepdf = new RooAddPdf("purepdf","cbconbw+exp",RooArgList(*newconvpdf,*expo),*fsig);



	for (int i=0; i<22; i++){
		RooFitResult* fitresult;
		RooFitResult* fitresult1;

		if (fittype == 1) fitresult = newconvpdf->fitTo(*roodataarray[i],RooFit::Save(true),RooFit::NumCPU(8),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
		if (fittype == 2) fitresult = purepdf->fitTo(*roodataarray[i],RooFit::Save(true),RooFit::NumCPU(8),RooFit::SumW2Error(true),RooFit::Minimizer("Minuit2","migrad"));
	
		frame = x->frame(RooFit::Title("Z mass fit"));
		framecheck = x->frame(RooFit::Title("Background"));

		roodataarray[i]->plotOn(frame,RooFit::Name("roodata"),RooFit::Binning(60),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(0.1));
		if (fittype == 1){
			newconvpdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			newconvpdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
		}
		if (fittype == 2){
			purepdf->plotOn(frame,RooFit::Name("fit"),RooFit::LineWidth(1));
			purepdf->paramOn(frame,RooFit::Format("NEU",RooFit::AutoPrecision(3)),RooFit::Layout(0.6,1,0.9),RooFit::ShowConstants(kTRUE));
			purepdf->plotOn(framecheck,RooFit::Components(*expo), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
		}
		residuals = frame->residHist("roodata","fit",true,false); //true = pull, false = center of the bin

		pullFrame = x->frame();
		pullFrame->addPlotable(residuals, "P");
		pullFrame->SetTitle("");
		pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
		frame->SetTitle("");
		frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");	

		c1 = new TCanvas("c1","",1600,1200);
		c1->Divide(2,2);
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
		textBox->AddText(Form("Width value: %.2f",widthofBW));
		textBox->AddText(Form("# of Entries: %i",numEntries));
		textBox->AddText("Raw,Asym");
		textBox->Draw();
		c1->cd(2);
		pullFrame->Draw();
		c1->cd(3);
		framecheck->Draw();

		dMass[i] = meanofBW;
		cout << "Mean of BW is " << meanofBW << endl;
		dMass_Err[i] = meanerrorofBW;
		dWidth[i] = widthofBW;
		dWidth_Err[i] = widtherrorofBW;

		if (fittype == 1)c1->SaveAs(Form("./fit/run_period_%i.png",i));
		if (fittype == 2)c1->SaveAs(Form("./fit/run_period_expbk_%i.png",i));
	}

	TGraphErrors *g_1 = new TGraphErrors(22,xposition,dMass,xposition_err,dMass_Err);
	TGraphErrors *g_2 = new TGraphErrors(22,xposition,dWidth,xposition_err,dWidth_Err);

	f1->cd();
		if (fittype == 1){
		g_1->Write("runperiod_22_dM",2);
		g_2->Write("runperiod_22_dWidth",2);
	}
	if (fittype == 2){
		g_1->Write("runperiod_22_expbk_dM",2);
		g_2->Write("runperiod_22_expbk_dWidth",2);
	}
	f1->Close();





}