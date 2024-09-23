#include "fit.h"
void fit_HI(int fittype = 2){

	TFile *f1 = new TFile("data_file.root","READ");

	RooDataSet* roodataarray[6];

	Double_t dMass[6];
	Double_t dMass_Err[6];
	Double_t dWidth[6];
	Double_t dWidth_Err[6];
	Double_t xposition[6];
	Double_t xposition_err[6];

	for (int i=1; i<=5; i++){
		roodataarray[i] = (RooDataSet*)f1->Get(Form("roodata_eta_%i",i-1));
		xposition[i] = i+1; 
		xposition_err[i] = 0;
	}

	xposition[0] = 1;
	xposition_err[0] = 0;
	roodataarray[0] = (RooDataSet*)f1->Get("roodata_eta_10");

	x = (RooRealVar*)roodataarray[0]->get()->find("roomass");
	x->setBinning(RooBinning(10000,60,120),"cache");
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



	for (int i=0; i<6; i++){
		RooFitResult* fitresult;
		RooFitResult* fitresult1;
		//RooFit::SumW2Error(true)

		//RooFit::AsymptoticError(false)
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

		if (fittype == 1)c1->SaveAs(Form("./fitHI/HI_18_etacut_%i.png",i));
		if (fittype == 2)c1->SaveAs(Form("./fitHI/HI_18_etacut_expbk_%i.png",i));
	}

	TGraphErrors *g_1 = new TGraphErrors(6,xposition,dMass,xposition_err,dMass_Err);
	TGraphErrors *g_2 = new TGraphErrors(6,xposition,dWidth,xposition_err,dWidth_Err);
	
	TFile *f2 = new TFile("pp_data_file_stability.root","UPDATE");
	f2->cd();
	if (fittype == 1){
		g_1->Write("HI_18_etacut_dM",2);
		g_2->Write("HI_18_etacut_dWidth",2);
	}
	if (fittype == 2){
		g_1->Write("HI_18_etacut_expbk_dM",2);
		g_2->Write("HI_18_etacut_expbk_dWidth",2);
	}
	f1->Close();
	f2->Close();

}