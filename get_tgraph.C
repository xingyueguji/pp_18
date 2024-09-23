	void changecosmetic(TGraphErrors* x1,int type){
		x1->SetTitle("");
		if (type == 1) x1->GetYaxis()->SetTitle("dM(GeV)");
		if (type == 2) x1->GetYaxis()->SetTitle("dWidth(GeV)");
		x1->GetXaxis()->SetTitle("Run Period");
		x1->GetXaxis()->CenterTitle();
		x1->GetYaxis()->CenterTitle();

		//if (type == 1) x1->GetYaxis()->SetRangeUser(-0.2,0.1);
		//if (type == 2) x1->GetYaxis()->SetRangeUser(-0.2,0.1);

		//x1->SetLineWidth(2);
		//x1->SetMarkerSize(1.3);

		x1->GetXaxis()->SetLimits(0, 30);
		x1->GetXaxis()->SetNdivisions(30,0,0,kFALSE);


		for (int i=1; i <= 31; ++i){
			x1->GetXaxis()->ChangeLabel(i,-1,0,-1,-1,-1,"");
		}
		x1->GetXaxis()->ChangeLabel(2,-1,0.03,-1,-1,-1,"2018-07-30");
		x1->GetXaxis()->ChangeLabel(21,-1,0.03,-1,-1,-1,"2018-10");
		x1->GetXaxis()->ChangeLabel(24,45,0.03,-1,-1,-1,"inclusive");
		x1->GetXaxis()->ChangeLabel(25,45,0.03,-1,-1,-1,"0-10%");
		x1->GetXaxis()->ChangeLabel(26,45,0.03,-1,-1,-1,"10-20%");
		x1->GetXaxis()->ChangeLabel(27,45,0.03,-1,-1,-1,"20-30%");
		x1->GetXaxis()->ChangeLabel(28,45,0.03,-1,-1,-1,"30-50%");
		x1->GetXaxis()->ChangeLabel(29,45,0.03,-1,-1,-1,"50-100%");

	}

void get_tgraph(int type = 1){

	//type == 1 means no bk sub
	//Type == 2 means bk sub

	TFile *f1 = new TFile("pp_data_file_stability.root","READ");
	TGraphErrors *HI_dm;
	TGraphErrors *HI_dw;
	TGraphErrors *HI_hist_dm;
	TGraphErrors *HI_hist_dw;
	TGraphErrors *HI_chi2_dm;
	TGraphErrors *HI_chi2_dw;

	TGraphErrors *pp_dm;
	TGraphErrors *pp_dw;
	TGraphErrors *pp_hist_dm;
	TGraphErrors *pp_hist_dw;
	TGraphErrors *pp_chi2_dm;
	TGraphErrors *pp_chi2_dw;


	if (type == 1){
		HI_dm = (TGraphErrors*)f1->Get("HI_18_etacut_dM");
		HI_dw = (TGraphErrors*)f1->Get("HI_18_etacut_dWidth");
		HI_hist_dm = (TGraphErrors*)f1->Get("HI_dM_hist");
		HI_hist_dw = (TGraphErrors*)f1->Get("HI_dWidth_hist");
		HI_chi2_dm = (TGraphErrors*)f1->Get("HI_dM_chi2");
		HI_chi2_dw = (TGraphErrors*)f1->Get("HI_dWidth_chi2");

		pp_dm = (TGraphErrors*)f1->Get("runperiod_22_dM");
		pp_dw = (TGraphErrors*)f1->Get("runperiod_22_dWidth");
		pp_hist_dm = (TGraphErrors*)f1->Get("run_period_dM_hist");
		pp_hist_dw = (TGraphErrors*)f1->Get("run_period_dWidth_hist");
		pp_chi2_dm = (TGraphErrors*)f1->Get("run_period_dM_chi2");
		pp_chi2_dw = (TGraphErrors*)f1->Get("run_period_dWidth_chi2");

	}
	if (type == 2){
		HI_dm = (TGraphErrors*)f1->Get("HI_18_etacut_expbk_dM");
		HI_dw = (TGraphErrors*)f1->Get("HI_18_etacut_expbk_dWidth");
		HI_hist_dm = (TGraphErrors*)f1->Get("HI_dM_bksub_hist");
		HI_hist_dw = (TGraphErrors*)f1->Get("HI_dWidth_bksub_hist");
		HI_chi2_dm = (TGraphErrors*)f1->Get("HI_dM_chi2_bksub");
		HI_chi2_dw = (TGraphErrors*)f1->Get("HI_dWidth_chi2_bksub");

		pp_dm = (TGraphErrors*)f1->Get("runperiod_22_expbk_dM");
		pp_dw = (TGraphErrors*)f1->Get("runperiod_22_expbk_dWidth");
		pp_hist_dm = (TGraphErrors*)f1->Get("run_period_dM_hist_bksub");
		pp_hist_dw = (TGraphErrors*)f1->Get("run_period_dWidth_hist_bksub");
		pp_chi2_dm = (TGraphErrors*)f1->Get("run_period_dM_chi2_bksub");
		pp_chi2_dw = (TGraphErrors*)f1->Get("run_period_dWidth_chi2_bksub");

	}

	TGraphErrors* Combined_fit_dm = new TGraphErrors(pp_dm->GetN() + HI_dm->GetN());
	TGraphErrors* Combined_hist_dm = new TGraphErrors(pp_dm->GetN() + HI_dm->GetN());
	TGraphErrors* Combined_chi2_dm = new TGraphErrors(pp_dm->GetN() + HI_dm->GetN());

	TGraphErrors* Combined_fit_dwidth = new TGraphErrors(HI_dw->GetN() + pp_dw->GetN());
	TGraphErrors* Combined_hist_dwidth = new TGraphErrors(HI_dw->GetN() + pp_dw->GetN());
	TGraphErrors* Combined_chi2_dwidth = new TGraphErrors(HI_dw->GetN() + pp_dw->GetN());


    for (int i = 0; i < pp_dm->GetN(); ++i) {
        Combined_fit_dm->SetPoint(i, pp_dm->GetX()[i], pp_dm->GetY()[i]);
        Combined_fit_dm->SetPointError(i, pp_dm->GetErrorX(i), pp_dm->GetErrorY(i));

		Combined_hist_dm->SetPoint(i, pp_hist_dm->GetX()[i], pp_hist_dm->GetY()[i]);
        Combined_hist_dm->SetPointError(i, pp_hist_dm->GetErrorX(i), pp_hist_dm->GetErrorY(i));

		Combined_chi2_dm->SetPoint(i, pp_chi2_dm->GetX()[i], pp_chi2_dm->GetY()[i]);
		Combined_chi2_dm->SetPointError(i, pp_chi2_dm->GetErrorX(i), pp_chi2_dm->GetErrorY(i));
		
		Combined_fit_dwidth->SetPoint(i, pp_dw->GetX()[i], pp_dw->GetY()[i]);
		Combined_fit_dwidth->SetPointError(i, pp_dw->GetErrorX(i), pp_dw->GetErrorY(i));

		Combined_hist_dwidth->SetPoint(i, pp_hist_dw->GetX()[i], pp_hist_dw->GetY()[i]);
		Combined_hist_dwidth->SetPointError(i, pp_hist_dw->GetErrorX(i), pp_hist_dw->GetErrorY(i));

		Combined_chi2_dwidth->SetPoint(i, pp_chi2_dw->GetX()[i], pp_chi2_dw->GetY()[i]);
		Combined_chi2_dwidth->SetPointError(i, pp_chi2_dw->GetErrorX(i), pp_chi2_dw->GetErrorY(i));

		
    }
    for (int i = 0; i < HI_dm->GetN(); ++i) {
        Combined_fit_dm->SetPoint(i + pp_dm->GetN(), HI_dm->GetX()[i] + 22, HI_dm->GetY()[i]);
        Combined_fit_dm->SetPointError(i + pp_dm->GetN(), HI_dm->GetErrorX(i), HI_dm->GetErrorY(i));

		Combined_hist_dm->SetPoint(i + pp_hist_dm->GetN(), HI_hist_dm->GetX()[i] + 22, HI_hist_dm->GetY()[i]);
        Combined_hist_dm->SetPointError(i + pp_hist_dm->GetN(), HI_hist_dm->GetErrorX(i), HI_hist_dm->GetErrorY(i));

		Combined_chi2_dm->SetPoint(i + pp_chi2_dm->GetN(), HI_chi2_dm->GetX()[i] + 22, HI_chi2_dm->GetY()[i]);
		Combined_chi2_dm->SetPointError(i + pp_chi2_dm->GetN(), HI_chi2_dm->GetErrorX(i), HI_chi2_dm->GetErrorY(i));

		Combined_fit_dwidth->SetPoint(i + pp_dw->GetN(), HI_dw->GetX()[i] + 22, HI_dw->GetY()[i]);
        Combined_fit_dwidth->SetPointError(i + pp_dw->GetN(), HI_dw->GetErrorX(i), HI_dw->GetErrorY(i));

		Combined_hist_dwidth->SetPoint(i + pp_hist_dw->GetN(), HI_hist_dw->GetX()[i] + 22, HI_hist_dw->GetY()[i]);
		Combined_hist_dwidth->SetPointError(i + pp_hist_dw->GetN(), HI_hist_dw->GetErrorX(i), HI_hist_dw->GetErrorY(i));

		Combined_chi2_dwidth->SetPoint(i + pp_chi2_dw->GetN(), HI_chi2_dw->GetX()[i] + 22, HI_chi2_dw->GetY()[i]);
		Combined_chi2_dwidth->SetPointError(i + pp_chi2_dw->GetN(), HI_chi2_dw->GetErrorX(i), HI_chi2_dw->GetErrorY(i));
    }

	gStyle->SetEndErrorSize(4);
	gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.11);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
	gStyle->SetTitleFont(42);
	gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.03, "XYZ");



	TCanvas *c1 = new TCanvas("","",1600,800);

	c1->cd();

	changecosmetic(Combined_fit_dm,1);
	changecosmetic(Combined_hist_dm,1);

	Combined_fit_dm->SetMarkerColor(2);
	Combined_fit_dm->SetLineColor(2);
	Combined_fit_dm->SetMarkerStyle(20);
	Combined_fit_dm->SetMarkerSize(1);

	Combined_hist_dm->SetMarkerColor(4);
	Combined_hist_dm->SetLineColor(4);
	Combined_hist_dm->SetMarkerStyle(21);
	Combined_hist_dm->SetMarkerSize(1);

	Combined_chi2_dm->SetMarkerColor(3);
	Combined_chi2_dm->SetLineColor(3);
	Combined_chi2_dm->SetMarkerStyle(4);
	Combined_chi2_dm->SetMarkerSize(1);
	

	Combined_fit_dm->GetYaxis()->SetRangeUser(-2.2,2);

	Combined_fit_dm->Draw("AP");
	Combined_hist_dm->Draw("PSAME");
	Combined_chi2_dm->Draw("PSAME");

	TF1 *fit_pp_fit = new TF1("fit_pp_fit", "[0]", 0, 22.5);
	TF1 *fit_pp_hist = new TF1("fit_pp_hist", "[0]", 0, 22.5);
	TF1 *fit_pp_chi2 = new TF1("fit_pp_chi2", "[0]", 0, 22.5);

	TF1 *fit_HI_fit = new TF1("fit_HI_fit", "[0]", 0, 80);
	TF1 *fit_HI_hist = new TF1("fit_HI_hist", "[0]", 0, 80);
	TF1 *fit_HI_chi2 = new TF1("fit_HI_chi2", "[0]", 0, 80);

	fit_pp_fit->SetLineColor(2);
	fit_pp_hist->SetLineColor(4);
	fit_pp_chi2->SetLineColor(3);
	fit_HI_fit->SetLineColor(2);
	fit_HI_hist->SetLineColor(4);
	fit_HI_chi2->SetLineColor(3);

	Combined_fit_dm->Fit(fit_pp_fit,"QR","",0.5,22.5);
	Combined_hist_dm->Fit(fit_pp_hist,"QR","",0.5,22.5);
	Combined_chi2_dm->Fit(fit_pp_chi2,"QR","",0.5,22.5);

	Combined_fit_dm->Fit(fit_HI_fit,"QR","",22.6,28.5);
	Combined_hist_dm->Fit(fit_HI_hist,"QR","",22.6,28.5);
	Combined_chi2_dm->Fit(fit_HI_chi2,"QR","",22.6,28.5);

	fit_pp_fit->Draw("SAME");
	fit_pp_hist->Draw("SAME");
	fit_pp_chi2->Draw("SAME");


	// Create a TLine to draw on the graph
	TLine *line = new TLine(22.5, -2.2, 22.5, 2); // Define the start (1, 4.0) and end (5, 4.0) points of the line

	// Set the line style to dashed
	line->SetLineStyle(2); // 2 corresponds to a dashed line style
	line->SetLineColor(kBlack); // Optional: Set the color of the line	
	line->SetLineWidth(2); // Optional: Set the width of the line

	// Draw the line on the canvas
	line->Draw("same"); // "same" ensures the line is drawn on the current canvas


	TLegend *legend = new TLegend(0.6, 0.75, 0.75, 0.9);
	legend->SetTextFont(40);
	legend->SetBorderSize(0);
    if (type == 1) legend->AddEntry(Combined_fit_dm, "unbinned_fit", "lep"); 
	if (type == 2) legend->AddEntry(Combined_fit_dm, "unbinned_fit_exp_bk", "lep"); 
	if (type == 1)legend->AddEntry(Combined_hist_dm, "hist_stat", "lep"); 
	if (type == 2)legend->AddEntry(Combined_hist_dm, "hist_stat_exp_bk", "lep"); 
	if (type == 1)legend->AddEntry(Combined_chi2_dm, "chi2", "lep"); 
	if (type == 2)legend->AddEntry(Combined_chi2_dm, "chi2_bksub", "lep"); 
    legend->SetFillColor(0);

	legend->Draw();

	TPaveText *pt = new TPaveText(0.2, 0.65, 0.55, 0.9, "brNDC"); // normalized coordinates
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextAlign(12); // Align left and vertically centered
	pt->SetTextFont(42);
	pt->SetTextSize(0.03);
	pt->SetMargin(0.01);
	pt->AddText(Form("pp_fit fit = %.4f #pm %.4f", fit_pp_fit->GetParameter(0),fit_pp_fit->GetParError(0)));
	pt->AddText(Form("HI_fit fit = %.4f #pm %.4f", fit_HI_fit->GetParameter(0),fit_HI_fit->GetParError(0)));
	pt->AddText(Form("pp_hist fit = %.4f #pm %.4f", fit_pp_hist->GetParameter(0),fit_pp_hist->GetParError(0)));
	pt->AddText(Form("HI_hist fit = %.4f #pm %.4f", fit_HI_hist->GetParameter(0),fit_HI_hist->GetParError(0)));
	pt->AddText(Form("pp_chi2 fit = %.4f #pm %.4f", fit_pp_chi2->GetParameter(0),fit_pp_hist->GetParError(0)));
	pt->AddText(Form("HI_chi2 fit = %.4f #pm %.4f", fit_HI_chi2->GetParameter(0),fit_HI_hist->GetParError(0)));
	pt->Draw();




	//Below are width

	TCanvas *c2 = new TCanvas("","",1600,800);
	
	c2->cd();

	changecosmetic(Combined_fit_dwidth,2);
	changecosmetic(Combined_hist_dwidth,2);

	Combined_fit_dwidth->GetYaxis()->SetRangeUser(-1,2);

	Combined_fit_dwidth->SetMarkerColor(2);
	Combined_fit_dwidth->SetLineColor(2);
	Combined_fit_dwidth->SetMarkerStyle(20);
	Combined_fit_dwidth->SetMarkerSize(1);

	Combined_hist_dwidth->SetMarkerColor(4);
	Combined_hist_dwidth->SetLineColor(4);
	Combined_hist_dwidth->SetMarkerStyle(21);
	Combined_hist_dwidth->SetMarkerSize(1);

	Combined_chi2_dwidth->SetMarkerColor(3);
	Combined_chi2_dwidth->SetLineColor(3);
	Combined_chi2_dwidth->SetMarkerStyle(4);
	Combined_chi2_dwidth->SetMarkerSize(1);

	Combined_fit_dwidth->Draw("AP");
	Combined_hist_dwidth->Draw("PSAME");
	Combined_chi2_dwidth->Draw("PSAME");

	TF1 *fit_pp_fit_width = new TF1("fit_pp_fit_width", "[0]", 0, 22.5);
	TF1 *fit_pp_hist_width = new TF1("fit_pp_hist_width", "[0]", 0, 22.5);
	TF1 *fit_pp_chi2_width = new TF1("fit_pp_chi2_width", "[0]", 0, 22.5);

	TF1 *fit_HI_fit_width = new TF1("fit_HI_fit_width", "[0]", 0, 80);
	TF1 *fit_HI_hist_width = new TF1("fit_HI_hist_width", "[0]", 0, 80);
	TF1 *fit_HI_chi2_width = new TF1("fit_HI_chi2_width", "[0]", 0, 80);


	fit_pp_fit_width->SetLineColor(2);
	fit_pp_hist_width->SetLineColor(4);
	fit_pp_chi2_width->SetLineColor(3);
	fit_HI_fit_width->SetLineColor(2);
	fit_HI_hist_width->SetLineColor(4);
	fit_HI_chi2_width->SetLineColor(3);


	Combined_fit_dwidth->Fit(fit_pp_fit_width,"QR","",0.5,22.5);
	Combined_hist_dwidth->Fit(fit_pp_hist_width,"QR","",0.5,22.5);
	Combined_chi2_dwidth->Fit(fit_pp_chi2_width,"QR","",0.5,22.5);

	Combined_fit_dwidth->Fit(fit_HI_fit_width,"QR","",22.6,28.5);
	Combined_hist_dwidth->Fit(fit_HI_hist_width,"QR","",22.6,28.5);
	Combined_chi2_dwidth->Fit(fit_HI_chi2_width,"QR","",22.6,28.5);

	fit_pp_fit_width->Draw("SAME");
	fit_pp_hist_width->Draw("SAME");
	fit_pp_chi2_width->Draw("SAME");


	TLegend *legend_wdith = new TLegend(0.6, 0.75, 0.75, 0.9);
	legend_wdith->SetTextFont(40);
	legend_wdith->SetBorderSize(0);
    if (type == 1) legend_wdith->AddEntry(Combined_fit_dwidth, "unbinned_fit", "lep"); 
	if (type == 2) legend_wdith->AddEntry(Combined_fit_dwidth, "unbinned_fit_exp", "lep"); 
	if (type == 1) legend_wdith->AddEntry(Combined_hist_dwidth, "hist_stat", "lep"); 
	if (type == 2) legend_wdith->AddEntry(Combined_hist_dwidth, "hist_stat_exp", "lep"); 
	if (type == 1) legend_wdith->AddEntry(Combined_chi2_dwidth, "chi2", "lep"); 
	if (type == 2) legend_wdith->AddEntry(Combined_chi2_dwidth, "chi2_bksub", "lep"); 
    legend_wdith->SetFillColor(0);

	TPaveText *pt_1;
	if (type == 1) pt_1 = new TPaveText(0.2, 0.75, 0.55, 0.9, "brNDC"); 
	else{
		pt_1 = new TPaveText(0.2, 0.65, 0.55, 0.9, "brNDC");
	}
	pt_1->SetBorderSize(0);
	pt_1->SetFillColor(0);
	pt_1->SetTextAlign(12); // Align left and vertically centered
	pt_1->SetTextFont(42);
	pt_1->SetTextSize(0.03);
	pt_1->SetMargin(0.01);
	pt_1->AddText(Form("pp_fit fit = %.4f #pm %.4f", fit_pp_fit_width->GetParameter(0),fit_pp_fit_width->GetParError(0)));
	pt_1->AddText(Form("HI_fit fit = %.4f #pm %.4f", fit_HI_fit_width->GetParameter(0),fit_HI_fit_width->GetParError(0)));
	pt_1->AddText(Form("pp_hist fit = %.4f #pm %.4f + 3", fit_pp_hist_width->GetParameter(0),fit_pp_hist_width->GetParError(0)));
	pt_1->AddText(Form("HI_hist fit = %.4f #pm %.4f + 3", fit_HI_hist_width->GetParameter(0),fit_HI_hist_width->GetParError(0)));
	pt_1->AddText(Form("pp_chi2 fit = %.4f #pm %.4f", fit_pp_chi2_width->GetParameter(0),fit_pp_chi2_width->GetParError(0)));
	pt_1->AddText(Form("HI_chi2 fit = %.4f #pm %.4f", fit_HI_chi2_width->GetParameter(0),fit_HI_chi2_width->GetParError(0)));
	pt_1->Draw();
	legend_wdith->Draw();

	TLine *line_width = new TLine(22.5, -1, 22.5, 2); // Define the start (1, 4.0) and end (5, 4.0) points of the line

	// Set the line style to dashed
	line_width->SetLineStyle(2); // 2 corresponds to a dashed line style
	line_width->SetLineColor(kBlack); // Optional: Set the color of the line	
	line_width->SetLineWidth(2); // Optional: Set the width of the line

	// Draw the line on the canvas
	line_width->Draw("same"); // "same" ensures the line is drawn on the current canvas
	



	



	if (type == 1){
		c1->SaveAs("Combined_dm.png");
		c2->SaveAs("Combined_dwidth.png");
	}
	else{
		c1->SaveAs("Combined_exp_dm.png");
		c2->SaveAs("Combined_exp_dwidth.png");
	}

	




}