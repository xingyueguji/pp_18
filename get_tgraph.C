#include "TH1.h"
#include "TCanvas.h"
#include "tdrStyle.C"
#include "CMS_lumi.C"

void getdiffcosmetic(TCanvas *c1, Double_t *arr, Double_t *arrerr, bool isdM)
{
	c1->cd();
	setTDRStyle();

	double y[3] = {1, 2, 3};
	double y_err[3] = {0, 0, 0};

	TGraphErrors *g1 = new TGraphErrors(3, arr, y, arrerr, y_err);

	c1->SetTitle("");
	g1->SetTitle("");

	if (isdM)
		g1->GetXaxis()->SetTitle("#DeltaM (MeV)");
	else
		g1->GetXaxis()->SetTitle("#DeltaWidth (MeV)");

	if (isdM)
		g1->GetXaxis()->SetLimits(-130, 10);

	// g1->GetXaxis()->CenterTitle();
	g1->GetYaxis()->SetLimits(0, 5);
	g1->GetYaxis()->SetRangeUser(0, 4);
	g1->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

	for (int i = 1; i <= 5; ++i)
	{
		g1->GetYaxis()->ChangeLabel(i, -1, 0, -1, -1, -1, "");
	}

	g1->GetYaxis()->ChangeLabel(2, -1, 0.03, -1, -1, -1, "Fit");
	g1->GetYaxis()->ChangeLabel(3, -1, 0.03, -1, -1, -1, "Hist");
	g1->GetYaxis()->ChangeLabel(4, -1, 0.03, -1, -1, -1, "Template");

	g1->SetMarkerColor(kBlack);
	g1->SetLineColor(kBlack);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(1);

	g1->Draw("AP");

	TLine *line = new TLine(0, 0, 0, 4); // Define the start (1, 4.0) and end (5, 4.0) points of the line

	// Set the line style to dashed
	line->SetLineStyle(2);		// 2 corresponds to a dashed line style
	line->SetLineColor(kBlack); // Optional: Set the color of the line
	line->SetLineWidth(2);		// Optional: Set the width of the line

	// Draw the line on the canvas
	line->Draw("same");

	TPaveText *pave = new TPaveText(0.65, 0.75, 0.9, 0.9, "NDC");
	pave->AddText("Inclusive HI - Fitted pp");
	pave->SetTextFont(42);
	pave->SetTextSize(0.04);
	pave->SetFillColor(0);
	pave->SetBorderSize(0);
	pave->Draw("same");

	CMS_lumi(c1, 13, 11);
}
void changecosmetic(TCanvas *c1, TGraphErrors *x1, int axistype, int bktype, TGraphErrors *hist = nullptr, TGraphErrors *chi2 = nullptr, Double_t *arr = nullptr, Double_t *arrerr = nullptr)
{
	c1->cd();
	x1->SetTitle("");
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	TString yaxistitle;
	Double_t yaxisupperlimit;
	Double_t yaxislowerlimit;

	if (axistype == 1)
	{
		yaxistitle = "#DeltaM (GeV)";
		yaxisupperlimit = 0.7;
		yaxislowerlimit = -1.5;
	}

	if (axistype == 2)
	{
		yaxistitle = "#DeltaWidth (GeV)";
		yaxisupperlimit = 2;
		yaxislowerlimit = -1;
	}
	if (axistype == 3)
	{
		yaxistitle = "Alpha";
		yaxisupperlimit = 2.2;
		yaxislowerlimit = 1.6;
	}
	if (axistype == 4)
	{
		yaxistitle = "N";
		yaxisupperlimit = 1.2;
		yaxislowerlimit = 0.7;
	}
	if (axistype == 5)
	{
		yaxistitle = "STD";
		yaxisupperlimit = 1.4;
		yaxislowerlimit = 0;
	}

	x1->GetYaxis()->SetTitle(yaxistitle);
	x1->GetXaxis()->SetTitle("Run Period");
	x1->GetXaxis()->CenterTitle();
	x1->GetYaxis()->CenterTitle();
	x1->GetYaxis()->SetRangeUser(yaxislowerlimit, yaxisupperlimit);
	x1->GetXaxis()->SetLimits(0, 30);
	x1->GetXaxis()->SetNdivisions(30, 0, 0, kFALSE);

	for (int i = 1; i <= 31; ++i)
	{
		x1->GetXaxis()->ChangeLabel(i, -1, 0, -1, -1, -1, "");
	}
	x1->GetXaxis()->ChangeLabel(2, -1, 0.03, -1, -1, -1, "2018-07-30");
	x1->GetXaxis()->ChangeLabel(21, -1, 0.03, -1, -1, -1, "2018-10");
	x1->GetXaxis()->ChangeLabel(24, 45, 0.03, -1, -1, -1, "0-100%");
	x1->GetXaxis()->ChangeLabel(25, 45, 0.03, -1, -1, -1, "0-10%");
	x1->GetXaxis()->ChangeLabel(26, 45, 0.03, -1, -1, -1, "10-20%");
	x1->GetXaxis()->ChangeLabel(27, 45, 0.03, -1, -1, -1, "20-30%");
	x1->GetXaxis()->ChangeLabel(28, 45, 0.03, -1, -1, -1, "30-100%");
	x1->GetXaxis()->SetLabelOffset(0.03);
	// x1->GetXaxis()->ChangeLabel(29, 45, 0.03, -1, -1, -1, "50-100%");

	x1->SetMarkerColor(2);
	x1->SetLineColor(2);
	x1->SetMarkerStyle(20);
	x1->SetMarkerSize(1);
	x1->SetFillColorAlpha(2, 0.5);

	if (hist != nullptr)
	{
		hist->SetMarkerColor(4);
		hist->SetLineColor(4);
		hist->SetMarkerStyle(21);
		hist->SetMarkerSize(1);
		hist->SetFillColorAlpha(4, 0.5);
	}

	if (chi2 != nullptr)
	{
		chi2->SetMarkerColor(kGreen + 3);
		chi2->SetLineColor(kGreen + 3);
		chi2->SetFillColorAlpha(kGreen + 3, 0.35);
		chi2->SetMarkerStyle(4);
		chi2->SetMarkerSize(1);
	}

	x1->Draw("AP");
	// x1->Draw("P SAME");
	// x1->SetFillStyle(3001);

	if (hist != nullptr)
	{
		// hist->SetFillStyle(3001);
		hist->Draw("P SAME");
		// hist->Draw("P SAME");
	}
	if (chi2 != nullptr)
	{
		// chi2->SetFillStyle(3005);
		chi2->Draw("P SAME");
		// chi2->Draw("P SAME");
	}

	TF1 *fit_pp_fit = new TF1("fit_pp_fit", "[0]", 0, 22.5);
	TF1 *fit_pp_hist = new TF1("fit_pp_hist", "[0]", 0, 22.5);
	TF1 *fit_pp_chi2 = new TF1("fit_pp_template", "[0]", 0, 22.5);

	TF1 *fit_HI_fit = new TF1("fit_HI_fit", "[0]", 0, 80);
	TF1 *fit_HI_hist = new TF1("fit_HI_hist", "[0]", 0, 80);
	TF1 *fit_HI_chi2 = new TF1("fit_HI_template", "[0]", 0, 80);

	fit_pp_fit->SetLineColor(2);
	fit_pp_hist->SetLineColor(4);
	fit_pp_chi2->SetLineColor(kGreen + 3);
	fit_HI_fit->SetLineColor(2);
	fit_HI_hist->SetLineColor(4);
	fit_HI_chi2->SetLineColor(kGreen + 3);

	x1->Fit(fit_pp_fit, "QR", "", 0.5, 22.5);
	if (hist != nullptr)
	{
		hist->Fit(fit_pp_hist, "QR", "", 0.5, 22.5);
	}
	if (chi2 != nullptr)
	{
		chi2->Fit(fit_pp_chi2, "QR", "", 0.5, 22.5);
	}

	x1->Fit(fit_HI_fit, "QR", "", 23.6, 28.5);
	if (hist != nullptr)
	{
		hist->Fit(fit_HI_hist, "QR", "", 23.6, 28.5);
	}
	if (chi2 != nullptr)
	{
		chi2->Fit(fit_HI_chi2, "QR", "", 23.6, 28.5);
	}

	fit_pp_fit->Draw("SAME");
	if (hist != nullptr)
	{
		fit_pp_hist->Draw("SAME");
	}
	if (chi2 != nullptr)
	{
		fit_pp_chi2->Draw("SAME");
	}

	// Create a TLine to draw on the graph
	TLine *line = new TLine(22.5, yaxislowerlimit, 22.5, yaxisupperlimit); // Define the start (1, 4.0) and end (5, 4.0) points of the line

	// Set the line style to dashed
	line->SetLineStyle(2);		// 2 corresponds to a dashed line style
	line->SetLineColor(kBlack); // Optional: Set the color of the line
	line->SetLineWidth(2);		// Optional: Set the width of the line

	// Draw the line on the canvas
	line->Draw("same"); // "same" ensures the line is drawn on the current canvas

	TLegend *legend;

	if (bktype == 1)
		legend = new TLegend(0.6, 0.75, 0.75, 0.9);
	if (bktype == 2)
		legend = new TLegend(0.6, 0.75, 0.85, 0.9);
	legend->SetTextFont(40);
	legend->SetBorderSize(0);
	if (bktype == 1)
		legend->AddEntry(x1, "unbinned_fit", "lep");
	if (bktype == 2)
		legend->AddEntry(x1, "unbinned_fit_exp_bk", "lep");
	if (bktype == 1 && hist != nullptr)
		legend->AddEntry(hist, "hist_stat", "lep");
	if (bktype == 2 && hist != nullptr)
		legend->AddEntry(hist, "hist_stat_exp_bk", "lep");
	if (bktype == 1 && chi2 != nullptr)
		legend->AddEntry(chi2, "template_fit", "lep");
	if (bktype == 2 && chi2 != nullptr)
		legend->AddEntry(chi2, "template_fit - template background", "lep");
	legend->SetFillColor(0);

	legend->Draw();

	TPaveText *pt = new TPaveText(0.2, 0.7, 0.55, 0.9, "brNDC"); // normalized coordinates
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextAlign(12); // Align left and vertically centered
	pt->SetTextFont(42);
	pt->SetTextSize(0.02);
	pt->SetMargin(0.01);
	pt->AddText(Form("pp_fit fit = %.4f #pm %.4f", fit_pp_fit->GetParameter(0), fit_pp_fit->GetParError(0)));
	pt->AddText(Form("HI_fit fit = %.4f #pm %.4f", fit_HI_fit->GetParameter(0), fit_HI_fit->GetParError(0)));

	if (arr != nullptr && arrerr != nullptr)
	{
		double pp_result_fit = fit_pp_fit->GetParameter(0);
		double pp_result_fit_err = fit_pp_fit->GetParError(0);

		double pp_result_temp = fit_pp_chi2->GetParameter(0);
		double pp_result_temp_err = fit_pp_chi2->GetParError(0);

		double pp_result_hist = fit_pp_hist->GetParameter(0);
		double pp_result_hist_err = fit_pp_hist->GetParError(0);

		double HI_inclusive_fit_x;
		double HI_inclusive_fit_y;
		double HI_inclusive_hist_x;
		double HI_inclusive_hist_y;
		double HI_inclusive_temp_x;
		double HI_inclusive_temp_y;

		x1->GetPoint(22, HI_inclusive_fit_x, HI_inclusive_fit_y);
		hist->GetPoint(22, HI_inclusive_hist_x, HI_inclusive_hist_y);
		chi2->GetPoint(22, HI_inclusive_temp_x, HI_inclusive_temp_y);

		if (HI_inclusive_fit_x != 23)
			cout << "You selected the wrong point for inclusive" << endl;

		double HI_inclusive_fit_err = x1->GetErrorY(22);
		double HI_inclusive_hist_err = hist->GetErrorY(22);
		double HI_inclusive_temp_err = chi2->GetErrorY(22);

		double fit_diff = HI_inclusive_fit_y - pp_result_fit;
		double fit_diff_err = sqrt(TMath::Power(HI_inclusive_fit_err, 2) + TMath::Power(pp_result_fit_err, 2));
		double hist_diff = HI_inclusive_hist_y - pp_result_hist;
		double hist_diff_err = sqrt(TMath::Power(HI_inclusive_hist_err, 2) + TMath::Power(pp_result_hist_err, 2));
		double temp_diff = HI_inclusive_temp_y - pp_result_temp;
		double temp_diff_err = sqrt(TMath::Power(HI_inclusive_temp_err, 2) + TMath::Power(pp_result_temp_err, 2));

		arr[0] = fit_diff * 1000;
		arr[1] = hist_diff * 1000;
		arr[2] = temp_diff * 1000;

		arrerr[0] = fit_diff_err * 1000;
		arrerr[1] = hist_diff_err * 1000;
		arrerr[2] = temp_diff_err * 1000;

		cout << "The shifted diff for fit is " << fit_diff << " Error is " << fit_diff_err << endl;
		cout << "The shifted diff for hist is " << hist_diff << " Error is " << hist_diff_err << endl;
		cout << "The shifted diff for temp is " << temp_diff << " Error is " << temp_diff_err << endl;
	}

	if (hist != nullptr)
	{
		pt->AddText(Form("pp_hist fit = %.4f #pm %.4f", fit_pp_hist->GetParameter(0), fit_pp_hist->GetParError(0)));
		pt->AddText(Form("HI_hist fit = %.4f #pm %.4f", fit_HI_hist->GetParameter(0), fit_HI_hist->GetParError(0)));
	}
	if (chi2 != nullptr)
	{
		pt->AddText(Form("pp_template fit = %.4f #pm %.4f", fit_pp_chi2->GetParameter(0), fit_pp_chi2->GetParError(0)));
		pt->AddText(Form("HI_template fit = %.4f #pm %.4f", fit_HI_chi2->GetParameter(0), fit_HI_chi2->GetParError(0)));
	}

	pt->Draw();
}

void get_tgraph(int type = 2)
{

	// setTDRStyle();

	// type == 1 means no bk sub
	// Type == 2 means bk sub
	TFile *f1 = new TFile("All_plots.root", "READ");

	TGraphErrors *Diff_dM_raw;
	TGraphErrors *Diff_dW_raw;

	TGraphErrors *Diff_dM_eta;
	TGraphErrors *Diff_dW_eta;

	Double_t Diff_dM_array_raw[3];
	Double_t Diff_dM_array_err_raw[3];

	Double_t Diff_dWidth_array_raw[3];
	Double_t Diff_dWidth_array_err_raw[3];

	Double_t Diff_dM_array_eta[3];
	Double_t Diff_dM_array_err_eta[3];

	Double_t Diff_dWidth_array_eta[3];
	Double_t Diff_dWidth_array_err_eta[3];

	TGraphErrors *HI_dm_raw;
	TGraphErrors *HI_dwidth_raw;
	TGraphErrors *HI_alpha_raw;
	TGraphErrors *HI_n_raw;
	TGraphErrors *HI_STD_raw;

	TGraphErrors *HI_dm_eta;
	TGraphErrors *HI_dwidth_eta;
	TGraphErrors *HI_alpha_eta;
	TGraphErrors *HI_n_eta;
	TGraphErrors *HI_STD_eta;

	TGraphErrors *HI_hist_dm_raw;
	TGraphErrors *HI_hist_dm_eta;
	TGraphErrors *HI_hist_dw_raw;
	TGraphErrors *HI_hist_dw_eta;

	TGraphErrors *HI_chi2_dm_raw;
	TGraphErrors *HI_chi2_dw_raw;
	TGraphErrors *HI_chi2_dm_eta;
	TGraphErrors *HI_chi2_dw_eta;

	TGraphErrors *pp_dm_raw;
	TGraphErrors *pp_dwidth_raw;
	TGraphErrors *pp_alpha_raw;
	TGraphErrors *pp_n_raw;
	TGraphErrors *pp_STD_raw;

	TGraphErrors *pp_dm_eta;
	TGraphErrors *pp_dwidth_eta;
	TGraphErrors *pp_alpha_eta;
	TGraphErrors *pp_n_eta;
	TGraphErrors *pp_STD_eta;

	TGraphErrors *pp_hist_dm_raw;
	TGraphErrors *pp_hist_dm_eta;
	TGraphErrors *pp_hist_dw_raw;
	TGraphErrors *pp_hist_dw_eta;

	TGraphErrors *pp_chi2_dm_raw;
	TGraphErrors *pp_chi2_dw_raw;
	TGraphErrors *pp_chi2_dm_eta;
	TGraphErrors *pp_chi2_dw_eta;

	if (type == 1)
	{
		// Here's HI_fit_raw and eta
		HI_dm_raw = (TGraphErrors *)f1->Get("HI_raw_dM");
		HI_dwidth_raw = (TGraphErrors *)f1->Get("HI_raw_dW");
		HI_alpha_raw = (TGraphErrors *)f1->Get("HI_raw_alpha");
		HI_n_raw = (TGraphErrors *)f1->Get("HI_raw_N");
		HI_STD_raw = (TGraphErrors *)f1->Get("HI_raw_sigma");

		HI_dm_eta = (TGraphErrors *)f1->Get("HI_eta_dM");
		HI_dwidth_eta = (TGraphErrors *)f1->Get("HI_eta_dW");
		HI_alpha_eta = (TGraphErrors *)f1->Get("HI_eta_alpha");
		HI_n_eta = (TGraphErrors *)f1->Get("HI_eta_N");
		HI_STD_eta = (TGraphErrors *)f1->Get("HI_eta_sigma");

		// Here's HI hist raw and eta

		HI_hist_dm_raw = (TGraphErrors *)f1->Get("hist_HI_dM_raw_nobksub");
		HI_hist_dm_eta = (TGraphErrors *)f1->Get("hist_HI_dM_eta_nobksub");
		HI_hist_dw_raw = (TGraphErrors *)f1->Get("hist_HI_dWidth_raw_nobksub");
		HI_hist_dw_eta = (TGraphErrors *)f1->Get("hist_HI_dWidth_eta_nobksub");

		// Here's HI chi2 raw and eta

		HI_chi2_dm_raw = (TGraphErrors *)f1->Get("HI_dM_chi2_raw");
		HI_chi2_dw_raw = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw");
		HI_chi2_dm_eta = (TGraphErrors *)f1->Get("HI_dM_chi2_eta");
		HI_chi2_dw_eta = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta");

		// Here's pp_fit_raw and eta
		pp_dm_raw = (TGraphErrors *)f1->Get("pp_raw_dM");
		pp_dwidth_raw = (TGraphErrors *)f1->Get("pp_raw_dW");
		pp_alpha_raw = (TGraphErrors *)f1->Get("pp_raw_alpha");
		pp_n_raw = (TGraphErrors *)f1->Get("pp_raw_N");
		pp_STD_raw = (TGraphErrors *)f1->Get("pp_raw_sigma");

		pp_dm_eta = (TGraphErrors *)f1->Get("pp_eta_dM");
		pp_dwidth_eta = (TGraphErrors *)f1->Get("pp_eta_dW");
		pp_alpha_eta = (TGraphErrors *)f1->Get("pp_eta_alpha");
		pp_n_eta = (TGraphErrors *)f1->Get("pp_eta_N");
		pp_STD_eta = (TGraphErrors *)f1->Get("pp_eta_sigma");

		pp_hist_dm_raw = (TGraphErrors *)f1->Get("hist_pp_dM_raw_nobksub");
		pp_hist_dm_eta = (TGraphErrors *)f1->Get("hist_pp_dM_eta_nobksub");
		pp_hist_dw_raw = (TGraphErrors *)f1->Get("hist_pp_dWidth_raw_nobksub");
		pp_hist_dw_eta = (TGraphErrors *)f1->Get("hist_pp_dWidth_eta_nobksub");

		pp_chi2_dm_raw = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nobksub");
		pp_chi2_dw_raw = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nobksub");
		pp_chi2_dm_eta = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nobksub");
		pp_chi2_dw_eta = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nobksub");
		/*pp_hist_dm = (TGraphErrors *)f1->Get("run_period_dM_hist");
		pp_hist_dw = (TGraphErrors *)f1->Get("run_period_dWidth_hist");
		pp_chi2_dm = (TGraphErrors *)f1->Get("run_period_dM_chi2");
		pp_chi2_dw = (TGraphErrors *)f1->Get("run_period_dWidth_chi2");*/
	}
	if (type == 2)
	{
		HI_dm_raw = (TGraphErrors *)f1->Get("HI_raw_dM_exp");
		HI_dwidth_raw = (TGraphErrors *)f1->Get("HI_raw_dW_exp");
		HI_alpha_raw = (TGraphErrors *)f1->Get("HI_raw_alpha_exp");
		HI_n_raw = (TGraphErrors *)f1->Get("HI_raw_N_exp");
		HI_STD_raw = (TGraphErrors *)f1->Get("HI_raw_sigma_exp");

		HI_dm_eta = (TGraphErrors *)f1->Get("HI_eta_dM_exp");
		HI_dwidth_eta = (TGraphErrors *)f1->Get("HI_eta_dW_exp");
		HI_alpha_eta = (TGraphErrors *)f1->Get("HI_eta_alpha_exp");
		HI_n_eta = (TGraphErrors *)f1->Get("HI_eta_N_exp");
		HI_STD_eta = (TGraphErrors *)f1->Get("HI_eta_sigma_exp");

		HI_hist_dm_raw = (TGraphErrors *)f1->Get("hist_HI_dM_raw_bksub");
		HI_hist_dm_eta = (TGraphErrors *)f1->Get("hist_HI_dM_eta_bksub");
		HI_hist_dw_raw = (TGraphErrors *)f1->Get("hist_HI_dWidth_raw_bksub");
		HI_hist_dw_eta = (TGraphErrors *)f1->Get("hist_HI_dWidth_eta_bksub");

		// Here's HI chi2 raw and eta

		HI_chi2_dm_raw = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_bksub");
		HI_chi2_dw_raw = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_bksub");
		HI_chi2_dm_eta = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_bksub");
		HI_chi2_dw_eta = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_bksub");

		// Below are pp

		pp_dm_raw = (TGraphErrors *)f1->Get("pp_raw_dM_exp");
		pp_dwidth_raw = (TGraphErrors *)f1->Get("pp_raw_dW_exp");
		pp_alpha_raw = (TGraphErrors *)f1->Get("pp_raw_alpha_exp");
		pp_n_raw = (TGraphErrors *)f1->Get("pp_raw_N_exp");
		pp_STD_raw = (TGraphErrors *)f1->Get("pp_raw_sigma_exp");

		pp_dm_eta = (TGraphErrors *)f1->Get("pp_eta_dM_exp");
		pp_dwidth_eta = (TGraphErrors *)f1->Get("pp_eta_dW_exp");
		pp_alpha_eta = (TGraphErrors *)f1->Get("pp_eta_alpha_exp");
		pp_n_eta = (TGraphErrors *)f1->Get("pp_eta_N_exp");
		pp_STD_eta = (TGraphErrors *)f1->Get("pp_eta_sigma_exp");

		pp_hist_dm_raw = (TGraphErrors *)f1->Get("hist_pp_dM_raw_bksub");
		pp_hist_dm_eta = (TGraphErrors *)f1->Get("hist_pp_dM_eta_bksub");
		pp_hist_dw_raw = (TGraphErrors *)f1->Get("hist_pp_dWidth_raw_bksub");
		pp_hist_dw_eta = (TGraphErrors *)f1->Get("hist_pp_dWidth_eta_bksub");

		pp_chi2_dm_raw = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_bksub");
		pp_chi2_dw_raw = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_bksub");
		pp_chi2_dm_eta = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_bksub");
		pp_chi2_dw_eta = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_bksub");
	}

	TGraphErrors *Combined_fit_dm_raw = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_fit_dwidth_raw = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());
	TGraphErrors *Combined_fit_alpha_raw = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());
	TGraphErrors *Combined_fit_n_raw = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());
	TGraphErrors *Combined_fit_STD_raw = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());

	TGraphErrors *Combined_fit_dm_eta = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_fit_dwidth_eta = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());
	TGraphErrors *Combined_fit_alpha_eta = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());
	TGraphErrors *Combined_fit_n_eta = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());
	TGraphErrors *Combined_fit_STD_eta = new TGraphErrors(HI_dwidth_raw->GetN() + pp_dwidth_raw->GetN());

	TGraphErrors *Combined_hist_dm_raw = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_hist_dwidth_raw = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());

	TGraphErrors *Combined_hist_dm_eta = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_hist_dwidth_eta = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());

	TGraphErrors *Combined_chi2_dm_raw = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_chi2_dwidth_raw = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_chi2_dm_eta = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());
	TGraphErrors *Combined_chi2_dwidth_eta = new TGraphErrors(pp_dm_raw->GetN() + HI_dm_raw->GetN());

	for (int i = 0; i < pp_dm_raw->GetN(); ++i)
	{
		Combined_fit_dm_raw->SetPoint(i, pp_dm_raw->GetX()[i], pp_dm_raw->GetY()[i]);
		Combined_fit_dwidth_raw->SetPoint(i, pp_dwidth_raw->GetX()[i], pp_dwidth_raw->GetY()[i]);
		Combined_fit_alpha_raw->SetPoint(i, pp_alpha_raw->GetX()[i], pp_alpha_raw->GetY()[i]);
		Combined_fit_n_raw->SetPoint(i, pp_n_raw->GetX()[i], pp_n_raw->GetY()[i]);
		Combined_fit_STD_raw->SetPoint(i, pp_STD_raw->GetX()[i], pp_STD_raw->GetY()[i]);

		Combined_fit_dm_raw->SetPointError(i, pp_dm_raw->GetErrorX(i), pp_dm_raw->GetErrorY(i));
		Combined_fit_dwidth_raw->SetPointError(i, pp_dwidth_raw->GetErrorX(i), pp_dwidth_raw->GetErrorY(i));
		Combined_fit_alpha_raw->SetPointError(i, pp_alpha_raw->GetErrorX(i), pp_alpha_raw->GetErrorY(i));
		Combined_fit_n_raw->SetPointError(i, pp_n_raw->GetErrorX(i), pp_n_raw->GetErrorY(i));
		Combined_fit_STD_raw->SetPointError(i, pp_STD_raw->GetErrorX(i), pp_STD_raw->GetErrorY(i));

		Combined_fit_dm_eta->SetPoint(i, pp_dm_raw->GetX()[i], pp_dm_eta->GetY()[i]);
		Combined_fit_dwidth_eta->SetPoint(i, pp_dwidth_eta->GetX()[i], pp_dwidth_eta->GetY()[i]);
		Combined_fit_alpha_eta->SetPoint(i, pp_alpha_eta->GetX()[i], pp_alpha_eta->GetY()[i]);
		Combined_fit_n_eta->SetPoint(i, pp_n_eta->GetX()[i], pp_n_eta->GetY()[i]);
		Combined_fit_STD_eta->SetPoint(i, pp_STD_eta->GetX()[i], pp_STD_eta->GetY()[i]);

		Combined_fit_dm_eta->SetPointError(i, pp_dm_raw->GetErrorX(i), pp_dm_eta->GetErrorY(i));
		Combined_fit_dwidth_eta->SetPointError(i, pp_dwidth_eta->GetErrorX(i), pp_dwidth_eta->GetErrorY(i));
		Combined_fit_alpha_eta->SetPointError(i, pp_alpha_eta->GetErrorX(i), pp_alpha_eta->GetErrorY(i));
		Combined_fit_n_eta->SetPointError(i, pp_n_eta->GetErrorX(i), pp_n_eta->GetErrorY(i));
		Combined_fit_STD_eta->SetPointError(i, pp_STD_eta->GetErrorX(i), pp_STD_eta->GetErrorY(i));

		Combined_hist_dm_raw->SetPoint(i, pp_hist_dm_raw->GetX()[i], pp_hist_dm_raw->GetY()[i]);
		Combined_hist_dwidth_raw->SetPoint(i, pp_hist_dw_raw->GetX()[i], pp_hist_dw_raw->GetY()[i]);

		Combined_hist_dm_raw->SetPointError(i, pp_hist_dm_raw->GetErrorX(i), pp_hist_dm_raw->GetErrorY(i));
		Combined_hist_dwidth_raw->SetPointError(i, pp_hist_dw_raw->GetErrorX(i), pp_hist_dw_raw->GetErrorY(i));

		Combined_hist_dm_eta->SetPoint(i, pp_hist_dm_eta->GetX()[i], pp_hist_dm_eta->GetY()[i]);
		Combined_hist_dwidth_eta->SetPoint(i, pp_hist_dw_eta->GetX()[i], pp_hist_dw_eta->GetY()[i]);

		Combined_hist_dm_eta->SetPointError(i, pp_hist_dm_eta->GetErrorX(i), pp_hist_dm_eta->GetErrorY(i));
		Combined_hist_dwidth_eta->SetPointError(i, pp_hist_dw_eta->GetErrorX(i), pp_hist_dw_eta->GetErrorY(i));

		Combined_chi2_dm_raw->SetPoint(i, pp_chi2_dm_raw->GetX()[i], pp_chi2_dm_raw->GetY()[i]);
		Combined_chi2_dwidth_raw->SetPoint(i, pp_chi2_dw_raw->GetX()[i], pp_chi2_dw_raw->GetY()[i]);

		Combined_chi2_dm_raw->SetPointError(i, pp_chi2_dm_raw->GetErrorX(i), pp_chi2_dm_raw->GetErrorY(i));
		Combined_chi2_dwidth_raw->SetPointError(i, pp_chi2_dw_raw->GetErrorX(i), pp_chi2_dw_raw->GetErrorY(i));

		Combined_chi2_dm_eta->SetPoint(i, pp_chi2_dm_eta->GetX()[i], pp_chi2_dm_eta->GetY()[i]);
		Combined_chi2_dwidth_eta->SetPoint(i, pp_chi2_dw_eta->GetX()[i], pp_chi2_dw_eta->GetY()[i]);

		Combined_chi2_dm_eta->SetPointError(i, pp_chi2_dm_eta->GetErrorX(i), pp_chi2_dm_eta->GetErrorY(i));
		Combined_chi2_dwidth_eta->SetPointError(i, pp_chi2_dw_eta->GetErrorX(i), pp_chi2_dw_eta->GetErrorY(i));
	}
	for (int i = 0; i < HI_dm_raw->GetN(); ++i)
	{

		Combined_fit_dm_raw->SetPoint(i + pp_dm_raw->GetN(), HI_dm_raw->GetX()[i] + 22, HI_dm_raw->GetY()[i]);
		Combined_fit_dwidth_raw->SetPoint(i + pp_dm_raw->GetN(), HI_dwidth_raw->GetX()[i] + 22, HI_dwidth_raw->GetY()[i]);
		Combined_fit_alpha_raw->SetPoint(i + pp_dm_raw->GetN(), HI_alpha_raw->GetX()[i] + 22, HI_alpha_raw->GetY()[i]);
		Combined_fit_n_raw->SetPoint(i + pp_dm_raw->GetN(), HI_n_raw->GetX()[i] + 22, HI_n_raw->GetY()[i]);
		Combined_fit_STD_raw->SetPoint(i + pp_dm_raw->GetN(), HI_STD_raw->GetX()[i] + 22, HI_STD_raw->GetY()[i]);

		Combined_fit_dm_raw->SetPointError(i + pp_dm_raw->GetN(), HI_dm_raw->GetErrorX(i), HI_dm_raw->GetErrorY(i));
		Combined_fit_dwidth_raw->SetPointError(i + pp_dm_raw->GetN(), HI_dwidth_raw->GetErrorX(i), HI_dwidth_raw->GetErrorY(i));
		Combined_fit_alpha_raw->SetPointError(i + pp_dm_raw->GetN(), HI_alpha_raw->GetErrorX(i), HI_alpha_raw->GetErrorY(i));
		Combined_fit_n_raw->SetPointError(i + pp_dm_raw->GetN(), HI_n_raw->GetErrorX(i), HI_n_raw->GetErrorY(i));
		Combined_fit_STD_raw->SetPointError(i + pp_dm_raw->GetN(), HI_STD_raw->GetErrorX(i), HI_STD_raw->GetErrorY(i));

		Combined_fit_dm_eta->SetPoint(i + pp_dm_raw->GetN(), HI_dm_eta->GetX()[i] + 22, HI_dm_eta->GetY()[i]);
		Combined_fit_dwidth_eta->SetPoint(i + pp_dm_raw->GetN(), HI_dwidth_eta->GetX()[i] + 22, HI_dwidth_eta->GetY()[i]);
		Combined_fit_alpha_eta->SetPoint(i + pp_dm_raw->GetN(), HI_alpha_eta->GetX()[i] + 22, HI_alpha_eta->GetY()[i]);
		Combined_fit_n_eta->SetPoint(i + pp_dm_raw->GetN(), HI_n_eta->GetX()[i] + 22, HI_n_eta->GetY()[i]);
		Combined_fit_STD_eta->SetPoint(i + pp_dm_raw->GetN(), HI_STD_eta->GetX()[i] + 22, HI_STD_eta->GetY()[i]);

		Combined_fit_dm_eta->SetPointError(i + pp_dm_raw->GetN(), HI_dm_eta->GetErrorX(i), HI_dm_eta->GetErrorY(i));
		Combined_fit_dwidth_eta->SetPointError(i + pp_dm_raw->GetN(), HI_dwidth_eta->GetErrorX(i), HI_dwidth_eta->GetErrorY(i));
		Combined_fit_alpha_eta->SetPointError(i + pp_dm_raw->GetN(), HI_alpha_eta->GetErrorX(i), HI_alpha_eta->GetErrorY(i));
		Combined_fit_n_eta->SetPointError(i + pp_dm_raw->GetN(), HI_n_eta->GetErrorX(i), HI_n_eta->GetErrorY(i));
		Combined_fit_STD_eta->SetPointError(i + pp_dm_raw->GetN(), HI_STD_eta->GetErrorX(i), HI_STD_eta->GetErrorY(i));

		Combined_hist_dm_raw->SetPoint(i + pp_dm_raw->GetN(), HI_hist_dm_raw->GetX()[i] + 22, HI_hist_dm_raw->GetY()[i]);
		Combined_hist_dwidth_raw->SetPoint(i + pp_dm_raw->GetN(), HI_hist_dw_raw->GetX()[i] + 22, HI_hist_dw_raw->GetY()[i]);

		Combined_hist_dm_raw->SetPointError(i + pp_dm_raw->GetN(), HI_hist_dm_raw->GetErrorX(i), HI_hist_dm_raw->GetErrorY(i));
		Combined_hist_dwidth_raw->SetPointError(i + pp_dm_raw->GetN(), HI_hist_dw_raw->GetErrorX(i), HI_hist_dw_raw->GetErrorY(i));

		Combined_hist_dm_eta->SetPoint(i + pp_dm_raw->GetN(), HI_hist_dm_eta->GetX()[i] + 22, HI_hist_dm_eta->GetY()[i]);
		Combined_hist_dwidth_eta->SetPoint(i + pp_dm_raw->GetN(), HI_hist_dw_eta->GetX()[i] + 22, HI_hist_dw_eta->GetY()[i]);

		Combined_hist_dm_eta->SetPointError(i + pp_dm_raw->GetN(), HI_hist_dm_eta->GetErrorX(i), HI_hist_dm_eta->GetErrorY(i));
		Combined_hist_dwidth_eta->SetPointError(i + pp_dm_raw->GetN(), HI_hist_dw_eta->GetErrorX(i), HI_hist_dw_eta->GetErrorY(i));

		Combined_chi2_dm_raw->SetPoint(i + pp_dm_raw->GetN(), HI_chi2_dm_raw->GetX()[i] + 22, HI_chi2_dm_raw->GetY()[i]);
		Combined_chi2_dwidth_raw->SetPoint(i + pp_dm_raw->GetN(), HI_chi2_dw_raw->GetX()[i] + 22, HI_chi2_dw_raw->GetY()[i]);

		Combined_chi2_dm_raw->SetPointError(i + pp_dm_raw->GetN(), HI_chi2_dm_raw->GetErrorX(i), HI_chi2_dm_raw->GetErrorY(i));
		Combined_chi2_dwidth_raw->SetPointError(i + pp_dm_raw->GetN(), HI_chi2_dw_raw->GetErrorX(i), HI_chi2_dw_raw->GetErrorY(i));

		Combined_chi2_dm_eta->SetPoint(i + pp_dm_raw->GetN(), HI_chi2_dm_eta->GetX()[i] + 22, HI_chi2_dm_eta->GetY()[i]);
		Combined_chi2_dwidth_eta->SetPoint(i + pp_dm_raw->GetN(), HI_chi2_dw_eta->GetX()[i] + 22, HI_chi2_dw_eta->GetY()[i]);

		Combined_chi2_dm_eta->SetPointError(i + pp_dm_raw->GetN(), HI_chi2_dm_eta->GetErrorX(i), HI_chi2_dm_eta->GetErrorY(i));
		Combined_chi2_dwidth_eta->SetPointError(i + pp_dm_raw->GetN(), HI_chi2_dw_eta->GetErrorX(i), HI_chi2_dw_eta->GetErrorY(i));
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
	gStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
	gStyle->SetPadTickY(1);
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "XYZ");
	gStyle->SetLabelSize(0.03, "XYZ");

	TCanvas *c_dm_raw = new TCanvas("c_dm_raw", "", 1600, 800);
	TCanvas *c_dwidth_raw = new TCanvas("c_dwidth_raw", "", 1600, 800);
	TCanvas *c_alpha_raw = new TCanvas("c_alpha_raw", "", 1600, 800);
	TCanvas *c_n_raw = new TCanvas("c_n_raw", "", 1600, 800);
	TCanvas *c_STD_raw = new TCanvas("c_STD_raw", "", 1600, 800);
	TCanvas *c_diff_dm_raw = new TCanvas("c_diff_dm_raw", "", 800, 800);
	TCanvas *c_diff_dwidth_raw = new TCanvas("c_diff_dwidth_raw", "", 800, 800);

	TCanvas *c_dm_eta = new TCanvas("c_dm_eta", "", 1600, 800);
	TCanvas *c_dwidth_eta = new TCanvas("c_dwidth_eta", "", 1600, 800);
	TCanvas *c_alpha_eta = new TCanvas("c_alpha_eta", "", 1600, 800);
	TCanvas *c_n_eta = new TCanvas("c_n_eta", "", 1600, 800);
	TCanvas *c_STD_eta = new TCanvas("c_STD_eta", "", 1600, 800);
	TCanvas *c_diff_dm_eta = new TCanvas("c_diff_dm_eta", "", 800, 800);
	TCanvas *c_diff_dwidth_eta = new TCanvas("c_diff_dwidth_eta", "", 800, 800);

	changecosmetic(c_dm_raw, Combined_fit_dm_raw, 1, type, Combined_hist_dm_raw, Combined_chi2_dm_raw, Diff_dM_array_raw, Diff_dM_array_err_raw);
	changecosmetic(c_dwidth_raw, Combined_fit_dwidth_raw, 2, type, Combined_hist_dwidth_raw, Combined_chi2_dwidth_raw, Diff_dWidth_array_raw, Diff_dWidth_array_err_raw);
	changecosmetic(c_alpha_raw, Combined_fit_alpha_raw, 3, type);
	changecosmetic(c_n_raw, Combined_fit_n_raw, 4, type);
	changecosmetic(c_STD_raw, Combined_fit_STD_raw, 5, type);

	changecosmetic(c_dm_eta, Combined_fit_dm_eta, 1, type, Combined_hist_dm_eta, Combined_chi2_dm_eta, Diff_dM_array_eta, Diff_dM_array_err_eta);
	changecosmetic(c_dwidth_eta, Combined_fit_dwidth_eta, 2, type, Combined_hist_dwidth_eta, Combined_chi2_dwidth_eta, Diff_dWidth_array_eta, Diff_dWidth_array_err_eta);
	changecosmetic(c_alpha_eta, Combined_fit_alpha_eta, 3, type);
	changecosmetic(c_n_eta, Combined_fit_n_eta, 4, type);
	changecosmetic(c_STD_eta, Combined_fit_STD_eta, 5, type);

	getdiffcosmetic(c_diff_dm_raw, Diff_dM_array_raw, Diff_dM_array_err_raw, 1);
	getdiffcosmetic(c_diff_dwidth_raw, Diff_dWidth_array_raw, Diff_dWidth_array_err_raw, 0);
	getdiffcosmetic(c_diff_dm_eta, Diff_dM_array_eta, Diff_dM_array_err_eta, 1);
	getdiffcosmetic(c_diff_dwidth_eta, Diff_dWidth_array_eta, Diff_dWidth_array_err_eta, 0);

	if (type == 1)
	{
		c_dm_raw->SaveAs("./combined/raw/Combined_dm_raw.png");
		c_dwidth_raw->SaveAs("./combined/raw/Combined_dwidth_raw.png");
		c_alpha_raw->SaveAs("./combined/raw/Combined_alpha_raw.png");
		c_n_raw->SaveAs("./combined/raw/Combined_n_raw.png");
		c_STD_raw->SaveAs("./combined/raw/Combined_STD_raw.png");
		c_diff_dm_raw->SaveAs("./diff/raw/dm.png");
		c_diff_dwidth_raw->SaveAs("./diff/raw/dw.png");

		c_dm_eta->SaveAs("./combined/eta/Combined_dm_eta.png");
		c_dwidth_eta->SaveAs("./combined/eta/Combined_dwidth_eta.png");
		c_alpha_eta->SaveAs("./combined/eta/Combined_alpha_eta.png");
		c_n_eta->SaveAs("./combined/eta/Combined_n_eta.png");
		c_STD_eta->SaveAs("./combined/eta/Combined_STD_eta.png");
		c_diff_dm_eta->SaveAs("./diff/eta/dm.png");
		c_diff_dwidth_eta->SaveAs("./diff/eta/dw.png");
	}
	else
	{
		c_dm_raw->SaveAs("./combined/raw/Combined_dm_exp_raw.png");
		c_dwidth_raw->SaveAs("./combined/raw/Combined_dwidth_exp_raw.png");
		c_alpha_raw->SaveAs("./combined/raw/Combined_alpha_exp_raw.png");
		c_n_raw->SaveAs("./combined/raw/Combined_n_exp_raw.png");
		c_STD_raw->SaveAs("./combined/raw/Combined_STD_exp_raw.png");
		c_diff_dm_raw->SaveAs("./diff/raw/dm_exp.png");
		c_diff_dwidth_raw->SaveAs("./diff/raw/dw_exp.png");

		c_dm_eta->SaveAs("./combined/eta/Combined_dm_exp_eta.png");
		c_dwidth_eta->SaveAs("./combined/eta/Combined_dwidth_exp_eta.png");
		c_alpha_eta->SaveAs("./combined/eta/Combined_alpha_exp_eta.png");
		c_n_eta->SaveAs("./combined/eta/Combined_n_exp_eta.png");
		c_STD_eta->SaveAs("./combined/eta/Combined_STD_exp_eta.png");
		c_diff_dm_eta->SaveAs("./diff/eta/dm_exp.png");
		c_diff_dwidth_eta->SaveAs("./diff/eta/dw_exp.png");
	}
}