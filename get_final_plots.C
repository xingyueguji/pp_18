void get_final_plots()
{
    TFile *f1 = new TFile("./contourtest/degen/contour.root", "READ");

    TEllipse *PbPb_stat[5];
    TEllipse *pp_stat[5];
    TEllipse *PbPb_sub_pp_stat[5];

    TEllipse *PbPb_syst[5];
    TEllipse *pp_syst[5];
    TEllipse *PbPb_sub_pp_syst[5];

    const int n_cent = 5;

    TGraphErrors *g_PbPb;
    TGraphErrors *g_pp;
    TGraphErrors *g_PbPb_sub_pp;

    TString cent_array[5] = {"0_10", "10_20", "20_30", "30_100", "0_100"};

    for (int i = 0; i < n_cent; i++)
    {
        f1->cd();
        PbPb_stat[i] = (TEllipse *)f1->Get("sig_1_" + cent_array[i] + "_stat_PbPb");
        pp_stat[i] = (TEllipse *)f1->Get("sig_1_" + cent_array[i] + "_stat_pp");
        PbPb_sub_pp_stat[i] = (TEllipse *)f1->Get("sig_1_" + cent_array[i] + "_stat_PbPb_sub_pp");

        PbPb_syst[i] = (TEllipse *)f1->Get("sig_1_" + cent_array[i] + "_syst_PbPb");
        pp_syst[i] = (TEllipse *)f1->Get("sig_1_" + cent_array[i] + "_syst_pp");
        PbPb_sub_pp_syst[i] = (TEllipse *)f1->Get("sig_1_" + cent_array[i] + "_syst_PbPb_sub_pp");
    }

    auto GetXYErrorsFromEllipse = [](TEllipse *e, double &xc, double &yc, double &sig_x, double &sig_y)
    {
        xc = yc = sig_x = sig_y = 0.0;
        if (!e)
            return;

        xc = e->GetX1();
        yc = e->GetY1();
        const double a = e->GetR1();
        const double b = e->GetR2();
        const double th = e->GetTheta() * TMath::DegToRad(); // radians
        const double c = std::cos(th), s = std::sin(th);

        // Sigma = R diag(a^2, b^2) R^T  =>  sig_x^2 = Sigma_xx, sig_y^2 = Sigma_yy
        const double a2 = a * a, b2 = b * b;
        const double sigxx = c * c * a2 + s * s * b2;
        const double sigyy = s * s * a2 + c * c * b2;
        sig_x = std::sqrt(std::max(0.0, sigxx));
        sig_y = std::sqrt(std::max(0.0, sigyy));
    };

    // ---- 2) Build data arrays
    std::vector<double> binX(n_cent), xval(n_cent), yval(n_cent), ex_stat(n_cent), ey_stat(n_cent);
    std::vector<TBox *> boxesX;
    boxesX.reserve(n_cent);
    std::vector<TBox *> boxesY;
    boxesY.reserve(n_cent);

    // choose a visual half-width for the syst boxes along the bin axis
    const double halfWidth = 0.35;

    double yminX = +1e30, ymaxX = -1e30;
    double yminY = +1e30, ymaxY = -1e30;

    for (int i = 0; i < n_cent; ++i)
    {
        // Skip if missing
        if (!PbPb_sub_pp_stat[i])
            continue;
        if (!PbPb_sub_pp_syst[i])
            continue;

        double xc_s, yc_s, sx_s, sy_s;
        double xc_sy, yc_sy, sx_sy, sy_sy;

        GetXYErrorsFromEllipse(PbPb_sub_pp_stat[i], xc_s, yc_s, sx_s, sy_s);
        GetXYErrorsFromEllipse(PbPb_sub_pp_syst[i], xc_sy, yc_sy, sx_sy, sy_sy);

        // Add this magical number for both stat and syst error
        sx_s = sx_s / 1.515;
        sy_s = sy_s / 1.515;
        sx_sy = sx_sy / 1.515;
        sy_sy = sy_sy / 1.515;

        // Positions along x = bin index 1..n_cent
        binX[i] = i + 1;

        // Points (centers) from STAT ellipse
        xval[i] = xc_s;
        yval[i] = yc_s;

        ex_stat[i] = sx_s; // X-graph uses σx(stat)
        ey_stat[i] = sy_s; // Y-graph uses σy(stat)

        // SYST boxes around the point (axis-aligned): ±σ_syst
        // For X graph: vertical box centered at xval[i] with half-height = sx_sy
        auto bX = new TBox(binX[i] - halfWidth, xval[i] - sx_sy, binX[i] + halfWidth, xval[i] + sx_sy);
        bX->SetFillColorAlpha(kGray + 1, 0.25);
        bX->SetLineColor(kGray + 2);
        boxesX.push_back(bX);

        // For Y graph: vertical box centered at yval[i] with half-height = sy_sy
        auto bY = new TBox(binX[i] - halfWidth, yval[i] - sy_sy, binX[i] + halfWidth, yval[i] + sy_sy);
        bY->SetFillColorAlpha(kGray + 1, 0.25);
        bY->SetLineColor(kGray + 2);
        boxesY.push_back(bY);

        // Track min/max to set frame ranges (include syst envelopes)
        yminX = std::min(yminX, xval[i] - sx_sy);
        ymaxX = std::max(ymaxX, xval[i] + sx_sy);
        yminY = std::min(yminY, yval[i] - sy_sy);
        ymaxY = std::max(ymaxY, yval[i] + sy_sy);
    }

    // Pad a bit
    auto padRange = [](double &mn, double &mx)
    {
        double span = mx - mn;
        if (span <= 0)
        {
            mn -= 1.0;
            mx += 1.0;
        } // fallback
        else
        {
            mn -= 0.10 * span;
            mx += 0.10 * span;
        }
    };
    padRange(yminX, ymaxX);
    padRange(yminY, ymaxY);

    // ---- 3) Frames with bin labels
    TH1F *frameX = new TH1F("dM", ";Centrality;#DeltaM (PbPb - pp)", n_cent, 0.5, n_cent + 0.5);
    TH1F *frameY = new TH1F("dWidth", ";Centrality;#Delta#Gamma (PbPb - pp)", n_cent, 0.5, n_cent + 0.5);
    for (int i = 0; i < n_cent; ++i)
    {
        frameX->GetXaxis()->SetBinLabel(i + 1, cent_array[i]); // cent_array[i] must be TString or const char*
        frameY->GetXaxis()->SetBinLabel(i + 1, cent_array[i]);
    }
    frameX->SetMinimum(-0.8);
    frameX->SetMaximum(0.8);
    frameY->SetMinimum(-0.8);
    frameY->SetMaximum(0.8);

    // ---- 4) Graphs with STAT errors
    TGraphErrors *gX = new TGraphErrors(n_cent, binX.data(), xval.data(), nullptr, ex_stat.data());
    gX->SetName("g_dx_stat");
    gX->SetMarkerStyle(20);
    gX->SetMarkerSize(1.0);
    gX->SetLineWidth(2);
    gX->SetLineColor(kBlack);

    TGraphErrors *gY = new TGraphErrors(n_cent, binX.data(), yval.data(), nullptr, ey_stat.data());
    gY->SetName("g_dy_stat");
    gY->SetMarkerStyle(21);
    gY->SetMarkerSize(1.0);
    gY->SetLineWidth(2);
    gY->SetLineColor(kBlack);

    // ---- 5) Draw X canvas
    TCanvas *cX = new TCanvas("cX", "Delta X vs centrality", 800, 600);
    frameX->Draw();
    for (auto b : boxesX)
        b->Draw("same"); // SYST boxes
    gX->Draw("PZ SAME"); // STAT error bars on points

    // Legend
    auto legX = new TLegend(0.15, 0.75, 0.48, 0.89);
    legX->SetFillStyle(0);
    legX->SetBorderSize(0);
    legX->SetTextFont(42);
    legX->SetTextSize(0.035);
    legX->AddEntry(gX, "Center #pm stat (projected #sigma_{x})", "lep");
    if (!boxesX.empty())
        legX->AddEntry(boxesX.front(), "Syst: #pm #sigma_{x}^{syst} (box)", "f");
    legX->Draw();

    // ---- 6) Draw Y canvas
    TCanvas *cY = new TCanvas("cY", "Delta Y vs centrality", 800, 600);
    frameY->Draw();
    for (auto b : boxesY)
        b->Draw("same");
    gY->Draw("PZ SAME");

    auto legY = new TLegend(0.15, 0.75, 0.48, 0.89);
    legY->SetFillStyle(0);
    legY->SetBorderSize(0);
    legY->SetTextFont(42);
    legY->SetTextSize(0.035);
    legY->AddEntry(gY, "Center #pm stat (projected #sigma_{y})", "lep");
    if (!boxesY.empty())
        legY->AddEntry(boxesY.front(), "Syst: #pm #sigma_{y}^{syst} (box)", "f");
    legY->Draw();

    gPad->Modified();
    gPad->Update();

    // ==================== EXTRA CANVASES: PbPb-only + last bin = pp ====================
    const int n_all = n_cent + 1; // 6 bins total (5 PbPb + 1 pp)

    std::vector<double> binX2(n_all), xval2(n_all), yval2(n_all), ex_stat2(n_all), ey_stat2(n_all);
    std::vector<TBox *> boxesX2;
    boxesX2.reserve(n_all);
    std::vector<TBox *> boxesY2;
    boxesY2.reserve(n_all);

    // labels: copy cent_array and append "pp"
    std::vector<TString> labels(n_all);
    for (int i = 0; i < n_cent; ++i)
        labels[i] = cent_array[i];
    labels[n_all - 1] = "pp";

    double yminX2 = +1e30, ymaxX2 = -1e30;
    double yminY2 = +1e30, ymaxY2 = -1e30;

    for (int i = 0; i < n_all; ++i)
    {
        binX2[i] = i + 1;

        TEllipse *e_stat = nullptr;
        TEllipse *e_syst = nullptr;

        if (i < n_cent)
        {
            // PbPb bins
            e_stat = PbPb_stat[i];
            e_syst = PbPb_syst[i];
        }
        else
        {
            // last bin = pp
            e_stat = pp_stat[0]; // choose the pp contour you want; often pp has only one bin
            e_syst = pp_syst[0];
        }

        if (!e_stat || !e_syst)
        {
            // leave defaults (zeros); or continue to skip drawing this bin
            continue;
        }

        double xc_s, yc_s, sx_s, sy_s;
        double xc_sy, yc_sy, sx_sy, sy_sy;

        GetXYErrorsFromEllipse(e_stat, xc_s, yc_s, sx_s, sy_s);
        GetXYErrorsFromEllipse(e_syst, xc_sy, yc_sy, sx_sy, sy_sy);

        // Add this magical number for both stat and syst error
        sx_s = sx_s / 1.515;
        sy_s = sy_s / 1.515;
        sx_sy = sx_sy / 1.515;
        sy_sy = sy_sy / 1.515;

        xval2[i] = xc_s + 91.1876;
        yval2[i] = yc_s + 2.4955;
        ex_stat2[i] = sx_s;
        ey_stat2[i] = sy_s;

        // syst boxes (vertical, ±syst)
        auto bX = new TBox(binX2[i] - halfWidth, xval2[i] - sx_sy, binX2[i] + halfWidth, xval2[i] + sx_sy);
        bX->SetFillColorAlpha(kGray + 1, 0.25);
        bX->SetLineColor(kGray + 2);
        boxesX2.push_back(bX);

        auto bY = new TBox(binX2[i] - halfWidth, yval2[i] - sy_sy, binX2[i] + halfWidth, yval2[i] + sy_sy);
        bY->SetFillColorAlpha(kGray + 1, 0.25);
        bY->SetLineColor(kGray + 2);
        boxesY2.push_back(bY);

        yminX2 = std::min(yminX2, xval2[i] - sx_sy);
        ymaxX2 = std::max(ymaxX2, xval2[i] + sx_sy);
        yminY2 = std::min(yminY2, yval2[i] - sy_sy);
        ymaxY2 = std::max(ymaxY2, yval2[i] + sy_sy);
    }

    // pad ranges
    padRange(yminX2, ymaxX2);
    padRange(yminY2, ymaxY2);

    // frames with labels (6 bins)
    TH1F *frameX2 = new TH1F("dM_pbpbpp", ";Centrality;#DeltaM", n_all, 0.5, n_all + 0.5);
    TH1F *frameY2 = new TH1F("dW_pbpbpp", ";Centrality;#Delta#Gamma", n_all, 0.5, n_all + 0.5);
    for (int i = 0; i < n_all; ++i)
    {
        frameX2->GetXaxis()->SetBinLabel(i + 1, labels[i]);
        frameY2->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    }
    frameX2->SetMinimum(89.5);
    frameX2->SetMaximum(92.5);
    frameY2->SetMinimum(0.5);
    frameY2->SetMaximum(5.5);

    // graphs
    TGraphErrors *gX2 = new TGraphErrors(n_all, binX2.data(), xval2.data(), nullptr, ex_stat2.data());
    gX2->SetName("g_dx_pbpbpp");
    gX2->SetMarkerStyle(20);
    gX2->SetMarkerSize(1.0);
    gX2->SetLineWidth(2);
    gX2->SetLineColor(kBlack);

    TGraphErrors *gY2 = new TGraphErrors(n_all, binX2.data(), yval2.data(), nullptr, ey_stat2.data());
    gY2->SetName("g_dy_pbpbpp");
    gY2->SetMarkerStyle(21);
    gY2->SetMarkerSize(1.0);
    gY2->SetLineWidth(2);
    gY2->SetLineColor(kBlack);

    // draw X canvas
    TCanvas *cX2 = new TCanvas("cX2", "Delta X: PbPb (5) + pp (last)", 800, 600);
    frameX2->Draw();
    for (auto b : boxesX2)
        b->Draw("same");
    gX2->Draw("PZ SAME");
    auto legX2 = new TLegend(0.15, 0.75, 0.50, 0.89);
    legX2->SetFillStyle(0);
    legX2->SetBorderSize(0);
    legX2->SetTextFont(42);
    legX2->SetTextSize(0.035);
    legX2->AddEntry(gX2, "Center #pm stat (projected #sigma_{x})", "lep");
    if (!boxesX2.empty())
        legX2->AddEntry(boxesX2.front(), "Syst: #pm #sigma_{x}^{syst} (box)", "f");
    legX2->Draw();

    // draw Y canvas
    TCanvas *cY2 = new TCanvas("cY2", "Delta Y: PbPb (5) + pp (last)", 800, 600);
    frameY2->Draw();
    for (auto b : boxesY2)
        b->Draw("same");
    gY2->Draw("PZ SAME");
    auto legY2 = new TLegend(0.15, 0.75, 0.50, 0.89);
    legY2->SetFillStyle(0);
    legY2->SetBorderSize(0);
    legY2->SetTextFont(42);
    legY2->SetTextSize(0.035);
    legY2->AddEntry(gY2, "Center #pm stat (projected #sigma_{y})", "lep");
    if (!boxesY2.empty())
        legY2->AddEntry(boxesY2.front(), "Syst: #pm #sigma_{y}^{syst} (box)", "f");
    legY2->Draw();

    gPad->Modified();
    gPad->Update();

    TFile *final_plot_root_file = new TFile("./send_to_frank.root", "UPDATE");
    final_plot_root_file->cd();
    cX->Write("dM_PbPb_sub_pp", 2);
    cY->Write("dW_PbPb_sub_pp", 2);
    cX2->Write("dM_PbPb_and_pp", 2);
    cY2->Write("dW_PbPb_and_pp", 2);
}