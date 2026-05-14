
void AnnotatePoints_ValStatSyst(
    TCanvas *c,
    TGraphErrors *g,
    const std::vector<TBox *> &systBoxes,
    TH1 *frame,
    double xOffset = 0.08,     // in "bin units" (since your x = 1..N)
    double yOffsetFrac = 0.00, // fraction of y-range (0 = same height; try 0.03 if overlaps)
    double textSize = 0.030,
    int textColor = kBlack,
    int align = 22,   // 12 = left-middle, 22 = center-middle, 32 = right-middle
    int precision = 3 // digits after decimal
)
{
    if (!c || !g || !frame)
        return;

    c->cd();

    const double yMin = frame->GetMinimum();
    const double yMax = frame->GetMaximum();
    const double ySpan = (yMax - yMin > 0) ? (yMax - yMin) : 1.0;

    // If systBoxes is not same size, we’ll try to match by x-bin proximity
    auto FindSystFromBoxes = [&](double xbin, double yval) -> double
    {
        if (systBoxes.empty())
            return 0.0;

        // 1) assume same ordering (common case)
        // If sizes match, use index mapping by nearest integer bin
        // but safer: find closest box center in x
        double best = 0.0;
        double bestDx = 1e9;

        for (auto b : systBoxes)
        {
            if (!b)
                continue;
            double xc = 0.5 * (b->GetX1() + b->GetX2());
            double dx = std::abs(xc - xbin);
            if (dx < bestDx)
            {
                bestDx = dx;
                best = 0.5 * std::abs(b->GetY2() - b->GetY1()); // half-height = syst
            }
        }
        return best;
    };

    TLatex latex;
    latex.SetNDC(false);
    latex.SetTextFont(42);
    latex.SetTextSize(textSize);
    latex.SetTextColor(textColor);
    latex.SetTextAlign(align);

    const int n = g->GetN();
    for (int i = 0; i < n; ++i)
    {
        double x, y;
        g->GetPoint(i, x, y);
        double estat = g->GetErrorY(i);
        double esyst = FindSystFromBoxes(x, y);

        // place text a bit to the right of the marker
        double xt = x + xOffset;
        double yt = y + yOffsetFrac * ySpan;

        TString fmt = Form("%%.%df #pm %%.%df #pm %%.%df", precision, precision, precision);
        TString s = Form(fmt.Data(), y, estat, esyst);

        latex.DrawLatex(xt, yt, s);
    }

    c->Modified();
    c->Update();
}

struct SubsetPlotOut
{
    TCanvas *c = nullptr;
    TH1F *frame = nullptr;
    TGraphErrors *g = nullptr;
    std::vector<TBox *> boxes;
    TLegend *leg = nullptr;
};

SubsetPlotOut DrawSubsetPlot(
    const char *cname,
    const char *ctitle,
    const char *ytitle,
    const std::vector<int> &idx,        // which original bins, in desired order
    const std::vector<TString> &labels, // same size as idx
    TEllipse *statEll[11],
    TEllipse *systEll[11],
    bool useX, // true: use ellipse X (dM); false: use ellipse Y (dGamma)
    double yMin, double yMax,
    double halfWidth = 0.35,
    double magicScale = std::sqrt(2.3))
{
    SubsetPlotOut out;
    const int n = (int)idx.size();
    if (n <= 0)
        return out;

    // frame (x bins are 1..n)
    out.frame = new TH1F(Form("frame_%s", cname), Form(";Centrality;%s", ytitle), n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i)
        out.frame->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    out.frame->SetMinimum(yMin);
    out.frame->SetMaximum(yMax);

    std::vector<double> bx(n), y(n), ey(n);
    out.boxes.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        const int j = idx[i];
        bx[i] = i + 1;

        if (!statEll[j] || !systEll[j])
        {
            y[i] = 0;
            ey[i] = 0;
            continue;
        }

        // --- your ellipse -> (center, projected sigmas)
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

            const double a2 = a * a, b2 = b * b;
            const double sigxx = c * c * a2 + s * s * b2;
            const double sigyy = s * s * a2 + c * c * b2;
            sig_x = std::sqrt(std::max(0.0, sigxx));
            sig_y = std::sqrt(std::max(0.0, sigyy));
        };

        double xc_s, yc_s, sx_s, sy_s;
        double xc_sy, yc_sy, sx_sy, sy_sy;
        GetXYErrorsFromEllipse(statEll[j], xc_s, yc_s, sx_s, sy_s);
        GetXYErrorsFromEllipse(systEll[j], xc_sy, yc_sy, sx_sy, sy_sy);

        // your “magical number”
        sx_s /= magicScale;
        sy_s /= magicScale;
        sx_sy /= magicScale;
        sy_sy /= magicScale;

        // pick X or Y dimension
        const double center = useX ? xc_s : yc_s;
        const double statSig = useX ? sx_s : sy_s;
        const double systSig = useX ? sx_sy : sy_sy;

        y[i] = center;
        ey[i] = statSig;

        // syst box (vertical, centered at y[i], half-height = systSig)
        auto b = new TBox(bx[i] - halfWidth, y[i] - systSig, bx[i] + halfWidth, y[i] + systSig);
        b->SetFillColorAlpha(kGray + 1, 0.25);
        b->SetLineColor(kGray + 2);
        out.boxes.push_back(b);
    }

    // graph (stat errors)
    out.g = new TGraphErrors(n, bx.data(), y.data(), nullptr, ey.data());
    out.g->SetMarkerStyle(useX ? 20 : 21);
    out.g->SetMarkerSize(1.0);
    out.g->SetLineWidth(2);
    out.g->SetLineColor(kBlack);

    // draw
    out.c = new TCanvas(cname, ctitle, 800, 800);
    out.frame->Draw();
    for (auto b : out.boxes)
        b->Draw("same");
    out.g->Draw("PZ SAME");

    // legend
    out.leg = new TLegend(0.15, 0.75, 0.50, 0.89);
    out.leg->SetFillStyle(0);
    out.leg->SetBorderSize(0);
    out.leg->SetTextFont(42);
    out.leg->SetTextSize(0.035);
    out.leg->AddEntry(out.g, Form("Center #pm stat (projected #sigma_{%s})", useX ? "x" : "y"), "lep");
    if (!out.boxes.empty())
        out.leg->AddEntry(out.boxes.front(), Form("Syst: #pm #sigma_{%s}^{syst} (box)", useX ? "x" : "y"), "f");
    out.leg->Draw();

    out.c->Modified();
    out.c->Update();
    return out;
}

SubsetPlotOut DrawSubsetPlot_PbPbPlusPP(
    const char *cname,
    const char *ctitle,
    const char *ytitle,
    const std::vector<int> &idx,        // PbPb bin indices; use -1 to mean "pp"
    const std::vector<TString> &labels, // same size as idx
    TEllipse *pbpb_stat[11],
    TEllipse *pbpb_syst[11],
    TEllipse *pp_stat, // pp ellipse (single)
    TEllipse *pp_syst, // pp ellipse (single)
    bool useX,         // true: X (mass); false: Y (width)
    double yMin, double yMax,
    double halfWidth = 0.35,
    double magicScale = 1.515,
    double addOffset = 0.0 // use 91.1876 for mass, 2.4955 for width if you want absolute values
)
{
    SubsetPlotOut out;
    const int n = (int)idx.size();
    if (n <= 0)
        return out;

    out.frame = new TH1F(Form("frame_%s", cname), Form(";Centrality;%s", ytitle), n, 0.5, n + 0.5);
    for (int i = 0; i < n; ++i)
        out.frame->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    out.frame->SetMinimum(yMin);
    out.frame->SetMaximum(yMax);

    std::vector<double> bx(n), y(n), ey(n);
    out.boxes.reserve(n);

    auto GetXYErrorsFromEllipse = [](TEllipse *e, double &xc, double &yc, double &sig_x, double &sig_y)
    {
        xc = yc = sig_x = sig_y = 0.0;
        if (!e)
            return;

        xc = e->GetX1();
        yc = e->GetY1();
        const double a = e->GetR1();
        const double b = e->GetR2();
        const double th = e->GetTheta() * TMath::DegToRad();
        const double c = std::cos(th), s = std::sin(th);

        const double a2 = a * a, b2 = b * b;
        const double sigxx = c * c * a2 + s * s * b2;
        const double sigyy = s * s * a2 + c * c * b2;
        sig_x = std::sqrt(std::max(0.0, sigxx));
        sig_y = std::sqrt(std::max(0.0, sigyy));
    };

    for (int i = 0; i < n; ++i)
    {
        bx[i] = i + 1;

        TEllipse *e_stat = nullptr;
        TEllipse *e_syst = nullptr;

        if (idx[i] >= 0)
        {
            // PbPb bin
            const int j = idx[i];
            e_stat = pbpb_stat[j];
            e_syst = pbpb_syst[j];
        }
        else
        {
            // pp
            e_stat = pp_stat;
            e_syst = pp_syst;
        }

        if (!e_stat || !e_syst)
        {
            y[i] = 0;
            ey[i] = 0;
            continue;
        }

        double xc_s, yc_s, sx_s, sy_s;
        double xc_sy, yc_sy, sx_sy, sy_sy;
        GetXYErrorsFromEllipse(e_stat, xc_s, yc_s, sx_s, sy_s);
        GetXYErrorsFromEllipse(e_syst, xc_sy, yc_sy, sx_sy, sy_sy);

        sx_s /= magicScale;
        sy_s /= magicScale;
        sx_sy /= magicScale;
        sy_sy /= magicScale;

        const double center = (useX ? xc_s : yc_s) + addOffset;
        const double statSig = (useX ? sx_s : sy_s);
        const double systSig = (useX ? sx_sy : sy_sy);

        y[i] = center;
        ey[i] = statSig;

        auto b = new TBox(bx[i] - halfWidth, y[i] - systSig, bx[i] + halfWidth, y[i] + systSig);
        b->SetFillColorAlpha(kGray + 1, 0.25);
        b->SetLineColor(kGray + 2);
        out.boxes.push_back(b);
    }

    out.g = new TGraphErrors(n, bx.data(), y.data(), nullptr, ey.data());
    out.g->SetMarkerStyle(useX ? 20 : 21);
    out.g->SetMarkerSize(1.0);
    out.g->SetLineWidth(2);
    out.g->SetLineColor(kBlack);

    out.c = new TCanvas(cname, ctitle, 800, 600);
    out.frame->Draw();
    for (auto b : out.boxes)
        b->Draw("same");
    out.g->Draw("PZ SAME");

    out.leg = new TLegend(0.15, 0.75, 0.50, 0.89);
    out.leg->SetFillStyle(0);
    out.leg->SetBorderSize(0);
    out.leg->SetTextFont(42);
    out.leg->SetTextSize(0.035);
    out.leg->AddEntry(out.g, Form("Center #pm stat (projected #sigma_{%s})", useX ? "x" : "y"), "lep");
    if (!out.boxes.empty())
        out.leg->AddEntry(out.boxes.front(), Form("Syst: #pm #sigma_{%s}^{syst} (box)", useX ? "x" : "y"), "f");
    out.leg->Draw();

    out.c->Modified();
    out.c->Update();
    return out;
}

void FitConstAndTestZero(TGraphErrors *g, int i0, int i1, const char *tag)
{
    if (!g || i1 < i0)
        return;

    const int N = i1 - i0 + 1;

    // ---- Extract x, y, ey
    std::vector<double> x(N), y(N), ey(N);
    for (int k = 0; k < N; ++k)
    {
        g->GetPoint(i0 + k, x[k], y[k]);
        ey[k] = g->GetErrorY(i0 + k);
    }

    // ---- Determine bin width automatically
    double binw = 0.0;
    if (N >= 2)
        binw = x[1] - x[0];
    else
        binw = 1.0; // fallback

    // (optional safety: average spacing)
    if (N >= 3)
    {
        double sum = 0;
        for (int k = 1; k < N; ++k)
            sum += (x[k] - x[k - 1]);
        binw = sum / (N - 1);
    }

    // ---- Fit range
    const double xmin = x.front() - 0.5 * binw;
    const double xmax = x.back() + 0.5 * binw;

    // ---- Constant fit
    TF1 fconst(Form("fconst_%s", tag), "[0]", xmin, xmax);
    fconst.SetParameter(0, 0.0);

    g->Fit(&fconst, "Q0R");

    const double mu = fconst.GetParameter(0);
    const double mu_err = fconst.GetParError(0);
    const double chi2fit = fconst.GetChisquare();
    const int ndffit = fconst.GetNDF();
    const double pfit = TMath::Prob(chi2fit, ndffit);

    // ---- Test vs zero
    double chi2_zero = 0.0;
    int used = 0;
    for (int k = 0; k < N; ++k)
    {
        if (ey[k] <= 0)
            continue;
        chi2_zero += (y[k] * y[k]) / (ey[k] * ey[k]);
        used++;
    }

    const int ndf_zero = used - 1;
    const double p_zero = TMath::Prob(chi2_zero, ndf_zero);

    // ---- Output
    std::cout << "==== " << tag << " (points " << i0 << "-" << i1 << ") ====\n";
    std::cout << "Auto bin width = " << binw << "\n";
    std::cout << "Fit range      = [" << xmin << ", " << xmax << "]\n";
    std::cout << "Const fit: mu = " << mu << " ± " << mu_err
              << "   chi2/ndf = " << chi2fit << "/" << ndffit
              << "   p = " << pfit << "\n";
    std::cout << "Test mu = 0: chi2/ndf = " << chi2_zero << "/" << ndf_zero
              << "   p = " << p_zero << "\n";
}

void get_final_plots()
{

    gStyle->SetOptStat(0);

    TFile *f1 = new TFile("./contourtest/degen/contour.root", "READ");

    TEllipse *PbPb_stat[11];
    TEllipse *pp_stat[11];
    TEllipse *PbPb_sub_pp_stat[11];

    TEllipse *PbPb_syst[11];
    TEllipse *pp_syst[11];
    TEllipse *PbPb_sub_pp_syst[11];

    const int n_cent = 11;

    TGraphErrors *g_PbPb;
    TGraphErrors *g_pp;
    TGraphErrors *g_PbPb_sub_pp;

    TString cent_array[11] = {"0_10", "10_20", "20_30", "30_100", "0_100", "0_6.5", "6.5_14.5", "14.5_26.5", "26.5_100", "0_15", "15_100"};

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

    FitConstAndTestZero(gX, 0, 3, "dM_first_4");
    FitConstAndTestZero(gY, 0, 3, "dW_first_4");
    FitConstAndTestZero(gX, 5, 8, "dM_second_4");
    FitConstAndTestZero(gY, 5, 8, "dW_second_4");

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

    // vertical separators
    TLine *l1 = new TLine(5.5, frameX->GetMinimum(), 5.5, frameX->GetMaximum());
    TLine *l2 = new TLine(9.5, frameX->GetMinimum(), 9.5, frameX->GetMaximum());

    l1->SetLineStyle(2); // dashed
    l2->SetLineStyle(2);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);

    l1->Draw("same");
    l2->Draw("same");

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

    TLine *l1y = new TLine(5.5, frameY->GetMinimum(), 5.5, frameY->GetMaximum());
    TLine *l2y = new TLine(9.5, frameY->GetMinimum(), 9.5, frameY->GetMaximum());

    l1y->SetLineStyle(2);
    l2y->SetLineStyle(2);
    l1y->SetLineWidth(2);
    l2y->SetLineWidth(2);

    l1y->Draw("same");
    l2y->Draw("same");

    gPad->Modified();
    gPad->Update();

    // ---- NEW: two subsets for PbPb_sub_pp
    std::vector<int> idx_first5 = {0, 1, 2, 3, 4};
    std::vector<TString> lab_first5 = {"0_10", "10_20", "20_30", "30_100", "0_100"};

    std::vector<int> idx_65scheme = {5, 6, 7, 8, 4}; // move 0_100 (index 4) to last
    std::vector<TString> lab_65scheme = {"0_6.5", "6.5_14.5", "14.5_26.5", "26.5_100", "0_100"};

    // dM (X of ellipse)
    auto dM_first5 = DrawSubsetPlot(
        "c_dM_first5", "dM: first 5 bins", "#DeltaM (PbPb - pp)",
        idx_first5, lab_first5,
        PbPb_sub_pp_stat, PbPb_sub_pp_syst,
        true, -0.8, 0.8);

    auto dM_65 = DrawSubsetPlot(
        "c_dM_65scheme", "dM: 6.5 scheme (0_100 last)", "#DeltaM (PbPb - pp)",
        idx_65scheme, lab_65scheme,
        PbPb_sub_pp_stat, PbPb_sub_pp_syst,
        true, -0.8, 0.8);

    // dGamma (Y of ellipse)
    auto dW_first5 = DrawSubsetPlot(
        "c_dW_first5", "dGamma: first 5 bins", "#Delta#Gamma (PbPb - pp)",
        idx_first5, lab_first5,
        PbPb_sub_pp_stat, PbPb_sub_pp_syst,
        false, -0.8, 0.8);

    auto dW_65 = DrawSubsetPlot(
        "c_dW_65scheme", "dGamma: 6.5 scheme (0_100 last)", "#Delta#Gamma (PbPb - pp)",
        idx_65scheme, lab_65scheme,
        PbPb_sub_pp_stat, PbPb_sub_pp_syst,
        false, -0.8, 0.8);

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

    auto DrawSep = [](double x, TH1 *frame)
    {
        TLine *l = new TLine(x, frame->GetMinimum(), x, frame->GetMaximum());
        l->SetLineStyle(7); // dashed
        l->SetLineWidth(2);
        l->SetLineColor(kGray + 2);
        l->Draw("same");
    };

    DrawSep(5.5, frameX2);
    DrawSep(9.5, frameX2);
    DrawSep(11.5, frameX2);

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
    DrawSep(5.5, frameY2);
    DrawSep(9.5, frameY2);
    DrawSep(11.5, frameY2);
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

    // pp ellipses you already loaded: pp_stat[0], pp_syst[0]
    TEllipse *ppStat = pp_stat[0];
    TEllipse *ppSyst = pp_syst[0];

    // --- (A) first 5 PbPb bins + pp at last
    std::vector<int> idx_pbpb_first5_pp = {0, 1, 2, 3, 4, -1};
    std::vector<TString> lab_pbpb_first5_pp = {"0_10", "10_20", "20_30", "30_100", "0_100", "pp"};

    // --- (B) 6.5 scheme, move 0_100 to last among PbPb, then pp last
    std::vector<int> idx_pbpb_65_pp = {5, 6, 7, 8, 4, -1};
    std::vector<TString> lab_pbpb_65_pp = {"0_6.5", "6.5_14.5", "14.5_26.5", "26.5_100", "0_100", "pp"};

    // If you want ABSOLUTE values like your screenshot:
    // mass offset = 91.1876, width offset = 2.4955
    const double mass0 = 91.1876;
    const double width0 = 2.4955;

    // dM (absolute)
    auto dM_pbpbpp_first5 = DrawSubsetPlot_PbPbPlusPP(
        "c_dM_pbpbpp_first5", "PbPb first 5 + pp (last)", "#DeltaM",
        idx_pbpb_first5_pp, lab_pbpb_first5_pp,
        PbPb_stat, PbPb_syst, ppStat, ppSyst,
        true, 89.5, 92.5,
        0.35, 1.515, mass0);

    auto dM_pbpbpp_65 = DrawSubsetPlot_PbPbPlusPP(
        "c_dM_pbpbpp_65", "PbPb 6.5 scheme + pp (last)", "#DeltaM",
        idx_pbpb_65_pp, lab_pbpb_65_pp,
        PbPb_stat, PbPb_syst, ppStat, ppSyst,
        true, 89.5, 92.5,
        0.35, 1.515, mass0);

    // dGamma (absolute)
    auto dW_pbpbpp_first5 = DrawSubsetPlot_PbPbPlusPP(
        "c_dW_pbpbpp_first5", "PbPb first 5 + pp (last)", "#Delta#Gamma",
        idx_pbpb_first5_pp, lab_pbpb_first5_pp,
        PbPb_stat, PbPb_syst, ppStat, ppSyst,
        false, 0.5, 5.5,
        0.35, 1.515, width0);

    auto dW_pbpbpp_65 = DrawSubsetPlot_PbPbPlusPP(
        "c_dW_pbpbpp_65", "PbPb 6.5 scheme + pp (last)", "#Delta#Gamma",
        idx_pbpb_65_pp, lab_pbpb_65_pp,
        PbPb_stat, PbPb_syst, ppStat, ppSyst,
        false, 0.5, 5.5,
        0.35, 1.515, width0);

    AnnotatePoints_ValStatSyst(cX, gX, boxesX, frameX, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(cY, gY, boxesY, frameY, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(cX2, gX2, boxesX2, frameX2, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(cY2, gY2, boxesY2, frameY2, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dM_first5.c, dM_first5.g, dM_first5.boxes, dM_first5.frame, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dW_first5.c, dW_first5.g, dW_first5.boxes, dW_first5.frame, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dM_65.c, dM_65.g, dM_65.boxes, dM_65.frame, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dW_65.c, dW_65.g, dW_65.boxes, dW_65.frame, 0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dM_pbpbpp_first5.c, dM_pbpbpp_first5.g,
                               dM_pbpbpp_first5.boxes, dM_pbpbpp_first5.frame,
                               0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dW_pbpbpp_first5.c, dW_pbpbpp_first5.g,
                               dW_pbpbpp_first5.boxes, dW_pbpbpp_first5.frame,
                               0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dM_pbpbpp_65.c, dM_pbpbpp_65.g,
                               dM_pbpbpp_65.boxes, dM_pbpbpp_65.frame,
                               0.0, 0.02, 0.017);
    AnnotatePoints_ValStatSyst(dW_pbpbpp_65.c, dW_pbpbpp_65.g,
                               dW_pbpbpp_65.boxes, dW_pbpbpp_65.frame,
                               0.0, 0.02, 0.017);

    TFile *final_plot_root_file = new TFile("./send_to_frank.root", "UPDATE");
    final_plot_root_file->cd();
    cX->Write("dM_PbPb_sub_pp", 2);
    cY->Write("dW_PbPb_sub_pp", 2);
    cX2->Write("dM_PbPb_and_pp", 2);
    cY2->Write("dW_PbPb_and_pp", 2);

    final_plot_root_file->cd();
    dM_first5.c->Write("dM_PbPb_sub_pp_first5", 2);
    dW_first5.c->Write("dW_PbPb_sub_pp_first5", 2);
    dM_65.c->Write("dM_PbPb_sub_pp_65scheme", 2);
    dW_65.c->Write("dW_PbPb_sub_pp_65scheme", 2);

    final_plot_root_file->cd();
    dM_pbpbpp_first5.c->Write("dM_PbPb_first5_plus_pp", 2);
    dW_pbpbpp_first5.c->Write("dW_PbPb_first5_plus_pp", 2);
    dM_pbpbpp_65.c->Write("dM_PbPb_65scheme_plus_pp", 2);
    dW_pbpbpp_65.c->Write("dW_PbPb_65scheme_plus_pp", 2);
}