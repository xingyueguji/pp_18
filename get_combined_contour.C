TVectorD fit_ellipse(TGraph *g)
{
    TVectorD ellipse;

    if (!g)
        return ellipse; // just a precaution
    if (g->GetN() < 6)
        return ellipse; // just a precaution

    Int_t i;
    Double_t tmp;

    Int_t N = g->GetN();
    Double_t xmin, xmax, ymin, ymax, X0, Y0;
    g->ComputeRange(xmin, ymin, xmax, ymax);
#if 1 /* 0 or 1 */
    X0 = (xmax + xmin) / 2.0;
    Y0 = (ymax + ymin) / 2.0;
#else  /* 0 or 1 */
    X0 = Y0 = 0.0;
#endif /* 0 or 1 */

    TMatrixD D1(N, 3); // quadratic part of the design matrix
    TMatrixD D2(N, 3); // linear part of the design matrix

    for (i = 0; i < N; i++)
    {
        Double_t x = (g->GetX())[i] - X0;
        Double_t y = (g->GetY())[i] - Y0;
        D1[i][0] = x * x;
        D1[i][1] = x * y;
        D1[i][2] = y * y;
        D2[i][0] = x;
        D2[i][1] = y;
        D2[i][2] = 1.0;
    }

    // quadratic part of the scatter matrix
    TMatrixD S1(TMatrixD::kAtA, D1);
    // combined part of the scatter matrix
    TMatrixD S2(D1, TMatrixD::kTransposeMult, D2);
    // linear part of the scatter matrix
    TMatrixD S3(TMatrixD::kAtA, D2);
    S3.Invert(&tmp);
    S3 *= -1.0;
    if (tmp == 0.0)
    {
        std::cout << "fit_ellipse : linear part of the scatter matrix is singular!" << std::endl;
        return ellipse;
    }
    // for getting a2 from a1
    TMatrixD T(S3, TMatrixD::kMultTranspose, S2);
    // reduced scatter matrix
    TMatrixD M(S2, TMatrixD::kMult, T);
    M += S1;
    // premultiply by inv(C1)
    for (i = 0; i < 3; i++)
    {
        tmp = M[0][i] / 2.0;
        M[0][i] = M[2][i] / 2.0;
        M[2][i] = tmp;
        M[1][i] *= -1.0;
    }
    // solve eigensystem
    TMatrixDEigen eig(M); // note: eigenvectors are not normalized
    const TMatrixD &evec = eig.GetEigenVectors();
    // const TVectorD &eval = eig.GetEigenValuesRe();
    if ((eig.GetEigenValuesIm()).Norm2Sqr() != 0.0)
    {
        std::cout << "fit_ellipse : eigenvalues have nonzero imaginary parts!" << std::endl;
        return ellipse;
    }
    // evaluate a’Ca (in order to find the eigenvector for min. pos. eigenvalue)
    for (i = 0; i < 3; i++)
    {
        tmp = 4.0 * evec[0][i] * evec[2][i] - evec[1][i] * evec[1][i];
        if (tmp > 0.0)
            break;
    }
    if (i > 2)
    {
        std::cout << "fit_ellipse : no min. pos. eigenvalue found!" << std::endl;
        // i = 2;
        return ellipse;
    }
    // eigenvector for min. pos. eigenvalue
    TVectorD a1(TMatrixDColumn_const(evec, i));
    tmp = a1.Norm2Sqr();
    if (tmp > 0.0)
    {
        a1 *= 1.0 / std::sqrt(tmp); // normalize this eigenvector
    }
    else
    {
        std::cout << "fit_ellipse : eigenvector for min. pos. eigenvalue is NULL!" << std::endl;
        return ellipse;
    }
    TVectorD a2(T * a1);

    // ellipse coefficients
    ellipse.ResizeTo(8);
    ellipse[0] = X0;    // "X0"
    ellipse[1] = Y0;    // "Y0"
    ellipse[2] = a1[0]; // "A"
    ellipse[3] = a1[1]; // "B"
    ellipse[4] = a1[2]; // "C"
    ellipse[5] = a2[0]; // "D"
    ellipse[6] = a2[1]; // "E"
    ellipse[7] = a2[2]; // "F"

    return ellipse;
}

//
// http://mathworld.wolfram.com/Ellipse.html
// http://mathworld.wolfram.com/QuadraticCurve.html
// http://mathworld.wolfram.com/ConicSection.html
//
// "Using the Ellipse to Fit and Enclose Data Points"
// Charles F. Van Loan
// http://www.cs.cornell.edu/cv/OtherPdf/Ellipse.pdf
//
// input parameter is a reference to a "TVectorD" which describes
// an ellipse according to the equation:
//   0 = A * (x - X0)^2 + B * (x - X0) * (y - Y0) + C * (y - Y0)^2
//     + D * (x - X0) + E * (y - Y0) + F
// conic[0] = "X0"
// conic[1] = "Y0"
// conic[2] = "A"
// conic[3] = "B"
// conic[4] = "C"
// conic[5] = "D"
// conic[6] = "E"
// conic[7] = "F"
//
// returns a "TVectorD" ("empty" in case any problems encountered):
// ellipse[0] = ellipse's "x" center ("x0")
// ellipse[1] = ellipse's "y" center ("y0")
// ellipse[2] = ellipse's "semimajor" axis along "x" ("a" > 0)
// ellipse[3] = ellipse's "semiminor" axis along "y" ("b" > 0)
// ellipse[4] = ellipse's axes rotation angle ("theta" = -45 ... 135 degrees)
//
TVectorD ConicToParametric(const TVectorD &conic)
{
    TVectorD ellipse;

    if (conic.GetNrows() != 8)
    {
        std::cout << "ConicToParametric : improper input vector length!" << std::endl;
        return ellipse;
    }

    Double_t a, b, theta;
    Double_t x0 = conic[0]; // = X0
    Double_t y0 = conic[1]; // = Y0

    // http://mathworld.wolfram.com/Ellipse.html
    Double_t A = conic[2];
    Double_t B = conic[3] / 2.0;
    Double_t C = conic[4];
    Double_t D = conic[5] / 2.0;
    Double_t F = conic[6] / 2.0;
    Double_t G = conic[7];

    Double_t J = B * B - A * C;
    Double_t Delta = A * F * F + C * D * D + J * G - 2.0 * B * D * F;
    Double_t I = -(A + C);

    // http://mathworld.wolfram.com/QuadraticCurve.html
    if (!((Delta != 0.0) && (J < 0.0) && (I != 0.0) && (Delta / I < 0.0)))
    {
        std::cout << "ConicToParametric : ellipse (real) specific constraints not met!" << std::endl;
        return ellipse;
    }

    x0 += (C * D - B * F) / J;
    y0 += (A * F - B * D) / J;

    Double_t tmp = std::sqrt((A - C) * (A - C) + 4.0 * B * B);
    a = std::sqrt(2.0 * Delta / J / (I + tmp));
    b = std::sqrt(2.0 * Delta / J / (I - tmp));

    theta = 0.0;
    if (B != 0.0)
    {
        tmp = (A - C) / 2.0 / B;
        theta = -45.0 * (std::atan(tmp) / TMath::PiOver2());
        if (tmp < 0.0)
        {
            theta -= 45.0;
        }
        else
        {
            theta += 45.0;
        }
        if (A > C)
            theta += 90.0;
    }
    else if (A > C)
        theta = 90.0;

    // try to keep "a" > "b"
    if (a < b)
    {
        tmp = a;
        a = b;
        b = tmp;
        theta -= 90.0;
    }
    // try to keep "theta" = -45 ... 135 degrees
    if (theta < -45.0)
        theta += 180.0;
    if (theta > 135.0)
        theta -= 180.0;

    // ellipse coefficients
    ellipse.ResizeTo(5);
    ellipse[0] = x0;    // ellipse's "x" center
    ellipse[1] = y0;    // ellipse's "y" center
    ellipse[2] = a;     // ellipse's "semimajor" axis along "x"
    ellipse[3] = b;     // ellipse's "semiminor" axis along "y"
    ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)

    cout << "x is " << ellipse[0] << " y is " << ellipse[1] << " semimajor x is " << ellipse[2] << " semimajor y is " << ellipse[3] << " angle is " << ellipse[4] << endl;

    return ellipse;
}

void savePbPb(const TVectorD &d1, const TVectorD &d2, TCanvas *c1, int cent)
{
    c1->cd();
    TEllipse *e = new TEllipse(d1[0], d1[1], // "x0", "y0"
                               d1[2], d1[3], // "a", "b"
                               0, 360,
                               d1[4]); // "theta" (in degrees)
    e->SetFillStyle(0);                // hollow
    e->SetLineColor(kGreen);
    e->SetLineWidth(5);
    e->Draw("SAME");

    TEllipse *e1 = new TEllipse(d2[0], d2[1], // "x0", "y0"
                                d2[2], d2[3], // "a", "b"
                                0, 360,
                                d2[4]); // "theta" (in degrees)
    e1->SetFillStyle(0);                // hollow
    e1->SetLineColor(kOrange);
    e1->SetLineWidth(5);
    e1->Draw("SAME");

    c1->SaveAs(Form("./contourtest/PbPb_%i.png", cent));
}
TVectorD doeverything(TString canvasname, TString savename, TFile *f1, TVector2 *local_min, int opt = 0)
{
    TCanvas *c1;
    c1 = (TCanvas *)f1->Get(canvasname);
    TList *list = c1->GetListOfPrimitives();
    TGraph *g1 = (TGraph *)list->At(1);
    TGraph *g2 = (TGraph *)list->At(2);
    TVectorD coinc_1 = fit_ellipse(g1);
    TVectorD conic_2 = fit_ellipse(g2);
    TVectorD ellipse_1 = ConicToParametric(coinc_1);
    TVectorD ellipse_2 = ConicToParametric(conic_2);
    c1->cd();
    TEllipse *e = new TEllipse(ellipse_1[0], ellipse_1[1], // "x0", "y0"
                               ellipse_1[2], ellipse_1[3], // "a", "b"
                               0, 360,
                               ellipse_1[4]); // "theta" (in degrees)
    e->SetFillStyle(0);                       // hollow
    e->SetLineColor(kBlue);
    e->SetLineWidth(5);
    e->Draw("SAME");

    TEllipse *e1 = new TEllipse(ellipse_2[0], ellipse_2[1], // "x0", "y0"
                                ellipse_2[2], ellipse_2[3], // "a", "b"
                                0, 360,
                                ellipse_2[4]); // "theta" (in degrees)
    e1->SetFillStyle(0);                       // hollow
    e1->SetLineColor(kOrange);
    e1->SetLineWidth(5);
    e1->Draw("SAME");

    double x = local_min->X();
    double y = local_min->Y();

    TGraph *g = new TGraph();
    g->SetPoint(0, x, y);
    g->SetMarkerStyle(20); // Solid circle
    g->SetMarkerSize(1.5);
    g->SetMarkerColor(kRed);
    g->Draw("P SAME");

    TLatex latex1;
    latex1.SetTextSize(0.02);
    latex1.SetTextAlign(12); // left-aligned, vertically centered

    for (int i = 0; i < g->GetN(); ++i)
    {
        double x, y;
        g->GetPoint(i, x, y);
        TString label = Form("%.3f, %.3f", x, y);
        latex1.DrawLatex(x + 0.01, y + 0.1, label);
    }

    TGraph *g_1 = new TGraph();
    g_1->SetPoint(0, ellipse_1[0], ellipse_1[1]);
    g_1->SetMarkerStyle(20); // Solid circle
    g_1->SetMarkerSize(1.5);
    g_1->SetMarkerColor(kBlue);
    g_1->Draw("P SAME");

    TGraph *g_2 = new TGraph();
    g_2->SetPoint(0, ellipse_2[0], ellipse_2[1]);
    g_2->SetMarkerStyle(20); // Solid circle
    g_2->SetMarkerSize(1.5);
    g_2->SetMarkerColor(kOrange);
    g_2->Draw("P SAME");

    /*double theta_deg = 84.143;

    double theta_rad = theta_deg * TMath::Pi() / 180.0;

    // Major axis direction (theta)
    double dx_major = 0.0947848 * cos(theta_rad);
    double dy_major = 0.0947848 * sin(theta_rad);

    // Minor axis direction (theta + 90 degrees)
    double dx_minor = 0.0547605 * cos(theta_rad + TMath::Pi()/2);
    double dy_minor = 0.0547605 * sin(theta_rad + TMath::Pi()/2);

    double x0 = -0.135428;
    double y0 = 0.0776946;

    TLine *major_axis = new TLine(x0 - dx_major, y0 - dy_major, x0 + dx_major, y0 + dy_major);
    TLine *minor_axis = new TLine(x0 - dx_minor, y0 - dy_minor, x0 + dx_minor, y0 + dy_minor);

    major_axis->Draw("SAME");
    minor_axis->Draw("SAME");*/

    c1->SaveAs(savename);

    if (opt == 1)
    {
        return ellipse_1;
    }
    if (opt == 2)
    {
        return ellipse_2;
    }
    else
    {
        return TVectorD();
    }
}
void Drawsubresult(TCanvas *c2, TEllipse *sig_1, TEllipse *sig_2, TString savename)
{
    c2->cd();

    // Create a dummy frame
    TH2F *frame = new TH2F("frame", "", 1, -0.3, 0.3, 1, -0.4, 0.4);
    frame->SetStats(0);
    frame->Draw();

    sig_1->SetFillStyle(0);
    sig_1->SetLineColor(kBlue);
    sig_1->SetLineWidth(5);
    sig_1->Draw("SAME");

    sig_2->SetFillStyle(0);
    sig_2->SetLineColor(kOrange);
    sig_2->SetLineWidth(5);
    sig_2->Draw("SAME");

    TGraph *g_1 = new TGraph();
    g_1->SetPoint(0, sig_1->GetX1(), sig_1->GetY1());
    g_1->SetMarkerStyle(20); // Solid circle
    g_1->SetMarkerSize(1.5);
    g_1->SetMarkerColor(kBlue);
    g_1->Draw("P SAME TEXT");

    TGraph *g_2 = new TGraph();
    g_2->SetPoint(0, 0, 0);
    g_2->SetMarkerStyle(20); // Solid circle
    g_2->SetMarkerSize(1.5);
    g_2->SetMarkerColor(kRed);
    g_2->Draw("P SAME TEXT");

    c2->SaveAs(savename);
    c2->Clear();
    delete frame;
    delete g_1;
    delete g_2;
}

TEllipse *CombineEllipsesFromVectors(const TVectorD &v1, const TVectorD &v2)
{
    auto GetCovariance = [](const TVectorD &v)
    {
        double a = v[2];
        double b = v[3];
        double theta = v[4] * TMath::DegToRad();

        double cosT = std::cos(theta);
        double sinT = std::sin(theta);

        TMatrixD D(2, 2);
        D(0, 0) = a * a;
        D(1, 1) = b * b;
        D(0, 1) = D(1, 0) = 0;

        TMatrixD R(2, 2);
        R(0, 0) = cosT;
        R(0, 1) = -sinT;
        R(1, 0) = sinT;
        R(1, 1) = cosT;

        TMatrixD Rt(TMatrixD::kTransposed, R);
        TMatrixD Sigma = R * D * Rt;

        return Sigma;
    };

    TMatrixD cov1 = GetCovariance(v1);
    TMatrixD cov2 = GetCovariance(v2);
    TMatrixD covSum = cov1;
    covSum += cov2;

    double epsilon = 1e-12;
    for (int i = 0; i < covSum.GetNrows(); ++i)
    {
        for (int j = 0; j < covSum.GetNcols(); ++j)
        {
            if (std::abs(covSum(i, j)) < epsilon)
            {
                covSum(i, j) = 0.0;
            }
        }
    }

    std::cout << "covSum matrix elements:" << std::endl;
    for (int i = 0; i < covSum.GetNrows(); ++i)
    {
        for (int j = 0; j < covSum.GetNcols(); ++j)
        {
            std::cout << covSum(i, j) << "\t";
        }
        std::cout << std::endl;
    }

    TMatrixDEigen eigen(covSum);
    TVectorD evals = eigen.GetEigenValuesRe();
    TMatrixD evecs = eigen.GetEigenVectors();

    double a_comb = std::sqrt(std::max(evals[0], evals[1]));
    double b_comb = std::sqrt(std::min(evals[0], evals[1]));

    int idx = (evals[0] > evals[1]) ? 0 : 1;
    double vx = evecs(0, idx);
    double vy = evecs(1, idx);

    double theta_comb = std::atan2(vy, vx) * TMath::RadToDeg();
    theta_comb = std::atan2(vy, vx) * TMath::RadToDeg();
    if (theta_comb < 0)
        theta_comb += 180.0;

    double xc = v1[0] - v2[0];
    double yc = v1[1] - v2[1];

    TEllipse *result = new TEllipse(xc, yc, a_comb, b_comb, 0, 360, theta_comb);
    result->SetLineColor(kRed);
    result->SetLineWidth(2);

    return result;
}

void DrawEllipsesComparison(const TVectorD &v1, const TVectorD &v2, TString Savename)
{
    auto CreateAxisLine = [](double angle_deg, double xmin, double xmax, double ymin, double ymax)
    {
        double angle_rad = angle_deg * TMath::DegToRad();
        double dx = std::cos(angle_rad);
        double dy = std::sin(angle_rad);
        double length = std::max(xmax - xmin, ymax - ymin) * 2.0;

        auto line = new TLine(-length * dx, -length * dy, length * dx, length * dy);
        line->SetLineStyle(2); // dashed
        line->SetLineColor(kGray + 2);
        line->SetLineWidth(1);
        return line;
    };

    auto CreatePerpendicularLine = [](double x0, double y0, double angle_deg)
    {
        double angle_rad = angle_deg * TMath::DegToRad();
        double cos_theta = std::cos(angle_rad);
        double sin_theta = std::sin(angle_rad);
        double proj = x0 * cos_theta + y0 * sin_theta;
        double x_proj = proj * cos_theta;
        double y_proj = proj * sin_theta;
        return new TLine(x0, y0, x_proj, y_proj);
    };

    // Reuse the previous function
    TEllipse *e1 = new TEllipse(v1[0], v1[1], v1[2], v1[3], 0, 360, v1[4]);
    TEllipse *e2 = new TEllipse(v2[0], v2[1], v2[2], v2[3], 0, 360, v2[4]);
    TEllipse *e3 = CombineEllipsesFromVectors(v1, v2);

    double x0 = e3->GetX1();
    double y0 = e3->GetY1();
    double theta = e3->GetTheta();

    // Cosmetics
    e1->SetLineColor(kBlack);
    e1->SetLineWidth(2);
    e1->SetLineStyle(1);
    e1->SetFillStyle(0);

    e2->SetLineColor(kBlue);
    e2->SetLineWidth(2);
    e2->SetLineStyle(2);
    e2->SetFillStyle(0);

    e3->SetLineColor(kRed);
    e3->SetLineWidth(2);
    e3->SetLineStyle(1);
    e3->SetFillStyle(0);

    // Find frame range
    // Find initial ranges
    double x_min = std::min({v1[0] - v1[2], v2[0] - v2[2], e3->GetX1() - e3->GetR1()});
    double x_max = std::max({v1[0] + v1[2], v2[0] + v2[2], e3->GetX1() + e3->GetR1()});
    double y_min = std::min({v1[1] - v1[3], v2[1] - v2[3], e3->GetY1() - e3->GetR2()});
    double y_max = std::max({v1[1] + v1[3], v2[1] + v2[3], e3->GetY1() + e3->GetR2()});

    // Force (0,0) to be included
    x_min = std::min(x_min, -0.1);
    x_max = std::max(x_max, 0.1);
    y_min = std::min(y_min, -0.1);
    y_max = std::max(y_max, 0.1);

    // Expand by 20%
    double x_center = 0.5 * (x_min + x_max);
    double y_center = 0.5 * (y_min + y_max);
    double x_half = 0.5 * (x_max - x_min);
    double y_half = 0.5 * (y_max - y_min);

    x_half *= 1.2;
    y_half *= 1.2;

    // Use the larger half-range for both x and y
    double range_half = std::max(x_half, y_half);

    x_min = x_center - range_half;
    x_max = x_center + range_half;
    y_min = y_center - range_half;
    y_max = y_center + range_half;

    TCanvas *c1 = new TCanvas("c1", "Ellipse Comparison", 800, 800);
    c1->DrawFrame(x_min, y_min, x_max, y_max, "Contours;Mass;Width");

    e1->Draw("same");
    e2->Draw("same");
    e3->Draw("same");

    TLine *line_a = CreateAxisLine(theta, x_min, x_max, y_min, y_max);
    TLine *line_b = CreateAxisLine(theta + 90.0, x_min, x_max, y_min, y_max);
    line_a->Draw("same");
    line_b->Draw("same");

    auto perp_a = CreatePerpendicularLine(x0, y0, theta);
    auto perp_b = CreatePerpendicularLine(x0, y0, theta + 90.0);
    perp_a->SetLineColor(kRed + 2);
    perp_a->SetLineStyle(3);
    perp_b->SetLineColor(kRed + 2);
    perp_b->SetLineStyle(3);
    perp_a->Draw("same");
    perp_b->Draw("same");

    TMarker *m1 = new TMarker(v1[0], v1[1], 20); // 20 = filled circle marker style
    m1->SetMarkerColor(kBlack);
    m1->SetMarkerSize(1.2);
    m1->Draw("same");

    TMarker *m2 = new TMarker(v2[0], v2[1], 20);
    m2->SetMarkerColor(kBlue);
    m2->SetMarkerSize(1.2);
    m2->Draw("same");

    TMarker *m3 = new TMarker(e3->GetX1(), e3->GetY1(), 20);
    m3->SetMarkerColor(kRed);
    m3->SetMarkerSize(1.2);
    m3->Draw("same");

    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->AddEntry(e1, "PbPb", "l");
    leg->AddEntry(e2, "pp", "l");
    leg->AddEntry(e3, "PbPb - pp", "l");
    leg->Draw();

    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(13); // left align, top align

    double ystart = 0.9;
    double dy = 0.06;

    latex.DrawLatexNDC(0.15, ystart, Form("Ellipse 1: x=%.3f, y=%.3f, a=%.3f, b=%.3f, #theta=%.1f^{#circ}", v1[0], v1[1], v1[2], v1[3], v1[4]));
    latex.DrawLatexNDC(0.15, ystart - dy, Form("Ellipse 2: x=%.3f, y=%.3f, a=%.3f, b=%.3f, #theta=%.1f^{#circ}", v2[0], v2[1], v2[2], v2[3], v2[4]));
    latex.DrawLatexNDC(0.15, ystart - 2 * dy, Form("Combined: x=%.3f, y=%.3f, a=%.3f, b=%.3f, #theta=%.1f^{#circ}", e3->GetX1(), e3->GetY1(), e3->GetR1(), e3->GetR2(), e3->GetTheta()));

    c1->Update();
    c1->SaveAs(Savename);
}

void get_combined_contour()
{
    TFile *f1 = new TFile("./contourrootfile/everything.root", "READ");
    TFile *f2 = new TFile("./bestfittemplaterootfile/template.root", "READ");

    TVector2 *PbPb_0 = (TVector2 *)f2->Get("PbPb_local_min_cent_0");
    TVector2 *PbPb_1 = (TVector2 *)f2->Get("PbPb_local_min_cent_1");
    TVector2 *PbPb_2 = (TVector2 *)f2->Get("PbPb_local_min_cent_2");
    TVector2 *PbPb_3 = (TVector2 *)f2->Get("PbPb_local_min_cent_3");
    TVector2 *PbPb_4 = (TVector2 *)f2->Get("PbPb_local_min_cent_10");

    TVector2 *pp_0 = (TVector2 *)f2->Get("pp_local_min_period_22");

    TVectorD pp_sig1 = doeverything("pp, |#eta| < 2.4, Nominal, Period: (22)", "./contourtest/pp.png", f1, pp_0, 1);
    TVectorD pp_sig2 = doeverything("pp, |#eta| < 2.4, Nominal, Period: (22)", "./contourtest/pp.png", f1, pp_0, 2);

    TVectorD PbPb_sig_1_0 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0-10)", "./contourtest/PbPb_1.png", f1, PbPb_0, 1);
    TVectorD PbPb_sig_2_0 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0-10)", "./contourtest/PbPb_1.png", f1, PbPb_0, 2);

    TVectorD PbPb_sig_1_1 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (10-20)", "./contourtest/PbPb_2.png", f1, PbPb_1, 1);
    TVectorD PbPb_sig_2_1 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (10-20)", "./contourtest/PbPb_2.png", f1, PbPb_1, 2);

    TVectorD PbPb_sig_1_2 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (20-30)", "./contourtest/PbPb_3.png", f1, PbPb_2, 1);
    TVectorD PbPb_sig_2_2 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (20-30)", "./contourtest/PbPb_3.png", f1, PbPb_2, 2);

    TVectorD PbPb_sig_1_3 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (30-100)", "./contourtest/PbPb_4.png", f1, PbPb_3, 1);
    TVectorD PbPb_sig_2_3 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (30-100)", "./contourtest/PbPb_4.png", f1, PbPb_3, 2);

    TVectorD PbPb_sig_1_4 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0-100)", "./contourtest/PbPb_5.png", f1, PbPb_4, 1);
    TVectorD PbPb_sig_2_4 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0-100)", "./contourtest/PbPb_5.png", f1, PbPb_4, 2);

    DrawEllipsesComparison(PbPb_sig_1_0, pp_sig1, "./contourtest/PbPb_0_sig_1.png");
    DrawEllipsesComparison(PbPb_sig_2_0, pp_sig2, "./contourtest/PbPb_0_sig_2.png");

    DrawEllipsesComparison(PbPb_sig_1_1, pp_sig1, "./contourtest/PbPb_1_sig_1.png");
    DrawEllipsesComparison(PbPb_sig_2_1, pp_sig2, "./contourtest/PbPb_1_sig_2.png");

    DrawEllipsesComparison(PbPb_sig_1_2, pp_sig1, "./contourtest/PbPb_2_sig_1.png");
    DrawEllipsesComparison(PbPb_sig_2_2, pp_sig2, "./contourtest/PbPb_2_sig_2.png");

    DrawEllipsesComparison(PbPb_sig_1_3, pp_sig1, "./contourtest/PbPb_3_sig_1.png");
    DrawEllipsesComparison(PbPb_sig_2_3, pp_sig2, "./contourtest/PbPb_3_sig_2.png");

    DrawEllipsesComparison(PbPb_sig_1_4, pp_sig1, "./contourtest/PbPb_4_sig_1.png");
    DrawEllipsesComparison(PbPb_sig_2_4, pp_sig2, "./contourtest/PbPb_4_sig_2.png");

    // Testing area
    // First I am going to test PbPb - pp case
    TVectorD *test_ellipse = new TVectorD(5);
    (*test_ellipse)[0] = 0.0;
    (*test_ellipse)[1] = 0.0;
    (*test_ellipse)[2] = 10.0;
    (*test_ellipse)[3] = 5.0;
    (*test_ellipse)[4] = 0.0;
    TVectorD *test_circle = new TVectorD(5);
    (*test_circle)[0] = 0.0;
    (*test_circle)[1] = 0.0;
    (*test_circle)[2] = 5.0;
    (*test_circle)[3] = 5.0;
    (*test_circle)[4] = 0.0;

    DrawEllipsesComparison(*test_ellipse, *test_circle, "./contourtest/ellipse_circle_test.png");

    TVectorD *test_ellipse_45deg = new TVectorD(5);
    (*test_ellipse_45deg)[0] = 0.0;
    (*test_ellipse_45deg)[1] = 0.0;
    (*test_ellipse_45deg)[2] = 10.0;
    (*test_ellipse_45deg)[3] = 5.0;
    (*test_ellipse_45deg)[4] = 45.0;
    TVectorD *test_ellipse_135deg = new TVectorD(5);
    (*test_ellipse_135deg)[0] = 0.0;
    (*test_ellipse_135deg)[1] = 0.0;
    (*test_ellipse_135deg)[2] = 10.0;
    (*test_ellipse_135deg)[3] = 5.0;
    (*test_ellipse_135deg)[4] = 135.0;

    DrawEllipsesComparison(*test_ellipse_45deg, *test_ellipse_135deg, "./contourtest/ellipse45_ellipse135_test.png");

    TVectorD *test_circle_1 = new TVectorD(5);
    (*test_circle_1)[0] = 0.0;
    (*test_circle_1)[1] = 0.0;
    (*test_circle_1)[2] = 5.0;
    (*test_circle_1)[3] = 5.0;
    (*test_circle_1)[4] = 0.0;
    TVectorD *test_circle_2 = new TVectorD(5);
    (*test_circle_2)[0] = 0.0;
    (*test_circle_2)[1] = 0.0;
    (*test_circle_2)[2] = 7.0;
    (*test_circle_2)[3] = 7.0;
    (*test_circle_2)[4] = 0.0;

    DrawEllipsesComparison(*test_circle_1, *test_circle_2, "./contourtest/circle1_circle_2_test.png");
}