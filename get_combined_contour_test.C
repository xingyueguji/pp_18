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
void get_combined_contour_test()
{
    TH1D *PbPb_Nominal[2];
    TCanvas *c_PbPb_Nominal[2];

    TFile *f1 = new TFile("./contourrootfile/everything.root", "READ");

    c_PbPb_Nominal[0] = (TCanvas *)f1->Get("PbPb, |#eta| < 2.4, Nominal, centrality: (0-10)");
    // Access the list of primitives
    TList *primList_PbPb[2];
    TGraph *g_1sig_pbpb[2];
    TGraph *g_2sig_pbpb[2];

    primList_PbPb[0] = c_PbPb_Nominal[0]->GetListOfPrimitives();
    g_1sig_pbpb[0] = (TGraph *)primList_PbPb[0]->At(1);
    g_2sig_pbpb[0] = (TGraph *)primList_PbPb[0]->At(2);

    TVectorD conic_pbpb_1sig_1 = fit_ellipse(g_1sig_pbpb[0]);
    TVectorD conic_pbpb_2sig_1 = fit_ellipse(g_2sig_pbpb[0]);

    TVectorD ellipse_1sig_pbpb_1 = ConicToParametric(conic_pbpb_1sig_1);
    TVectorD ellipse_2sig_pbpb_1 = ConicToParametric(conic_pbpb_2sig_1);

    c_PbPb_Nominal[0]->cd();
    TEllipse *e = new TEllipse(ellipse_1sig_pbpb_1[0], ellipse_1sig_pbpb_1[1], // "x0", "y0"
                               ellipse_1sig_pbpb_1[2], ellipse_1sig_pbpb_1[3], // "a", "b"
                               0, 360,
                               ellipse_1sig_pbpb_1[4]); // "theta" (in degrees)
    e->SetFillStyle(0);                                 // hollow
    e->SetLineColor(kGreen);
    e->SetLineWidth(5);
    e->Draw("SAME");

    TEllipse *e1 = new TEllipse(ellipse_2sig_pbpb_1[0], ellipse_2sig_pbpb_1[1], // "x0", "y0"
                                ellipse_2sig_pbpb_1[2], ellipse_2sig_pbpb_1[3], // "a", "b"
                                0, 360,
                                ellipse_2sig_pbpb_1[4]); // "theta" (in degrees)
    e1->SetFillStyle(0);                                 // hollow
    e1->SetLineColor(kOrange);
    e1->SetLineWidth(5);
    e1->Draw("SAME");
    c_PbPb_Nominal[0]->GetListOfPrimitives()->ls();
    c_PbPb_Nominal[0]->SaveAs("./contourtest/test.png");

    c_PbPb_Nominal[1] = (TCanvas *)f1->Get("PbPb, |#eta| < 2.4, Nominal, centrality: (10-20)");
    primList_PbPb[1] = c_PbPb_Nominal[1]->GetListOfPrimitives();
    g_1sig_pbpb[1] = (TGraph *)primList_PbPb[1]->At(1);
    g_2sig_pbpb[1] = (TGraph *)primList_PbPb[1]->At(2);
    TVectorD conic_pbpb_1sig_2 = fit_ellipse(g_1sig_pbpb[1]);
    TVectorD conic_pbpb_2sig_2 = fit_ellipse(g_2sig_pbpb[1]);
    TVectorD ellipse_1sig_pbpb_2 = ConicToParametric(conic_pbpb_1sig_2);
    TVectorD ellipse_2sig_pbpb_2 = ConicToParametric(conic_pbpb_2sig_2);
    
    c_PbPb_Nominal[1]->cd();

    TEllipse *e3 = new TEllipse(ellipse_1sig_pbpb_2[0], ellipse_1sig_pbpb_2[1], // "x0", "y0"
                                ellipse_1sig_pbpb_2[2], ellipse_1sig_pbpb_2[3], // "a", "b"
                                0, 360,
                                ellipse_1sig_pbpb_2[4]); // "theta" (in degrees)
    e3->SetFillStyle(0);                                 // hollow
    e3->SetLineColor(kGreen);
    e3->SetLineWidth(5);
    e3->Draw("SAME");

    TEllipse *e4 = new TEllipse(ellipse_2sig_pbpb_2[0], ellipse_2sig_pbpb_2[1], // "x0", "y0"
                                ellipse_2sig_pbpb_2[2], ellipse_2sig_pbpb_2[3], // "a", "b"
                                0, 360,
                                ellipse_2sig_pbpb_2[4]); // "theta" (in degrees)
    e4->SetFillStyle(0);                                 // hollow
    e4->SetLineColor(kOrange);
    e4->SetLineWidth(5);
    e4->Draw("SAME");

    c_PbPb_Nominal[1]->SaveAs("./contourtest/test_1.png");
}