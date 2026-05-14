#include "tdrStyle.C"

TEllipse *CreateScaledEllipseFromSystematic(const TVector2 &syst_point,
                                            TEllipse *nominal_ellipse,
                                            const TString &type)
{
    // Extract nominal center and axes
    double x0 = nominal_ellipse->GetX1();
    double y0 = nominal_ellipse->GetY1();
    double a0 = nominal_ellipse->GetR1();
    double b0 = nominal_ellipse->GetR2();
    double theta = nominal_ellipse->GetTheta();

    if (type == "normal")
    {
        // Shift systematic point relative to nominal center
        double dx = syst_point.X() - x0;
        double dy = syst_point.Y() - y0;

        // Convert angle to radians
        double theta_rad = theta * TMath::DegToRad();

        // Rotate the point into the ellipse's local frame
        double x_rot = dx * std::cos(theta_rad) + dy * std::sin(theta_rad);
        double y_rot = -dx * std::sin(theta_rad) + dy * std::cos(theta_rad);

        // Compute distance in ellipse units
        double u = x_rot / a0;
        double v = y_rot / b0;
        double scaling = std::sqrt(u * u + v * v);

        auto *new_ellipse = new TEllipse(x0, y0, a0 * scaling, b0 * scaling, 0, 360, theta);
        new_ellipse->SetLineStyle(2);
        new_ellipse->SetLineColor(kBlue + 2);
        new_ellipse->SetLineWidth(2);
        new_ellipse->SetFillStyle(0);

        return new_ellipse;
    }
    else if (type == "degen")
    {
        // Direction vector from nominal to systematic
        double dx = syst_point.X() - x0;
        double dy = syst_point.Y() - y0;

        double half_length = std::sqrt(dx * dx + dy * dy);
        double angle_deg = std::atan2(dy, dx) * TMath::RadToDeg();

        auto *new_ellipse = new TEllipse(x0, y0, half_length, 0.001, 0, 360, angle_deg);
        new_ellipse->SetLineStyle(2);
        new_ellipse->SetLineColor(kRed + 1);
        new_ellipse->SetLineWidth(2);
        new_ellipse->SetFillStyle(0);

        return new_ellipse;
    }
    else
    {
        Error("CreateScaledEllipseFromSystematic", "Unknown type '%s'", type.Data());
        return nullptr;
    }
}

double GetSigmaXFromEllipse(const TVectorD &v, double scale = std::sqrt(2.30))
{
    double a = v[2];
    double b = v[3];
    double th = v[4] * TMath::DegToRad();

    double c = std::cos(th);
    double s = std::sin(th);

    double sigxx = c * c * a * a + s * s * b * b;
    return std::sqrt(std::max(0.0, sigxx)) / scale;
}

double GetSigmaYFromEllipse(const TVectorD &v, double scale = std::sqrt(2.30))
{
    double a = v[2];
    double b = v[3];
    double th = v[4] * TMath::DegToRad();

    double c = std::cos(th);
    double s = std::sin(th);

    double sigyy = s * s * a * a + c * c * b * b;
    return std::sqrt(std::max(0.0, sigyy)) / scale;
}

TVectorD EllipseToVector(TEllipse *e)
{
    TVectorD v(5);
    v[0] = e->GetX1();
    v[1] = e->GetY1();
    v[2] = e->GetR1();
    v[3] = e->GetR2();
    v[4] = e->GetTheta();
    return v;
}

TGraphErrors *MakeCentralityTGraph(const char *name,
                                   const std::vector<TVectorD> &statEllipses,
                                   const std::vector<TVectorD> &systEllipses,
                                   bool useX)
{
    static const char *centLabels[5] = {
        "0-100%", "0-10%", "10-20%", "20-30%", "30-100%"};

    int n = statEllipses.size();
    auto *g = new TGraphErrors(n);
    g->SetName(name);

    std::cout << "==============================" << std::endl;
    std::cout << "This is " << name << std::endl;

    for (int i = 0; i < n; ++i)
    {
        double stat = useX ? GetSigmaXFromEllipse(statEllipses[i])
                           : GetSigmaYFromEllipse(statEllipses[i]);
        double syst = useX ? GetSigmaXFromEllipse(systEllipses[i])
                           : GetSigmaYFromEllipse(systEllipses[i]);

        double val = syst; // <-- this is the key change

        g->SetPoint(i, i, val);
        g->SetPointError(i, 0.0, stat);

        std::cout << centLabels[i] << " : value = " << val << std::endl;
        std::cout << centLabels[i] << " : stat  = " << stat << std::endl;
        std::cout << centLabels[i] << " : syst  = " << syst << std::endl;
    }

    std::cout << std::endl;
    return g;
}

std::vector<TVectorD> BuildSystematicEllipseSet(const std::vector<TVector2 *> &systPts,
                                                const std::vector<TVectorD> &nominalEllipses,
                                                const TString &type)
{
    std::vector<TVectorD> out;
    out.reserve(systPts.size());

    for (size_t i = 0; i < systPts.size(); ++i)
    {
        TEllipse nominal(nominalEllipses[i][0],
                         nominalEllipses[i][1],
                         nominalEllipses[i][2],
                         nominalEllipses[i][3],
                         0, 360,
                         nominalEllipses[i][4]);

        TEllipse *esys = CreateScaledEllipseFromSystematic(*systPts[i], &nominal, type);
        out.push_back(EllipseToVector(esys));
        delete esys;
    }

    return out;
}

std::vector<TVectorD> BuildZeroSystematicEllipseSet(const std::vector<TVectorD> &nominalEllipses)
{
    std::vector<TVectorD> out;
    out.reserve(nominalEllipses.size());

    for (const auto &v : nominalEllipses)
    {
        TVectorD z(5);
        z[0] = v[0];
        z[1] = v[1];
        z[2] = 0.0;
        z[3] = 0.0;
        z[4] = v[4];
        out.push_back(z);
    }

    return out;
}

void SaveAllCentralityTGraphs(const TString &outfile,
                              const TString &type,

                              const std::vector<TVectorD> &pbpbEll,
                              const std::vector<TVectorD> &ppEll,
                              const std::vector<TVectorD> &subEll,

                              const std::vector<TVector2 *> &pbpb_nominal,
                              const std::vector<TVector2 *> &pbpb_tnpU,
                              const std::vector<TVector2 *> &pbpb_tnpD,
                              const std::vector<TVector2 *> &pbpb_acoup,
                              const std::vector<TVector2 *> &pbpb_acodown,
                              const std::vector<TVector2 *> &pbpb_nobk,
                              const std::vector<TVector2 *> &pbpb_massrange,
                              const std::vector<TVector2 *> &pbpb_HFup,
                              const std::vector<TVector2 *> &pbpb_HFdown,
                              const std::vector<TVector2 *> &pbpb_1DpT,

                              const std::vector<TVector2 *> &pp_nominal,
                              const std::vector<TVector2 *> &pp_tnpU,
                              const std::vector<TVector2 *> &pp_tnpD,
                              const std::vector<TVector2 *> &pp_acoup,
                              const std::vector<TVector2 *> &pp_acodown,
                              const std::vector<TVector2 *> &pp_nobk,
                              const std::vector<TVector2 *> &pp_massrange,
                              const std::vector<TVector2 *> &pp_HFup,
                              const std::vector<TVector2 *> &pp_HFdown,
                              const std::vector<TVector2 *> &pp_1DpT,

                              const std::vector<TVector2 *> &sub_nominal,
                              const std::vector<TVector2 *> &sub_tnpU,
                              const std::vector<TVector2 *> &sub_tnpD,
                              const std::vector<TVector2 *> &sub_acoup,
                              const std::vector<TVector2 *> &sub_acodown,
                              const std::vector<TVector2 *> &sub_nobk,
                              const std::vector<TVector2 *> &sub_massrange,
                              const std::vector<TVector2 *> &sub_HFup,
                              const std::vector<TVector2 *> &sub_HFdown,
                              const std::vector<TVector2 *> &sub_1DpT)
{
    TFile *fout = new TFile(outfile, "RECREATE");

    // ---- build systematic ellipse sets ----
    auto pbpb_nominal_syst = BuildZeroSystematicEllipseSet(pbpbEll);
    auto pbpb_tnpU_syst = BuildSystematicEllipseSet(pbpb_tnpU, pbpbEll, type);
    auto pbpb_tnpD_syst = BuildSystematicEllipseSet(pbpb_tnpD, pbpbEll, type);
    auto pbpb_acoup_syst = BuildSystematicEllipseSet(pbpb_acoup, pbpbEll, type);
    auto pbpb_acodown_syst = BuildSystematicEllipseSet(pbpb_acodown, pbpbEll, type);
    auto pbpb_nobk_syst = BuildSystematicEllipseSet(pbpb_nobk, pbpbEll, type);
    auto pbpb_massrange_syst = BuildSystematicEllipseSet(pbpb_massrange, pbpbEll, type);
    auto pbpb_HFup_syst = BuildSystematicEllipseSet(pbpb_HFup, pbpbEll, type);
    auto pbpb_HFdown_syst = BuildSystematicEllipseSet(pbpb_HFdown, pbpbEll, type);
    auto pbpb_1DpT_syst = BuildSystematicEllipseSet(pbpb_1DpT, pbpbEll, type);

    auto pp_nominal_syst = BuildZeroSystematicEllipseSet(ppEll);
    auto pp_tnpU_syst = BuildSystematicEllipseSet(pp_tnpU, ppEll, type);
    auto pp_tnpD_syst = BuildSystematicEllipseSet(pp_tnpD, ppEll, type);
    auto pp_acoup_syst = BuildSystematicEllipseSet(pp_acoup, ppEll, type);
    auto pp_acodown_syst = BuildSystematicEllipseSet(pp_acodown, ppEll, type);
    auto pp_nobk_syst = BuildSystematicEllipseSet(pp_nobk, ppEll, type);
    auto pp_massrange_syst = BuildSystematicEllipseSet(pp_massrange, ppEll, type);
    auto pp_HFup_syst = BuildSystematicEllipseSet(pp_HFup, ppEll, type);
    auto pp_HFdown_syst = BuildSystematicEllipseSet(pp_HFdown, ppEll, type);
    auto pp_1DpT_syst = BuildSystematicEllipseSet(pp_1DpT, ppEll, type);

    auto sub_nominal_syst = BuildZeroSystematicEllipseSet(subEll);
    auto sub_tnpU_syst = BuildSystematicEllipseSet(sub_tnpU, subEll, type);
    auto sub_tnpD_syst = BuildSystematicEllipseSet(sub_tnpD, subEll, type);
    auto sub_acoup_syst = BuildSystematicEllipseSet(sub_acoup, subEll, type);
    auto sub_acodown_syst = BuildSystematicEllipseSet(sub_acodown, subEll, type);
    auto sub_nobk_syst = BuildSystematicEllipseSet(sub_nobk, subEll, type);
    auto sub_massrange_syst = BuildSystematicEllipseSet(sub_massrange, subEll, type);
    auto sub_HFup_syst = BuildSystematicEllipseSet(sub_HFup, subEll, type);
    auto sub_HFdown_syst = BuildSystematicEllipseSet(sub_HFdown, subEll, type);
    auto sub_1DpT_syst = BuildSystematicEllipseSet(sub_1DpT, subEll, type);

    // ---------- dM ----------
    MakeCentralityTGraph("HI_dM_chi2_raw_nominal", pbpbEll, pbpb_nominal_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_tnpU", pbpbEll, pbpb_tnpU_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_tnpD", pbpbEll, pbpb_tnpD_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_acoup", pbpbEll, pbpb_acoup_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_acodown", pbpbEll, pbpb_acodown_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_nominal_no_bk", pbpbEll, pbpb_nobk_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_nominal_mass_range", pbpbEll, pbpb_massrange_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_HF_up", pbpbEll, pbpb_HFup_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_HF_down", pbpbEll, pbpb_HFdown_syst, true)->Write();
    MakeCentralityTGraph("HI_dM_chi2_raw_1D_pT", pbpbEll, pbpb_1DpT_syst, true)->Write();

    MakeCentralityTGraph("HI_pp_dM_chi2_raw_nominal", ppEll, pp_nominal_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_tnpU", ppEll, pp_tnpU_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_tnpD", ppEll, pp_tnpD_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_acoup", ppEll, pp_acoup_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_acodown", ppEll, pp_acodown_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_nominal_no_bk", ppEll, pp_nobk_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_nominal_mass_range", ppEll, pp_massrange_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_HF_up", ppEll, pp_HFup_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_HF_down", ppEll, pp_HFdown_syst, true)->Write();
    MakeCentralityTGraph("HI_pp_dM_chi2_raw_1D_pT", ppEll, pp_1DpT_syst, true)->Write();

    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_nominal", subEll, sub_nominal_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_tnpU", subEll, sub_tnpU_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_tnpD", subEll, sub_tnpD_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_acoup", subEll, sub_acoup_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_acodown", subEll, sub_acodown_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_nominal_no_bk", subEll, sub_nobk_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_nominal_mass_range", subEll, sub_massrange_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_HF_up", subEll, sub_HFup_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_HF_down", subEll, sub_HFdown_syst, true)->Write();
    MakeCentralityTGraph("HI_sub_pp_dM_chi2_raw_1D_pT", subEll, sub_1DpT_syst, true)->Write();

    // ---------- dWidth ----------
    MakeCentralityTGraph("HI_dWidth_chi2_raw_nominal", pbpbEll, pbpb_nominal_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_tnpU", pbpbEll, pbpb_tnpU_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_tnpD", pbpbEll, pbpb_tnpD_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_acoup", pbpbEll, pbpb_acoup_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_acodown", pbpbEll, pbpb_acodown_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_nominal_no_bk", pbpbEll, pbpb_nobk_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_nominal_mass_range", pbpbEll, pbpb_massrange_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_HF_up", pbpbEll, pbpb_HFup_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_HF_down", pbpbEll, pbpb_HFdown_syst, false)->Write();
    MakeCentralityTGraph("HI_dWidth_chi2_raw_1D_pT", pbpbEll, pbpb_1DpT_syst, false)->Write();

    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_nominal", ppEll, pp_nominal_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_tnpU", ppEll, pp_tnpU_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_tnpD", ppEll, pp_tnpD_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_acoup", ppEll, pp_acoup_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_acodown", ppEll, pp_acodown_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_nominal_no_bk", ppEll, pp_nobk_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_nominal_mass_range", ppEll, pp_massrange_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_HF_up", ppEll, pp_HFup_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_HF_down", ppEll, pp_HFdown_syst, false)->Write();
    MakeCentralityTGraph("HI_pp_dWidth_chi2_raw_1D_pT", ppEll, pp_1DpT_syst, false)->Write();

    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_nominal", subEll, sub_nominal_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_tnpU", subEll, sub_tnpU_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_tnpD", subEll, sub_tnpD_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_acoup", subEll, sub_acoup_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_acodown", subEll, sub_acodown_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_nominal_no_bk", subEll, sub_nobk_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_nominal_mass_range", subEll, sub_massrange_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_HF_up", subEll, sub_HFup_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_HF_down", subEll, sub_HFdown_syst, false)->Write();
    MakeCentralityTGraph("HI_sub_pp_dWidth_chi2_raw_1D_pT", subEll, sub_1DpT_syst, false)->Write();

    fout->Close();
}

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

    // cout << "x is " << ellipse[0] << " y is " << ellipse[1] << " semimajor x is " << ellipse[2] << " semimajor y is " << ellipse[3] << " angle is " << ellipse[4] << endl;

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

    /*for (int i = 0; i < g->GetN(); ++i)
    {
        double x, y;
        g->GetPoint(i, x, y);
        TString label = Form("%.3f, %.3f", x, y);
        latex1.DrawLatex(x + 0.01, y + 0.1, label);
    }*/

    TGraph *g_1 = new TGraph();
    g_1->SetPoint(0, ellipse_1[0], ellipse_1[1]);
    g_1->SetMarkerStyle(20); // Solid circle
    g_1->SetMarkerSize(1.5);
    g_1->SetMarkerColor(kBlue);
    g_1->Draw("P SAME");

    /*TGraph *g_2 = new TGraph();
    g_2->SetPoint(0, ellipse_2[0], ellipse_2[1]);
    g_2->SetMarkerStyle(20); // Solid circle
    g_2->SetMarkerSize(1.5);
    g_2->SetMarkerColor(kOrange);
    g_2->Draw("P SAME");*/

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

    /*std::cout << "covSum matrix elements:" << std::endl;
    for (int i = 0; i < covSum.GetNrows(); ++i)
    {
        for (int j = 0; j < covSum.GetNcols(); ++j)
        {
            std::cout << covSum(i, j) << "\t";
        }
        std::cout << std::endl;
    }*/

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

    double xc, yc = 0;

    if (v1[0] == v2[0] && v1[1] == v2[1])
    {
        xc = v1[0];
        yc = v1[1];
    }
    else
    {
        xc = v1[0] - v2[0];
        yc = v1[1] - v2[1];
    }

    TEllipse *result = new TEllipse(xc, yc, a_comb, b_comb, 0, 360, theta_comb);
    result->SetLineColor(kRed);
    result->SetLineWidth(2);

    return result;
}

void DrawEllipsesComparison(const TVectorD &v1, const TVectorD &v2, TString Savename, const TVector2 &sys_1, const TVector2 &sys_2, const TVector2 &sys_3,
                            const TVector2 &sys_4, const TVector2 &sys_5, const TVector2 &sys_6, const TVector2 &sys_7, const TVector2 &sys_8, const TVector2 &pp_sys_1, const TVector2 &pp_sys_2, const TVector2 &pp_sys_3, const TVector2 &pp_sys_4, const TVector2 &pp_sys_5,
                            const TVector2 &pp_sys_6, const TVector2 &pp_sys_7, const TString &type, const TString &name)
{
    auto EllipseToVector = [](TEllipse *e) -> TVectorD
    {
        TVectorD v(5);
        v[0] = e->GetX1();
        v[1] = e->GetY1();
        v[2] = e->GetR1();
        v[3] = e->GetR2();
        v[4] = e->GetTheta();
        return v;
    };
    auto CreateAxisLine = [](double angle_deg, double xmin, double xmax, double ymin, double ymax)
    {
        double angle_rad = angle_deg * TMath::DegToRad();
        double dx = std::cos(angle_rad);
        double dy = std::sin(angle_rad);
        double length = std::max(xmax - xmin, ymax - ymin) * 2.0;

        auto line = new TLine(-length * dx, -length * dy, length * dx, length * dy);
        line->SetLineStyle(2); // dashed
        line->SetLineColor(kGray + 2);
        line->SetLineWidth(2);
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

    auto CreateScaledEllipseFromSystematic = [](const TVector2 &syst_point, TEllipse *nominal_ellipse, const TString &type) -> TEllipse *
    {
        // Extract nominal center and axes
        double x0 = nominal_ellipse->GetX1();
        double y0 = nominal_ellipse->GetY1();
        double a0 = nominal_ellipse->GetR1();       // semi-major axis
        double b0 = nominal_ellipse->GetR2();       // semi-minor axis
        double theta = nominal_ellipse->GetTheta(); // degrees

        if (type == "normal")
        {
            // Shift systematic point relative to nominal center
            double dx = syst_point.X() - x0;
            double dy = syst_point.Y() - y0;

            // Convert angle to radians
            double theta_rad = theta * TMath::DegToRad();

            // Rotate the point into the ellipse's local frame
            double x_rot = dx * std::cos(theta_rad) + dy * std::sin(theta_rad);
            double y_rot = -dx * std::sin(theta_rad) + dy * std::cos(theta_rad);

            // Compute distance in "ellipse units"
            double u = x_rot / a0;
            double v = y_rot / b0;
            double scaling = std::sqrt(u * u + v * v);

            // Create new ellipse with scaled axes
            auto new_ellipse = new TEllipse(x0, y0, a0 * scaling, b0 * scaling, 0, 360, theta);
            new_ellipse->SetLineStyle(2); // dashed
            new_ellipse->SetLineColor(kBlue + 2);
            new_ellipse->SetLineWidth(2);
            new_ellipse->SetFillStyle(0); // hollow

            return new_ellipse;
        }
        if (type == "degen")
        {
            // Direction vector from nominal to systematic
            double dx = syst_point.X() - x0;
            double dy = syst_point.Y() - y0;

            // Length from center to syst point
            double half_length = std::sqrt(dx * dx + dy * dy);

            // Angle of the line in degrees
            double angle_deg = std::atan2(dy, dx) * TMath::RadToDeg();

            // Create a degenerate ellipse (minor axis = 0) oriented along that direction
            auto new_ellipse = new TEllipse(x0, y0, half_length, 0.001, 0, 360, angle_deg);
            new_ellipse->SetLineStyle(2); // dashed
            new_ellipse->SetLineColor(kRed + 1);
            new_ellipse->SetLineWidth(2);
            new_ellipse->SetFillStyle(0); // hollow
            return new_ellipse;
        }
        else
        {
            Error("CreateScaledEllipseFromSystematic", "Unknown type '%s'", type.Data());
            return nullptr;
        }
    };

    auto CreateCombinedEnvelope = [](const TVectorD &nominal,
                                     const std::vector<std::pair<TEllipse *, TEllipse *>> &updown_pairs,
                                     const std::vector<TEllipse *> &single_systs) -> TEllipse *
    {
        auto EllipseToVector = [](TEllipse *e) -> TVectorD
        {
            TVectorD v(5);
            v[0] = e->GetX1();    // x center
            v[1] = e->GetY1();    // y center
            v[2] = e->GetR1();    // major axis
            v[3] = e->GetR2();    // minor axis
            v[4] = e->GetTheta(); // angle (deg)
            return v;
        };

        std::vector<TVectorD> syst_vectors;

        // Handle up/down pairs
        for (const auto &pair : updown_pairs)
        {
            if (!pair.first || !pair.second)
                continue;

            TVectorD v_up = EllipseToVector(pair.first);
            TVectorD v_down = EllipseToVector(pair.second);

            TVectorD v_max(5);
            v_max[0] = nominal[0];                   // center x
            v_max[1] = nominal[1];                   // center y
            v_max[2] = std::max(v_up[2], v_down[2]); // max major
            v_max[3] = std::max(v_up[3], v_down[3]); // max minor
            v_max[4] = nominal[4];                   // angle from up (or average if needed)

            syst_vectors.push_back(v_max);
        }

        // Handle single-sided systematics
        for (TEllipse *e : single_systs)
        {
            if (!e)
                continue;
            syst_vectors.push_back(EllipseToVector(e));
        }

        // Combine all systematics with the nominal
        TVectorD combined = nominal;

        for (const auto &v_syst : syst_vectors)
        {
            // cout << "BEfore : x is " << v_syst[0] << " y is " << v_syst[1] << " semimajor x is " << v_syst[2] << " semimajor y is " << v_syst[3] << " angle is " << v_syst[4] << endl;

            TEllipse *combined_ellipse = CombineEllipsesFromVectors(combined, v_syst);

            // cout << "after : x is " << combined_ellipse->GetX1() << " y is " << combined_ellipse->GetY1() << " semimajor x is " << combined_ellipse->GetR1() << " semimajor y is " << combined_ellipse->GetR2() << " angle is " << combined_ellipse->GetTheta() << endl;

            // Update combined vector
            TVectorD temp(5);
            temp[0] = combined_ellipse->GetX1();
            temp[1] = combined_ellipse->GetY1();
            temp[2] = combined_ellipse->GetR1();
            temp[3] = combined_ellipse->GetR2();
            temp[4] = combined_ellipse->GetTheta();
            combined = temp;
            delete combined_ellipse;
        }

        // Final ellipse
        TEllipse *result = new TEllipse(combined[0], combined[1], combined[2], combined[3], 0, 360, combined[4]);
        result->SetLineColor(kMagenta + 2);
        result->SetLineStyle(2);
        result->SetLineWidth(3);
        result->SetFillStyle(0);
        return result;
    };

    auto CombiningSystematicOnly = [&](const std::vector<std::pair<TEllipse *, TEllipse *>> &pairs,
                                       const std::vector<TEllipse *> &singles) -> TEllipse *
    {
        std::vector<TVectorD> selected_vectors;

        // 1. From up/down pairs, take max axes
        for (const auto &pair : pairs)
        {
            if (!pair.first || !pair.second)
                continue;

            TVectorD up = EllipseToVector(pair.first);
            TVectorD down = EllipseToVector(pair.second);

            TVectorD v_max(5);
            v_max[0] = up[0];                    // x-center
            v_max[1] = up[1];                    // y-center
            v_max[2] = std::max(up[2], down[2]); // max semi-major
            v_max[3] = std::max(up[3], down[3]); // max semi-minor
            v_max[4] = up[4];                    // use up’s angle (you could average if desired)

            selected_vectors.push_back(v_max);
        }

        // 2. From single systematics
        for (TEllipse *e : singles)
        {
            if (!e)
                continue;
            selected_vectors.push_back(EllipseToVector(e));
        }

        // 3. Combine all selected systematics
        if (selected_vectors.empty())
            return nullptr;

        TVectorD combined = selected_vectors[0];
        for (size_t i = 1; i < selected_vectors.size(); ++i)
        {
            TEllipse *temp = CombineEllipsesFromVectors(combined, selected_vectors[i]);

            // Convert combined TEllipse back to vector
            TVectorD temp_v(5);
            temp_v[0] = temp->GetX1();
            temp_v[1] = temp->GetY1();
            temp_v[2] = temp->GetR1();
            temp_v[3] = temp->GetR2();
            temp_v[4] = temp->GetTheta();

            combined = temp_v;
            delete temp;
        }

        // 4. Return final TEllipse
        TEllipse *final = new TEllipse(combined[0], combined[1], combined[2], combined[3], 0, 360, combined[4]);
        final->SetLineStyle(2);
        final->SetLineColor(kRed + 1);
        final->SetLineWidth(3);
        final->SetFillStyle(0);

        return final;
    };

    auto PrintXYErrorsFromEllipse = [](TEllipse *e, TString x)
    {
        if (!e)
        {
            std::cout << "Null ellipse pointer!" << std::endl;
            return;
        }

        double xc = e->GetX1();
        double yc = e->GetY1();
        const double a = e->GetR1();
        const double b = e->GetR2();
        const double th = e->GetTheta() * TMath::DegToRad(); // radians

        const double c = std::cos(th);
        const double s = std::sin(th);

        // Covariance matrix in (x,y) basis:
        // Σ = R diag(a², b²) Rᵀ
        const double a2 = a * a;
        const double b2 = b * b;
        const double sigxx = c * c * a2 + s * s * b2;
        const double sigyy = s * s * a2 + c * c * b2;

        double sig_x = std::sqrt(std::max(0.0, sigxx));
        double sig_y = std::sqrt(std::max(0.0, sigyy));

        cout << "This is " << x << endl;
        std::cout << "Ellipse center (x,y) = (" << xc << ", " << yc << ")\n";
        std::cout << "Uncertainty in x (sigx) = " << sig_x / 1.515 << "\n";
        std::cout << "Uncertainty in y (sigy) = " << sig_y / 1.515 << "\n";
        std::cout << "\n"
                  << endl;
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

    x_half *= 1.5;
    y_half *= 1.5;

    // Use the larger half-range for both x and y
    double range_half = std::max(x_half, y_half);

    x_min = x_center - range_half;
    x_max = x_center + range_half;
    y_min = y_center - range_half;
    y_max = y_center + range_half;

    TCanvas *c1 = new TCanvas("c1", "Ellipse Comparison", 800, 800);
    c1->DrawFrame(x_min, y_min, x_max, y_max, "Contours;#DeltaM(GeV);#DeltaWidth(GeV)");

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
    perp_a->SetLineWidth(3);
    perp_b->SetLineColor(kRed + 2);
    perp_b->SetLineStyle(3);
    perp_b->SetLineWidth(3);
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

    // PbPb systematic

    // tnpU

    TMarker *PbPb_1 = new TMarker(sys_1.X(), sys_1.Y(), 20);
    PbPb_1->SetMarkerColor(kRed);
    PbPb_1->SetMarkerSize(1.2);
    PbPb_1->Draw("same");

    TEllipse *PbPb_1_sys_contour = CreateScaledEllipseFromSystematic(sys_1, e1, type);
    PbPb_1_sys_contour->SetLineColor(kRed);
    PbPb_1_sys_contour->SetLineWidth(2);
    PbPb_1_sys_contour->SetLineStyle(1);
    PbPb_1_sys_contour->SetFillStyle(0);
    PbPb_1_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_1_sys_contour, "PbPb TnP U");

    // tnpD

    TMarker *PbPb_2 = new TMarker(sys_2.X(), sys_2.Y(), 20);
    PbPb_2->SetMarkerColor(kBlue);
    PbPb_2->SetMarkerSize(1.2);
    PbPb_2->Draw("same");

    TEllipse *PbPb_2_sys_contour = CreateScaledEllipseFromSystematic(sys_2, e1, type);
    PbPb_2_sys_contour->SetLineColor(kBlue);
    PbPb_2_sys_contour->SetLineWidth(2);
    PbPb_2_sys_contour->SetLineStyle(1);
    PbPb_2_sys_contour->SetFillStyle(0);
    PbPb_2_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_2_sys_contour, "PbPb TnP D");

    // AcoUp

    TMarker *PbPb_3 = new TMarker(sys_3.X(), sys_3.Y(), 20);
    PbPb_3->SetMarkerColor(kGreen + 2);
    PbPb_3->SetMarkerSize(1.2);
    PbPb_3->Draw("same");

    TEllipse *PbPb_3_sys_contour = CreateScaledEllipseFromSystematic(sys_3, e1, type);
    PbPb_3_sys_contour->SetLineColor(kGreen + 2);
    PbPb_3_sys_contour->SetLineWidth(2);
    PbPb_3_sys_contour->SetLineStyle(1);
    PbPb_3_sys_contour->SetFillStyle(0);
    PbPb_3_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_3_sys_contour, "PbPb Aco Up");

    // AcoDown

    TMarker *PbPb_4 = new TMarker(sys_4.X(), sys_4.Y(), 20);
    PbPb_4->SetMarkerColor(kMagenta);
    PbPb_4->SetMarkerSize(1.2);
    PbPb_4->Draw("same");

    TEllipse *PbPb_4_sys_contour = CreateScaledEllipseFromSystematic(sys_4, e1, type);
    PbPb_4_sys_contour->SetLineColor(kMagenta);
    PbPb_4_sys_contour->SetLineWidth(2);
    PbPb_4_sys_contour->SetLineStyle(1);
    PbPb_4_sys_contour->SetFillStyle(0);
    PbPb_4_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_4_sys_contour, "PbPb Aco Down");

    // nobk

    TMarker *PbPb_5 = new TMarker(sys_5.X(), sys_5.Y(), 20);
    PbPb_5->SetMarkerColor(kGray + 2);
    PbPb_5->SetMarkerSize(1.2);
    PbPb_5->Draw("same");

    TEllipse *PbPb_5_sys_contour = CreateScaledEllipseFromSystematic(sys_5, e1, type);
    PbPb_5_sys_contour->SetLineColor(kGray + 2);
    PbPb_5_sys_contour->SetLineWidth(2);
    PbPb_5_sys_contour->SetLineStyle(1);
    PbPb_5_sys_contour->SetFillStyle(0);
    PbPb_5_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_5_sys_contour, "PbPb nobk");

    // massrange

    TMarker *PbPb_6 = new TMarker(sys_6.X(), sys_6.Y(), 20);
    PbPb_6->SetMarkerColor(kOrange + 7);
    PbPb_6->SetMarkerSize(1.2);
    PbPb_6->Draw("same");

    TEllipse *PbPb_6_sys_contour = CreateScaledEllipseFromSystematic(sys_6, e1, type);
    PbPb_6_sys_contour->SetLineColor(kOrange + 7);
    PbPb_6_sys_contour->SetLineWidth(2);
    PbPb_6_sys_contour->SetLineStyle(1);
    PbPb_6_sys_contour->SetFillStyle(0);
    PbPb_6_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_6_sys_contour, "PbPb mass range");

    // HFup

    TMarker *PbPb_7 = new TMarker(sys_7.X(), sys_7.Y(), 20);
    PbPb_7->SetMarkerColor(kCyan);
    PbPb_7->SetMarkerSize(1.2);
    PbPb_7->Draw("same");

    TEllipse *PbPb_7_sys_contour = CreateScaledEllipseFromSystematic(sys_7, e1, type);
    PbPb_7_sys_contour->SetLineColor(kCyan);
    PbPb_7_sys_contour->SetLineWidth(2);
    PbPb_7_sys_contour->SetLineStyle(1);
    PbPb_7_sys_contour->SetFillStyle(0);
    PbPb_7_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_7_sys_contour, "PbPb HF up");

    // HFdown

    TMarker *PbPb_8 = new TMarker(sys_8.X(), sys_8.Y(), 20);
    PbPb_8->SetMarkerColor(kBlack);
    PbPb_8->SetMarkerSize(1.2);
    PbPb_8->Draw("same");

    TEllipse *PbPb_8_sys_contour = CreateScaledEllipseFromSystematic(sys_8, e1, type);
    PbPb_8_sys_contour->SetLineColor(kBlack);
    PbPb_8_sys_contour->SetLineWidth(2);
    PbPb_8_sys_contour->SetLineStyle(1);
    PbPb_8_sys_contour->SetFillStyle(0);
    PbPb_8_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_8_sys_contour, "PbPb HF down");

    // TnpU

    TMarker *pp_1 = new TMarker(pp_sys_1.X(), pp_sys_1.Y(), 20);
    pp_1->SetMarkerColor(kRed);
    pp_1->SetMarkerSize(1.2);
    pp_1->Draw("same");

    TEllipse *pp_1_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_1, e2, type);
    pp_1_sys_contour->SetLineColor(kRed);
    pp_1_sys_contour->SetLineWidth(2);
    pp_1_sys_contour->SetLineStyle(1);
    pp_1_sys_contour->SetFillStyle(0);
    pp_1_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_1_sys_contour, "pp TnP U");

    // TnpD

    TMarker *pp_2 = new TMarker(pp_sys_2.X(), pp_sys_2.Y(), 20);
    pp_2->SetMarkerColor(kBlue);
    pp_2->SetMarkerSize(1.2);
    pp_2->Draw("same");

    TEllipse *pp_2_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_2, e2, type);
    pp_2_sys_contour->SetLineColor(kBlue);
    pp_2_sys_contour->SetLineWidth(2);
    pp_2_sys_contour->SetLineStyle(1);
    pp_2_sys_contour->SetFillStyle(0);
    pp_2_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_2_sys_contour, "pp TnP D");

    // AcoUp

    TMarker *pp_3 = new TMarker(pp_sys_3.X(), pp_sys_3.Y(), 20);
    pp_3->SetMarkerColor(kGreen + 2);
    pp_3->SetMarkerSize(1.2);
    pp_3->Draw("same");

    TEllipse *pp_3_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_3, e2, type);
    pp_3_sys_contour->SetLineColor(kGreen + 2);
    pp_3_sys_contour->SetLineWidth(2);
    pp_3_sys_contour->SetLineStyle(1);
    pp_3_sys_contour->SetFillStyle(0);
    pp_3_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_3_sys_contour, "pp Aco Up");

    // AcoDown

    TMarker *pp_4 = new TMarker(pp_sys_4.X(), pp_sys_4.Y(), 20);
    pp_4->SetMarkerColor(kMagenta);
    pp_4->SetMarkerSize(1.2);
    pp_4->Draw("same");

    TEllipse *pp_4_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_4, e2, type);
    pp_4_sys_contour->SetLineColor(kMagenta);
    pp_4_sys_contour->SetLineWidth(2);
    pp_4_sys_contour->SetLineStyle(1);
    pp_4_sys_contour->SetFillStyle(0);
    pp_4_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_4_sys_contour, "pp Aco Down");

    // nobk

    TMarker *pp_5 = new TMarker(pp_sys_5.X(), pp_sys_5.Y(), 20);
    pp_5->SetMarkerColor(kGray + 2);
    pp_5->SetMarkerSize(1.2);
    pp_5->Draw("same");

    TEllipse *pp_5_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_5, e2, type);
    pp_5_sys_contour->SetLineColor(kGray + 2);
    pp_5_sys_contour->SetLineWidth(2);
    pp_5_sys_contour->SetLineStyle(1);
    pp_5_sys_contour->SetFillStyle(0);
    pp_5_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_5_sys_contour, "pp nobk");

    // massrange

    TMarker *pp_6 = new TMarker(pp_sys_6.X(), pp_sys_6.Y(), 20);
    pp_6->SetMarkerColor(kOrange + 7);
    pp_6->SetMarkerSize(1.2);
    pp_6->Draw("same");

    TEllipse *pp_6_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_6, e2, type);
    pp_6_sys_contour->SetLineColor(kOrange + 7);
    pp_6_sys_contour->SetLineWidth(2);
    pp_6_sys_contour->SetLineStyle(1);
    pp_6_sys_contour->SetFillStyle(0);
    pp_6_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_6_sys_contour, "pp mass range");

    // 1D_pT

    TMarker *pp_7 = new TMarker(pp_sys_7.X(), pp_sys_7.Y(), 20);
    pp_7->SetMarkerColor(kCyan + 3);
    pp_7->SetMarkerSize(1.2);
    pp_7->Draw("same");

    TEllipse *pp_7_sys_contour = CreateScaledEllipseFromSystematic(pp_sys_7, e2, type);
    pp_7_sys_contour->SetLineColor(kCyan + 3);
    pp_7_sys_contour->SetLineWidth(2);
    pp_7_sys_contour->SetLineStyle(1);
    pp_7_sys_contour->SetFillStyle(0);
    pp_7_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(pp_7_sys_contour, "pp 1D pT");

    // Create legend in upper-left (adjust coordinates if needed)
    TLegend *leg1 = new TLegend(0.2, 0.65, 0.45, 0.9); // x1,y1,x2,y2 in NDC
    leg1->SetFillStyle(0);                             // Transparent background
    leg1->SetBorderSize(0);                            // No border
    leg1->SetTextFont(42);                             // CMS-like font
    leg1->SetTextSize(0.035);                          // Adjust for readability

    // Add your systematic contour entry
    leg1->AddEntry(PbPb_1_sys_contour, "TnP Up", "l");
    leg1->AddEntry(PbPb_2_sys_contour, "TnP Down", "l");
    leg1->AddEntry(PbPb_3_sys_contour, "Aco Up", "l");
    leg1->AddEntry(PbPb_4_sys_contour, "Aco Down", "l");
    leg1->AddEntry(PbPb_5_sys_contour, "Bk removed", "l");
    leg1->AddEntry(PbPb_6_sys_contour, "Mass range", "l");
    leg1->AddEntry(PbPb_7_sys_contour, "HF Up", "l");
    leg1->AddEntry(PbPb_8_sys_contour, "HF Down", "l");
    leg1->AddEntry(pp_7_sys_contour, "1D pT", "l");

    // Draw legend
    leg1->Draw();

    // Draw systematics

    std::vector<std::pair<TEllipse *, TEllipse *>> pbpb_pairs = {
        {PbPb_1_sys_contour, PbPb_2_sys_contour}, // TnP
        {PbPb_3_sys_contour, PbPb_4_sys_contour}  // Aco
    };

    std::vector<TEllipse *> pbpb_singles = {
        PbPb_5_sys_contour, // nobk
        PbPb_6_sys_contour, // massrange
        PbPb_7_sys_contour, // HFup
        PbPb_8_sys_contour  // HFdown
    };

    std::vector<std::pair<TEllipse *, TEllipse *>> pp_pairs = {
        {pp_1_sys_contour, pp_2_sys_contour}, // TnP
        {pp_3_sys_contour, pp_4_sys_contour}  // Aco
    };

    std::vector<TEllipse *> pp_singles = {
        pp_5_sys_contour, // nobk
        pp_6_sys_contour, // massrange
        pp_7_sys_contour  // 1D_pT
    };

    TFile *saveTEllipse;
    if (type == "normal")
    {
        saveTEllipse = new TFile("./contourtest/contour.root", "UPDATE");
    }
    else
    {
        saveTEllipse = new TFile("./contourtest/degen/contour.root", "UPDATE");
    }

    saveTEllipse->cd();

    TEllipse *pbpb_envelope = CreateCombinedEnvelope(v1, pbpb_pairs, pbpb_singles);
    pbpb_envelope->SetLineColor(kBlack);
    pbpb_envelope->SetLineStyle(2);
    pbpb_envelope->SetLineWidth(1);
    pbpb_envelope->Draw("SAME");

    pbpb_envelope->Write(name + "_total_PbPb", 2);
    e1->Write(name + "_stat_PbPb", 2);

    TEllipse *pp_envelope = CreateCombinedEnvelope(v2, pp_pairs, pp_singles);
    pp_envelope->SetLineColor(kBlue);
    pp_envelope->SetLineStyle(2);
    pp_envelope->SetLineWidth(1);
    pp_envelope->Draw("SAME");

    pp_envelope->Write(name + "_total_pp", 2);
    e2->Write(name + "_stat_pp", 2);

    TEllipse *pbpb_syst_only = CombiningSystematicOnly(pbpb_pairs, pbpb_singles);
    pbpb_syst_only->SetLineColor(kBlack);
    pbpb_syst_only->SetLineStyle(2);
    pbpb_syst_only->SetLineWidth(1);
    pbpb_syst_only->Draw("SAME");

    PrintXYErrorsFromEllipse(pbpb_syst_only, "PbPb total syst");

    pbpb_syst_only->Write(name + "_syst_PbPb", 2);
    e3->Write(name + "_stat_PbPb_sub_pp");

    TEllipse *pp_syst_only = CombiningSystematicOnly(pp_pairs, pp_singles);
    pp_syst_only->SetLineColor(kBlue);
    pp_syst_only->SetLineStyle(2);
    pp_syst_only->SetLineWidth(1);
    pp_syst_only->Draw("SAME");

    pp_syst_only->Write(name + "_syst_pp", 2);

    PrintXYErrorsFromEllipse(pp_syst_only, "pp total syst");

    // Here's doing combined version

    TVector2 PbPb_sub_pp_sys_1 = sys_1 - pp_sys_1;
    TVector2 PbPb_sub_pp_sys_2 = sys_2 - pp_sys_2;
    TVector2 PbPb_sub_pp_sys_3 = sys_3 - pp_sys_3;
    TVector2 PbPb_sub_pp_sys_4 = sys_4 - pp_sys_4;
    TVector2 PbPb_sub_pp_sys_5 = sys_5 - pp_sys_5;
    TVector2 PbPb_sub_pp_sys_6 = sys_6 - pp_sys_6;
    TVector2 v2_xy(v2[0], v2[1]);
    TVector2 PbPb_sub_pp_sys_7 = sys_7 - v2_xy;
    TVector2 PbPb_sub_pp_sys_8 = sys_8 - v2_xy;
    TVector2 v1_xy(v1[0], v1[1]);
    TVector2 PbPb_sub_pp_sys_9 = v1_xy - pp_sys_7;

    // Combined Systematic

    // tnpU

    TMarker *PbPb_sub_pp_1 = new TMarker(PbPb_sub_pp_sys_1.X(), PbPb_sub_pp_sys_1.Y(), 20);
    PbPb_sub_pp_1->SetMarkerColor(kRed);
    PbPb_sub_pp_1->SetMarkerSize(1.2);
    PbPb_sub_pp_1->Draw("same");

    TEllipse *PbPb_sub_pp_1_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_1, e3, type);
    PbPb_sub_pp_1_sys_contour->SetLineColor(kRed);
    PbPb_sub_pp_1_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_1_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_1_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_1_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_1_sys_contour, "PbPb - pp TnP U");

    // tnpD

    TMarker *PbPb_sub_pp_2 = new TMarker(PbPb_sub_pp_sys_2.X(), PbPb_sub_pp_sys_2.Y(), 20);
    PbPb_sub_pp_2->SetMarkerColor(kBlue);
    PbPb_sub_pp_2->SetMarkerSize(1.2);
    PbPb_sub_pp_2->Draw("same");

    TEllipse *PbPb_sub_pp_2_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_2, e3, type);
    PbPb_sub_pp_2_sys_contour->SetLineColor(kBlue);
    PbPb_sub_pp_2_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_2_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_2_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_2_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_2_sys_contour, "PbPb - pp TnP D");

    // AcoUp

    TMarker *PbPb_sub_pp_3 = new TMarker(PbPb_sub_pp_sys_3.X(), PbPb_sub_pp_sys_3.Y(), 20);
    PbPb_sub_pp_3->SetMarkerColor(kGreen + 2);
    PbPb_sub_pp_3->SetMarkerSize(1.2);
    PbPb_sub_pp_3->Draw("same");

    TEllipse *PbPb_sub_pp_3_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_3, e3, type);
    PbPb_sub_pp_3_sys_contour->SetLineColor(kGreen + 2);
    PbPb_sub_pp_3_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_3_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_3_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_3_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_3_sys_contour, "PbPb - pp Aco Up");

    // AcoDown

    TMarker *PbPb_sub_pp_4 = new TMarker(PbPb_sub_pp_sys_4.X(), PbPb_sub_pp_sys_4.Y(), 20);
    PbPb_sub_pp_4->SetMarkerColor(kMagenta);
    PbPb_sub_pp_4->SetMarkerSize(1.2);
    PbPb_sub_pp_4->Draw("same");

    TEllipse *PbPb_sub_pp_4_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_4, e3, type);
    PbPb_sub_pp_4_sys_contour->SetLineColor(kMagenta);
    PbPb_sub_pp_4_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_4_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_4_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_4_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_4_sys_contour, "PbPb - pp Aco Down");

    // nobk

    TMarker *PbPb_sub_pp_5 = new TMarker(PbPb_sub_pp_sys_5.X(), PbPb_sub_pp_sys_5.Y(), 20);
    PbPb_sub_pp_5->SetMarkerColor(kGray + 2);
    PbPb_sub_pp_5->SetMarkerSize(1.2);
    PbPb_sub_pp_5->Draw("same");

    TEllipse *PbPb_sub_pp_5_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_5, e3, type);
    PbPb_sub_pp_5_sys_contour->SetLineColor(kGray + 2);
    PbPb_sub_pp_5_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_5_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_5_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_5_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_5_sys_contour, "PbPb - pp no bk");

    // massrange

    TMarker *PbPb_sub_pp_6 = new TMarker(PbPb_sub_pp_sys_6.X(), PbPb_sub_pp_sys_6.Y(), 20);
    PbPb_sub_pp_6->SetMarkerColor(kOrange + 7);
    PbPb_sub_pp_6->SetMarkerSize(1.2);
    PbPb_sub_pp_6->Draw("same");

    TEllipse *PbPb_sub_pp_6_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_6, e3, type);
    PbPb_sub_pp_6_sys_contour->SetLineColor(kOrange + 7);
    PbPb_sub_pp_6_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_6_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_6_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_6_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_6_sys_contour, "PbPb - pp mass range");

    // HFup

    TMarker *PbPb_sub_pp_7 = new TMarker(PbPb_sub_pp_sys_7.X(), PbPb_sub_pp_sys_7.Y(), 20);
    PbPb_sub_pp_7->SetMarkerColor(kCyan);
    PbPb_sub_pp_7->SetMarkerSize(1.2);
    PbPb_sub_pp_7->Draw("same");

    TEllipse *PbPb_sub_pp_7_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_7, e3, type);
    PbPb_sub_pp_7_sys_contour->SetLineColor(kCyan);
    PbPb_sub_pp_7_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_7_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_7_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_7_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_7_sys_contour, "PbPb - pp HFup");

    // HFdown

    TMarker *PbPb_sub_pp_8 = new TMarker(PbPb_sub_pp_sys_8.X(), PbPb_sub_pp_sys_8.Y(), 20);
    PbPb_sub_pp_8->SetMarkerColor(kBlack);
    PbPb_sub_pp_8->SetMarkerSize(1.2);
    PbPb_sub_pp_8->Draw("same");

    TEllipse *PbPb_sub_pp_8_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_8, e3, type);
    PbPb_sub_pp_8_sys_contour->SetLineColor(kBlack);
    PbPb_sub_pp_8_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_8_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_8_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_8_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_8_sys_contour, "PbPb - pp HFdown");

    // 1D_pT

    TMarker *PbPb_sub_pp_9 = new TMarker(PbPb_sub_pp_sys_9.X(), PbPb_sub_pp_sys_9.Y(), 20);
    PbPb_sub_pp_9->SetMarkerColor(kCyan + 3);
    PbPb_sub_pp_9->SetMarkerSize(1.2);
    PbPb_sub_pp_9->Draw("same");

    TEllipse *PbPb_sub_pp_9_sys_contour = CreateScaledEllipseFromSystematic(PbPb_sub_pp_sys_9, e3, type);
    PbPb_sub_pp_9_sys_contour->SetLineColor(kCyan + 3);
    PbPb_sub_pp_9_sys_contour->SetLineWidth(2);
    PbPb_sub_pp_9_sys_contour->SetLineStyle(1);
    PbPb_sub_pp_9_sys_contour->SetFillStyle(0);
    PbPb_sub_pp_9_sys_contour->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_9_sys_contour, "PbPb - pp 1D_pT");

    std::vector<std::pair<TEllipse *, TEllipse *>> PbPb_sub_pp_pairs = {
        {PbPb_sub_pp_1_sys_contour, PbPb_sub_pp_2_sys_contour}, // TnP
        {PbPb_sub_pp_3_sys_contour, PbPb_sub_pp_4_sys_contour}  // Aco
    };

    std::vector<TEllipse *> PbPb_sub_pp_singles = {
        PbPb_sub_pp_5_sys_contour, // nobk
        PbPb_sub_pp_6_sys_contour, // massrange
        PbPb_sub_pp_7_sys_contour, // HFup
        PbPb_sub_pp_8_sys_contour, // HFdown
        PbPb_sub_pp_9_sys_contour  // 1D_pT
    };

    TEllipse *PbPb_sub_pp_envelope = CreateCombinedEnvelope(EllipseToVector(e3), PbPb_sub_pp_pairs, PbPb_sub_pp_singles);
    PbPb_sub_pp_envelope->SetLineColor(kRed);
    PbPb_sub_pp_envelope->SetLineStyle(2);
    PbPb_sub_pp_envelope->SetLineWidth(2);
    PbPb_sub_pp_envelope->Draw("SAME");

    PbPb_sub_pp_envelope->Write(name + "_total_PbPb_sub_pp", 2);

    TEllipse *PbPb_sub_pp_syst_only = CombiningSystematicOnly(PbPb_sub_pp_pairs, PbPb_sub_pp_singles);
    PbPb_sub_pp_syst_only->SetLineColor(kRed);
    PbPb_sub_pp_syst_only->SetLineStyle(2);
    PbPb_sub_pp_syst_only->SetLineWidth(2);
    PbPb_sub_pp_syst_only->Draw("SAME");

    PrintXYErrorsFromEllipse(PbPb_sub_pp_syst_only, "PbPb - pp total syst");

    PbPb_sub_pp_syst_only->Write(name + "_syst_PbPb_sub_pp", 2);

    /*TVectorD pbpb_total = EllipseToVector(pbpb_envelope);
    TVectorD pp_total = EllipseToVector(pp_envelope);

    TEllipse *total_contour_subtracted = CombineEllipsesFromVectors(pbpb_total, pp_total);
    total_contour_subtracted->SetLineColor(kRed);
    total_contour_subtracted->SetLineWidth(2);
    total_contour_subtracted->SetLineStyle(2);
    total_contour_subtracted->SetFillStyle(0);
    total_contour_subtracted->Draw("SAME");*/

    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->SetFillStyle(0);  // Transparent fill
    leg->SetBorderSize(0); // No border box
    leg->SetLineColor(0);  // No outline line
    leg->SetTextFont(42);  // CMS-style font
    leg->AddEntry(e1, "PbPb", "l");
    leg->AddEntry(e2, "pp", "l");
    leg->AddEntry(e3, "PbPb - pp", "l");
    leg->Draw();

    TLatex latex;
    latex.SetTextSize(0.02);
    latex.SetTextAlign(13); // left align, top align

    double ystart = 0.3;
    double dy = 0.035;

    // latex.DrawLatexNDC(0.2, ystart, Form("Ellipse 1: x=%.3f, y=%.3f, a=%.3f, b=%.3f, #theta=%.1f^{#circ}", v1[0], v1[1], v1[2], v1[3], v1[4]));
    // latex.DrawLatexNDC(0.2, ystart - dy, Form("Ellipse 2: x=%.3f, y=%.3f, a=%.3f, b=%.3f, #theta=%.1f^{#circ}", v2[0], v2[1], v2[2], v2[3], v2[4]));
    // latex.DrawLatexNDC(0.2, ystart - 2 * dy, Form("Combined: x=%.3f, y=%.3f, a=%.3f, b=%.3f, #theta=%.1f^{#circ}", e3->GetX1(), e3->GetY1(), e3->GetR1(), e3->GetR2(), e3->GetTheta()));

    c1->Update();
    c1->SaveAs(Savename);

    saveTEllipse->Close();

    // Save all result
}

void get_combined_contour(TString type = "degen")
{

    TFile *f1 = new TFile("./contourrootfile/everything.root", "READ");
    TFile *f2 = new TFile("./bestfittemplaterootfile/template.root", "READ");

    TString savepath = "";

    if (type == "normal")
    {
        savepath = "./contourtest/";
    }
    else
    {
        savepath = "./contourtest/degen/";
    }

    // Place to read Noimnal

    TVector2 *PbPb_0 = (TVector2 *)f2->Get("PbPb_local_min_cent_0");
    TVector2 *PbPb_1 = (TVector2 *)f2->Get("PbPb_local_min_cent_1");
    TVector2 *PbPb_2 = (TVector2 *)f2->Get("PbPb_local_min_cent_2");
    TVector2 *PbPb_3 = (TVector2 *)f2->Get("PbPb_local_min_cent_3");
    TVector2 *PbPb_4 = (TVector2 *)f2->Get("PbPb_local_min_cent_4");
    TVector2 *PbPb_5 = (TVector2 *)f2->Get("PbPb_local_min_cent_5");
    TVector2 *PbPb_6 = (TVector2 *)f2->Get("PbPb_local_min_cent_6");
    TVector2 *PbPb_7 = (TVector2 *)f2->Get("PbPb_local_min_cent_7");
    TVector2 *PbPb_8 = (TVector2 *)f2->Get("PbPb_local_min_cent_8");
    TVector2 *PbPb_9 = (TVector2 *)f2->Get("PbPb_local_min_cent_9");
    TVector2 *PbPb_10 = (TVector2 *)f2->Get("PbPb_local_min_cent_10");

    TVector2 *pp_0 = (TVector2 *)f2->Get("pp_local_min_period_22");

    // Place to read systematic points
    TVector2 *pp_0_tnpU = (TVector2 *)f2->Get("pp_local_min_tnpU_period_22");
    TVector2 *pp_0_tnpD = (TVector2 *)f2->Get("pp_local_min_tnpD_period_22");
    TVector2 *pp_0_AcoUp = (TVector2 *)f2->Get("pp_local_min_AcoUp_period_22");
    TVector2 *pp_0_AcoDown = (TVector2 *)f2->Get("pp_local_min_AcoDown_period_22");
    TVector2 *pp_0_nobk = (TVector2 *)f2->Get("pp_local_min_no_bk_period_22");
    TVector2 *pp_0_massrange = (TVector2 *)f2->Get("pp_local_min_massrange_period_22");
    TVector2 *pp_0_1D_pT = (TVector2 *)f2->Get("pp_local_min_1D_pT_period_22");

    // Here need to add pT 2D and pT 1D

    TVector2 *PbPb_0_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_0");
    TVector2 *PbPb_1_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_1");
    TVector2 *PbPb_2_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_2");
    TVector2 *PbPb_3_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_3");
    TVector2 *PbPb_4_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_4");
    TVector2 *PbPb_5_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_5");
    TVector2 *PbPb_6_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_6");
    TVector2 *PbPb_7_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_7");
    TVector2 *PbPb_8_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_8");
    TVector2 *PbPb_9_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_9");
    TVector2 *PbPb_10_tnpU = (TVector2 *)f2->Get("PbPb_local_min_tnpU_cent_10");

    TVector2 *PbPb_0_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_0");
    TVector2 *PbPb_1_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_1");
    TVector2 *PbPb_2_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_2");
    TVector2 *PbPb_3_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_3");
    TVector2 *PbPb_4_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_4");
    TVector2 *PbPb_5_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_5");
    TVector2 *PbPb_6_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_6");
    TVector2 *PbPb_7_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_7");
    TVector2 *PbPb_8_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_8");
    TVector2 *PbPb_9_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_9");
    TVector2 *PbPb_10_tnpD = (TVector2 *)f2->Get("PbPb_local_min_tnpD_cent_10");

    TVector2 *PbPb_0_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_0");
    TVector2 *PbPb_1_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_1");
    TVector2 *PbPb_2_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_2");
    TVector2 *PbPb_3_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_3");
    TVector2 *PbPb_4_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_4");
    TVector2 *PbPb_5_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_5");
    TVector2 *PbPb_6_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_6");
    TVector2 *PbPb_7_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_7");
    TVector2 *PbPb_8_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_8");
    TVector2 *PbPb_9_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_9");
    TVector2 *PbPb_10_AcoUp = (TVector2 *)f2->Get("PbPb_local_min_AcoUp_cent_10");

    TVector2 *PbPb_0_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_0");
    TVector2 *PbPb_1_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_1");
    TVector2 *PbPb_2_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_2");
    TVector2 *PbPb_3_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_3");
    TVector2 *PbPb_4_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_4");
    TVector2 *PbPb_5_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_5");
    TVector2 *PbPb_6_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_6");
    TVector2 *PbPb_7_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_7");
    TVector2 *PbPb_8_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_8");
    TVector2 *PbPb_9_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_9");
    TVector2 *PbPb_10_AcoDown = (TVector2 *)f2->Get("PbPb_local_min_AcoDown_cent_10");

    TVector2 *PbPb_0_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_0");
    TVector2 *PbPb_1_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_1");
    TVector2 *PbPb_2_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_2");
    TVector2 *PbPb_3_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_3");
    TVector2 *PbPb_4_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_4");
    TVector2 *PbPb_5_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_5");
    TVector2 *PbPb_6_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_6");
    TVector2 *PbPb_7_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_7");
    TVector2 *PbPb_8_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_8");
    TVector2 *PbPb_9_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_9");
    TVector2 *PbPb_10_nobk = (TVector2 *)f2->Get("PbPb_local_min_no_bk_cent_10");

    TVector2 *PbPb_0_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_0");
    TVector2 *PbPb_1_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_1");
    TVector2 *PbPb_2_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_2");
    TVector2 *PbPb_3_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_3");
    TVector2 *PbPb_4_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_4");
    TVector2 *PbPb_5_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_5");
    TVector2 *PbPb_6_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_6");
    TVector2 *PbPb_7_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_7");
    TVector2 *PbPb_8_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_8");
    TVector2 *PbPb_9_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_9");
    TVector2 *PbPb_10_massrange = (TVector2 *)f2->Get("PbPb_local_min_massrange_cent_10");

    TVector2 *PbPb_0_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_0");
    TVector2 *PbPb_1_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_1");
    TVector2 *PbPb_2_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_2");
    TVector2 *PbPb_3_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_3");
    TVector2 *PbPb_4_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_4");
    TVector2 *PbPb_5_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_5");
    TVector2 *PbPb_6_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_6");
    TVector2 *PbPb_7_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_7");
    TVector2 *PbPb_8_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_8");
    TVector2 *PbPb_9_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_9");
    TVector2 *PbPb_10_HFup = (TVector2 *)f2->Get("PbPb_local_min_HFup_cent_10");

    TVector2 *PbPb_0_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_0");
    TVector2 *PbPb_1_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_1");
    TVector2 *PbPb_2_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_2");
    TVector2 *PbPb_3_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_3");
    TVector2 *PbPb_4_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_4");
    TVector2 *PbPb_5_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_5");
    TVector2 *PbPb_6_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_6");
    TVector2 *PbPb_7_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_7");
    TVector2 *PbPb_8_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_8");
    TVector2 *PbPb_9_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_9");
    TVector2 *PbPb_10_HFdown = (TVector2 *)f2->Get("PbPb_local_min_HFdown_cent_10");

    TVectorD pp_sig1 = doeverything("pp, |#eta| < 2.4, Nominal, Period: (22)", "./contourtest/pp.png", f1, pp_0, 1);
    TVectorD pp_sig2 = doeverything("pp, |#eta| < 2.4, Nominal, Period: (22)", "./contourtest/pp.png", f1, pp_0, 2);

    TVectorD PbPb_sig_1_0 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-10.0)", "./contourtest/PbPb_0.png", f1, PbPb_0, 1);
    TVectorD PbPb_sig_2_0 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-10.0)", "./contourtest/PbPb_0.png", f1, PbPb_0, 2);

    TVectorD PbPb_sig_1_1 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (10.0-20.0)", "./contourtest/PbPb_1.png", f1, PbPb_1, 1);
    TVectorD PbPb_sig_2_1 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (10.0-20.0)", "./contourtest/PbPb_1.png", f1, PbPb_1, 2);

    TVectorD PbPb_sig_1_2 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (20.0-30.0)", "./contourtest/PbPb_2.png", f1, PbPb_2, 1);
    TVectorD PbPb_sig_2_2 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (20.0-30.0)", "./contourtest/PbPb_2.png", f1, PbPb_2, 2);

    TVectorD PbPb_sig_1_3 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (30.0-100.0)", "./contourtest/PbPb_3.png", f1, PbPb_3, 1);
    TVectorD PbPb_sig_2_3 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (30.0-100.0)", "./contourtest/PbPb_3.png", f1, PbPb_3, 2);

    TVectorD PbPb_sig_1_4 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-100.0)", "./contourtest/PbPb_4.png", f1, PbPb_4, 1);
    TVectorD PbPb_sig_2_4 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-100.0)", "./contourtest/PbPb_4.png", f1, PbPb_4, 2);

    TVectorD PbPb_sig_1_5 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-6.5)", "./contourtest/PbPb_5.png", f1, PbPb_5, 1);
    TVectorD PbPb_sig_2_5 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-6.5)", "./contourtest/PbPb_5.png", f1, PbPb_5, 2);

    TVectorD PbPb_sig_1_6 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (6.5-14.5)", "./contourtest/PbPb_6.png", f1, PbPb_6, 1);
    TVectorD PbPb_sig_2_6 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (6.5-14.5)", "./contourtest/PbPb_6.png", f1, PbPb_6, 2);

    TVectorD PbPb_sig_1_7 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (14.5-26.5)", "./contourtest/PbPb_7.png", f1, PbPb_7, 1);
    TVectorD PbPb_sig_2_7 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (14.5-26.5)", "./contourtest/PbPb_7.png", f1, PbPb_7, 2);

    TVectorD PbPb_sig_1_8 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (26.5-100.0)", "./contourtest/PbPb_8.png", f1, PbPb_8, 1);
    TVectorD PbPb_sig_2_8 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (26.5-100.0)", "./contourtest/PbPb_8.png", f1, PbPb_8, 2);

    TVectorD PbPb_sig_1_9 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-15.0)", "./contourtest/PbPb_9.png", f1, PbPb_9, 1);
    TVectorD PbPb_sig_2_9 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (0.0-15.0)", "./contourtest/PbPb_9.png", f1, PbPb_9, 2);

    TVectorD PbPb_sig_1_10 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (15.0-100.0)", "./contourtest/PbPb_10.png", f1, PbPb_10, 1);
    TVectorD PbPb_sig_2_10 = doeverything("PbPb, |#eta| < 2.4, Nominal, centrality: (15.0-100.0)", "./contourtest/PbPb_10.png", f1, PbPb_10, 2);

    setTDRStyle();
    gStyle->SetFrameLineWidth(3);

    DrawEllipsesComparison(PbPb_sig_1_0, pp_sig1, savepath + "PbPb_0_sig_1.png", *PbPb_0_tnpU, *PbPb_0_tnpD, *PbPb_0_AcoUp, *PbPb_0_AcoDown, *PbPb_0_nobk, *PbPb_0_massrange, *PbPb_0_HFup, *PbPb_0_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_0_10");
    DrawEllipsesComparison(PbPb_sig_2_0, pp_sig2, savepath + "PbPb_0_sig_2.png", *PbPb_0_tnpU, *PbPb_0_tnpD, *PbPb_0_AcoUp, *PbPb_0_AcoDown, *PbPb_0_nobk, *PbPb_0_massrange, *PbPb_0_HFup, *PbPb_0_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_0_10");

    DrawEllipsesComparison(PbPb_sig_1_1, pp_sig1, savepath + "PbPb_1_sig_1.png", *PbPb_1_tnpU, *PbPb_1_tnpD, *PbPb_1_AcoUp, *PbPb_1_AcoDown, *PbPb_1_nobk, *PbPb_1_massrange, *PbPb_1_HFup, *PbPb_1_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_10_20");
    DrawEllipsesComparison(PbPb_sig_2_1, pp_sig2, savepath + "PbPb_1_sig_2.png", *PbPb_1_tnpU, *PbPb_1_tnpD, *PbPb_1_AcoUp, *PbPb_1_AcoDown, *PbPb_1_nobk, *PbPb_1_massrange, *PbPb_1_HFup, *PbPb_1_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_10_20");

    DrawEllipsesComparison(PbPb_sig_1_2, pp_sig1, savepath + "PbPb_2_sig_1.png", *PbPb_2_tnpU, *PbPb_2_tnpD, *PbPb_2_AcoUp, *PbPb_2_AcoDown, *PbPb_2_nobk, *PbPb_2_massrange, *PbPb_2_HFup, *PbPb_2_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_20_30");
    DrawEllipsesComparison(PbPb_sig_2_2, pp_sig2, savepath + "PbPb_2_sig_2.png", *PbPb_2_tnpU, *PbPb_2_tnpD, *PbPb_2_AcoUp, *PbPb_2_AcoDown, *PbPb_2_nobk, *PbPb_2_massrange, *PbPb_2_HFup, *PbPb_2_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_20_30");

    DrawEllipsesComparison(PbPb_sig_1_3, pp_sig1, savepath + "PbPb_3_sig_1.png", *PbPb_3_tnpU, *PbPb_3_tnpD, *PbPb_3_AcoUp, *PbPb_3_AcoDown, *PbPb_3_nobk, *PbPb_3_massrange, *PbPb_3_HFup, *PbPb_3_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_30_100");
    DrawEllipsesComparison(PbPb_sig_2_3, pp_sig2, savepath + "PbPb_3_sig_2.png", *PbPb_3_tnpU, *PbPb_3_tnpD, *PbPb_3_AcoUp, *PbPb_3_AcoDown, *PbPb_3_nobk, *PbPb_3_massrange, *PbPb_3_HFup, *PbPb_3_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_30_100");

    DrawEllipsesComparison(PbPb_sig_1_4, pp_sig1, savepath + "PbPb_4_sig_1.png", *PbPb_4_tnpU, *PbPb_4_tnpD, *PbPb_4_AcoUp, *PbPb_4_AcoDown, *PbPb_4_nobk, *PbPb_4_massrange, *PbPb_4_HFup, *PbPb_4_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_0_100");
    DrawEllipsesComparison(PbPb_sig_2_4, pp_sig2, savepath + "PbPb_4_sig_2.png", *PbPb_4_tnpU, *PbPb_4_tnpD, *PbPb_4_AcoUp, *PbPb_4_AcoDown, *PbPb_4_nobk, *PbPb_4_massrange, *PbPb_4_HFup, *PbPb_4_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_0_100");

    DrawEllipsesComparison(PbPb_sig_1_5, pp_sig1, savepath + "PbPb_5_sig_1.png", *PbPb_5_tnpU, *PbPb_5_tnpD, *PbPb_5_AcoUp, *PbPb_5_AcoDown, *PbPb_5_nobk, *PbPb_5_massrange, *PbPb_5_HFup, *PbPb_5_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_0_6.5");
    DrawEllipsesComparison(PbPb_sig_2_5, pp_sig2, savepath + "PbPb_5_sig_2.png", *PbPb_5_tnpU, *PbPb_5_tnpD, *PbPb_5_AcoUp, *PbPb_5_AcoDown, *PbPb_5_nobk, *PbPb_5_massrange, *PbPb_5_HFup, *PbPb_5_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_0_6.5");

    DrawEllipsesComparison(PbPb_sig_1_6, pp_sig1, savepath + "PbPb_6_sig_1.png", *PbPb_6_tnpU, *PbPb_6_tnpD, *PbPb_6_AcoUp, *PbPb_6_AcoDown, *PbPb_6_nobk, *PbPb_6_massrange, *PbPb_6_HFup, *PbPb_6_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_6.5_14.5");
    DrawEllipsesComparison(PbPb_sig_2_6, pp_sig2, savepath + "PbPb_6_sig_2.png", *PbPb_6_tnpU, *PbPb_6_tnpD, *PbPb_6_AcoUp, *PbPb_6_AcoDown, *PbPb_6_nobk, *PbPb_6_massrange, *PbPb_6_HFup, *PbPb_6_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_6.5_14.5");

    DrawEllipsesComparison(PbPb_sig_1_7, pp_sig1, savepath + "PbPb_7_sig_1.png", *PbPb_7_tnpU, *PbPb_7_tnpD, *PbPb_7_AcoUp, *PbPb_7_AcoDown, *PbPb_7_nobk, *PbPb_7_massrange, *PbPb_7_HFup, *PbPb_7_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_14.5_26.5");
    DrawEllipsesComparison(PbPb_sig_2_7, pp_sig2, savepath + "PbPb_7_sig_2.png", *PbPb_7_tnpU, *PbPb_7_tnpD, *PbPb_7_AcoUp, *PbPb_7_AcoDown, *PbPb_7_nobk, *PbPb_7_massrange, *PbPb_7_HFup, *PbPb_7_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_14.5_26.5");

    DrawEllipsesComparison(PbPb_sig_1_8, pp_sig1, savepath + "PbPb_8_sig_1.png", *PbPb_8_tnpU, *PbPb_8_tnpD, *PbPb_8_AcoUp, *PbPb_8_AcoDown, *PbPb_8_nobk, *PbPb_8_massrange, *PbPb_8_HFup, *PbPb_8_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_26.5_100");
    DrawEllipsesComparison(PbPb_sig_2_8, pp_sig2, savepath + "PbPb_8_sig_2.png", *PbPb_8_tnpU, *PbPb_8_tnpD, *PbPb_8_AcoUp, *PbPb_8_AcoDown, *PbPb_8_nobk, *PbPb_8_massrange, *PbPb_8_HFup, *PbPb_8_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_26.5_100");

    DrawEllipsesComparison(PbPb_sig_1_9, pp_sig1, savepath + "PbPb_9_sig_1.png", *PbPb_9_tnpU, *PbPb_9_tnpD, *PbPb_9_AcoUp, *PbPb_9_AcoDown, *PbPb_9_nobk, *PbPb_9_massrange, *PbPb_9_HFup, *PbPb_9_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_0_15");
    DrawEllipsesComparison(PbPb_sig_2_9, pp_sig2, savepath + "PbPb_9_sig_2.png", *PbPb_9_tnpU, *PbPb_9_tnpD, *PbPb_9_AcoUp, *PbPb_9_AcoDown, *PbPb_9_nobk, *PbPb_9_massrange, *PbPb_9_HFup, *PbPb_9_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_0_15");

    DrawEllipsesComparison(PbPb_sig_1_10, pp_sig1, savepath + "PbPb_10_sig_1.png", *PbPb_10_tnpU, *PbPb_10_tnpD, *PbPb_10_AcoUp, *PbPb_10_AcoDown, *PbPb_10_nobk, *PbPb_10_massrange, *PbPb_10_HFup, *PbPb_10_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_1_15_100");
    DrawEllipsesComparison(PbPb_sig_2_10, pp_sig2, savepath + "PbPb_10_sig_2.png", *PbPb_10_tnpU, *PbPb_10_tnpD, *PbPb_10_AcoUp, *PbPb_10_AcoDown, *PbPb_10_nobk, *PbPb_10_massrange, *PbPb_10_HFup, *PbPb_10_HFdown, *pp_0_tnpU, *pp_0_tnpD, *pp_0_AcoUp, *pp_0_AcoDown, *pp_0_nobk, *pp_0_massrange, *pp_0_1D_pT, type, "sig_2_15_100");

    std::vector<TVectorD> pbpbEll5 = {
        PbPb_sig_1_4, PbPb_sig_1_0, PbPb_sig_1_1, PbPb_sig_1_2, PbPb_sig_1_3};

    std::vector<TVectorD> ppEll5 = {
        pp_sig1, pp_sig1, pp_sig1, pp_sig1, pp_sig1};

    auto EllipseToVector = [](TEllipse *e)
    {
        TVectorD v(5);
        v[0] = e->GetX1();
        v[1] = e->GetY1();
        v[2] = e->GetR1();
        v[3] = e->GetR2();
        v[4] = e->GetTheta();
        return v;
    };

    std::vector<TVectorD> subEll5;
    {
        TEllipse *e0 = CombineEllipsesFromVectors(PbPb_sig_1_4, pp_sig1);
        TEllipse *e1 = CombineEllipsesFromVectors(PbPb_sig_1_0, pp_sig1);
        TEllipse *e2 = CombineEllipsesFromVectors(PbPb_sig_1_1, pp_sig1);
        TEllipse *e3 = CombineEllipsesFromVectors(PbPb_sig_1_2, pp_sig1);
        TEllipse *e4 = CombineEllipsesFromVectors(PbPb_sig_1_3, pp_sig1);

        subEll5 = {
            EllipseToVector(e0),
            EllipseToVector(e1),
            EllipseToVector(e2),
            EllipseToVector(e3),
            EllipseToVector(e4)};

        delete e0;
        delete e1;
        delete e2;
        delete e3;
        delete e4;
    }

    std::vector<TVector2 *> pbpb_nominal5 = {PbPb_4, PbPb_0, PbPb_1, PbPb_2, PbPb_3};
    std::vector<TVector2 *> pbpb_tnpU5 = {PbPb_4_tnpU, PbPb_0_tnpU, PbPb_1_tnpU, PbPb_2_tnpU, PbPb_3_tnpU};
    std::vector<TVector2 *> pbpb_tnpD5 = {PbPb_4_tnpD, PbPb_0_tnpD, PbPb_1_tnpD, PbPb_2_tnpD, PbPb_3_tnpD};
    std::vector<TVector2 *> pbpb_acoup5 = {PbPb_4_AcoUp, PbPb_0_AcoUp, PbPb_1_AcoUp, PbPb_2_AcoUp, PbPb_3_AcoUp};
    std::vector<TVector2 *> pbpb_acodown5 = {PbPb_4_AcoDown, PbPb_0_AcoDown, PbPb_1_AcoDown, PbPb_2_AcoDown, PbPb_3_AcoDown};
    std::vector<TVector2 *> pbpb_nobk5 = {PbPb_4_nobk, PbPb_0_nobk, PbPb_1_nobk, PbPb_2_nobk, PbPb_3_nobk};
    std::vector<TVector2 *> pbpb_massrange5 = {PbPb_4_massrange, PbPb_0_massrange, PbPb_1_massrange, PbPb_2_massrange, PbPb_3_massrange};
    std::vector<TVector2 *> pbpb_HFup5 = {PbPb_4_HFup, PbPb_0_HFup, PbPb_1_HFup, PbPb_2_HFup, PbPb_3_HFup};
    std::vector<TVector2 *> pbpb_HFdown5 = {PbPb_4_HFdown, PbPb_0_HFdown, PbPb_1_HFdown, PbPb_2_HFdown, PbPb_3_HFdown};
    std::vector<TVector2 *> pbpb_1DpT5 = {PbPb_4, PbPb_0, PbPb_1, PbPb_2, PbPb_3}; // no PbPb 1D_pT in your current inputs, so keep nominal or replace later

    std::vector<TVector2 *> pp_nominal5 = {pp_0, pp_0, pp_0, pp_0, pp_0};
    std::vector<TVector2 *> pp_tnpU5 = {pp_0_tnpU, pp_0_tnpU, pp_0_tnpU, pp_0_tnpU, pp_0_tnpU};
    std::vector<TVector2 *> pp_tnpD5 = {pp_0_tnpD, pp_0_tnpD, pp_0_tnpD, pp_0_tnpD, pp_0_tnpD};
    std::vector<TVector2 *> pp_acoup5 = {pp_0_AcoUp, pp_0_AcoUp, pp_0_AcoUp, pp_0_AcoUp, pp_0_AcoUp};
    std::vector<TVector2 *> pp_acodown5 = {pp_0_AcoDown, pp_0_AcoDown, pp_0_AcoDown, pp_0_AcoDown, pp_0_AcoDown};
    std::vector<TVector2 *> pp_nobk5 = {pp_0_nobk, pp_0_nobk, pp_0_nobk, pp_0_nobk, pp_0_nobk};
    std::vector<TVector2 *> pp_massrange5 = {pp_0_massrange, pp_0_massrange, pp_0_massrange, pp_0_massrange, pp_0_massrange};
    std::vector<TVector2 *> pp_HFup5 = {pp_0, pp_0, pp_0, pp_0, pp_0}; // pp has no HF syst in current inputs
    std::vector<TVector2 *> pp_HFdown5 = {pp_0, pp_0, pp_0, pp_0, pp_0};
    std::vector<TVector2 *> pp_1DpT5 = {pp_0_1D_pT, pp_0_1D_pT, pp_0_1D_pT, pp_0_1D_pT, pp_0_1D_pT};

    std::vector<TVector2 *> sub_nominal5;
    std::vector<TVector2 *> sub_tnpU5;
    std::vector<TVector2 *> sub_tnpD5;
    std::vector<TVector2 *> sub_acoup5;
    std::vector<TVector2 *> sub_acodown5;
    std::vector<TVector2 *> sub_nobk5;
    std::vector<TVector2 *> sub_massrange5;
    std::vector<TVector2 *> sub_HFup5;
    std::vector<TVector2 *> sub_HFdown5;
    std::vector<TVector2 *> sub_1DpT5;

    auto makeDiff = [](TVector2 *a, TVector2 *b)
    {
        return new TVector2(a->X() - b->X(), a->Y() - b->Y());
    };

    sub_nominal5 = {makeDiff(PbPb_4, pp_0), makeDiff(PbPb_0, pp_0), makeDiff(PbPb_1, pp_0), makeDiff(PbPb_2, pp_0), makeDiff(PbPb_3, pp_0)};
    sub_tnpU5 = {makeDiff(PbPb_4_tnpU, pp_0_tnpU), makeDiff(PbPb_0_tnpU, pp_0_tnpU), makeDiff(PbPb_1_tnpU, pp_0_tnpU), makeDiff(PbPb_2_tnpU, pp_0_tnpU), makeDiff(PbPb_3_tnpU, pp_0_tnpU)};
    sub_tnpD5 = {makeDiff(PbPb_4_tnpD, pp_0_tnpD), makeDiff(PbPb_0_tnpD, pp_0_tnpD), makeDiff(PbPb_1_tnpD, pp_0_tnpD), makeDiff(PbPb_2_tnpD, pp_0_tnpD), makeDiff(PbPb_3_tnpD, pp_0_tnpD)};
    sub_acoup5 = {makeDiff(PbPb_4_AcoUp, pp_0_AcoUp), makeDiff(PbPb_0_AcoUp, pp_0_AcoUp), makeDiff(PbPb_1_AcoUp, pp_0_AcoUp), makeDiff(PbPb_2_AcoUp, pp_0_AcoUp), makeDiff(PbPb_3_AcoUp, pp_0_AcoUp)};
    sub_acodown5 = {makeDiff(PbPb_4_AcoDown, pp_0_AcoDown), makeDiff(PbPb_0_AcoDown, pp_0_AcoDown), makeDiff(PbPb_1_AcoDown, pp_0_AcoDown), makeDiff(PbPb_2_AcoDown, pp_0_AcoDown), makeDiff(PbPb_3_AcoDown, pp_0_AcoDown)};
    sub_nobk5 = {makeDiff(PbPb_4_nobk, pp_0_nobk), makeDiff(PbPb_0_nobk, pp_0_nobk), makeDiff(PbPb_1_nobk, pp_0_nobk), makeDiff(PbPb_2_nobk, pp_0_nobk), makeDiff(PbPb_3_nobk, pp_0_nobk)};
    sub_massrange5 = {makeDiff(PbPb_4_massrange, pp_0_massrange), makeDiff(PbPb_0_massrange, pp_0_massrange), makeDiff(PbPb_1_massrange, pp_0_massrange), makeDiff(PbPb_2_massrange, pp_0_massrange), makeDiff(PbPb_3_massrange, pp_0_massrange)};
    sub_HFup5 = {makeDiff(PbPb_4_HFup, pp_0), makeDiff(PbPb_0_HFup, pp_0), makeDiff(PbPb_1_HFup, pp_0), makeDiff(PbPb_2_HFup, pp_0), makeDiff(PbPb_3_HFup, pp_0)};
    sub_HFdown5 = {makeDiff(PbPb_4_HFdown, pp_0), makeDiff(PbPb_0_HFdown, pp_0), makeDiff(PbPb_1_HFdown, pp_0), makeDiff(PbPb_2_HFdown, pp_0), makeDiff(PbPb_3_HFdown, pp_0)};
    sub_1DpT5 = {makeDiff(PbPb_4, pp_0_1D_pT), makeDiff(PbPb_0, pp_0_1D_pT), makeDiff(PbPb_1, pp_0_1D_pT), makeDiff(PbPb_2, pp_0_1D_pT), makeDiff(PbPb_3, pp_0_1D_pT)};

    SaveAllCentralityTGraphs("./nominal_syst_tgraph.root",
                             type,
                             pbpbEll5, ppEll5, subEll5,

                             pbpb_nominal5, pbpb_tnpU5, pbpb_tnpD5, pbpb_acoup5, pbpb_acodown5,
                             pbpb_nobk5, pbpb_massrange5, pbpb_HFup5, pbpb_HFdown5, pbpb_1DpT5,

                             pp_nominal5, pp_tnpU5, pp_tnpD5, pp_acoup5, pp_acodown5,
                             pp_nobk5, pp_massrange5, pp_HFup5, pp_HFdown5, pp_1DpT5,

                             sub_nominal5, sub_tnpU5, sub_tnpD5, sub_acoup5, sub_acodown5,
                             sub_nobk5, sub_massrange5, sub_HFup5, sub_HFdown5, sub_1DpT5);
}