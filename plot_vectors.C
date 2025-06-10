#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TAxis3D.h>
#include <TLatex.h>

void plot_vectors()
{
    // Define vector properties
    double pT1 = 10.0, theta1 = TMath::DegToRad() * 0.0, phi1 = 0.0;
    double pT2 = 10.0, theta2 = TMath::DegToRad() * 180.0, phi2 = TMath::Pi() / 2.0;

    // Convert to Cartesian coordinates
    double px1 = pT1 * TMath::Cos(phi1);
    double py1 = pT1 * TMath::Sin(phi1);
    double pz1 = pT1 * TMath::Cos(theta1);

    double px2 = pT2 * TMath::Cos(phi2);
    double py2 = pT2 * TMath::Sin(phi2);
    double pz2 = pT2 * TMath::Cos(theta2);

    // Create a canvas
    TCanvas *c = new TCanvas("c", "3D Vector Visualization", 800, 600);

    // Create 3D graph for visualization
    TPolyLine3D *vector1 = new TPolyLine3D(2);
    vector1->SetPoint(0, 0, 0, 0); // Origin
    vector1->SetPoint(1, px1, py1, pz1);
    vector1->SetLineColor(kRed);
    vector1->SetLineWidth(3);

    TPolyLine3D *vector2 = new TPolyLine3D(2);
    vector2->SetPoint(0, 0, 0, 0); // Origin
    vector2->SetPoint(1, px2, py2, pz2);
    vector2->SetLineColor(kBlue);
    vector2->SetLineWidth(3);

    // Graph 3D axes for context
    TGraph2D *axis = new TGraph2D();
    axis->SetPoint(0, 0, 0, 0);
    axis->SetPoint(1, 12, 0, 0); // X-axis
    axis->SetPoint(2, 0, 12, 0); // Y-axis
    axis->SetPoint(3, 0, 0, 12); // Z-axis

    axis->SetMarkerStyle(20);
    axis->SetTitle("Back-to-Back in #theta, Different in #phi");
    axis->GetXaxis()->SetTitle("p_{x}");
    axis->GetYaxis()->SetTitle("p_{y}");
    axis->GetZaxis()->SetTitle("p_{z}");

    // Draw
    axis->Draw("P");
    vector1->Draw();
    vector2->Draw();

    // Add labels
    /*TLatex text;
    text.SetTextSize(0.03);
    text.DrawLatex(px1, py1, pz1, "#vec{p}_{1} (#theta=60^{#circ}, #phi=0)");
    text.DrawLatex(px2, py2, pz2, "#vec{p}_{2} (#theta=120^{#circ}, #phi=90^{#circ})");*/

    c->Update();
}