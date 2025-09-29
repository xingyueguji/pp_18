void Add_theory_prediction(const char *txtfile = "dM-dsigma.txt")
{
    // --- open old file ---
    TFile *f1 = new TFile("./send_to_frank.root", "READ");
    TCanvas *c_dM_sub = (TCanvas *)f1->Get("dM_PbPb_sub_pp");
    TCanvas *c_dW_sub = (TCanvas *)f1->Get("dW_PbPb_sub_pp");

    // --- extract stat graph and syst boxes (ΔM canvas) ---
    TGraphErrors *g_dM_stat = nullptr;
    std::vector<TBox *> dM_systBoxes;
    {
        TIter next(c_dM_sub->GetListOfPrimitives());
        TObject *obj;
        while ((obj = next()))
        {
            if (obj->InheritsFrom(TGraphErrors::Class()))
            {
                g_dM_stat = (TGraphErrors *)obj->Clone("g_dM_stat_clone");
            }
            else if (obj->InheritsFrom(TBox::Class()))
            {
                TBox *b = (TBox *)obj;
                // keep only systematic boxes (exclude frame/legend)
                if (b->GetFillColor() == 1179 &&
                    !(b->GetX1() <= 0.5 && b->GetX2() >= 5.5))
                {
                    dM_systBoxes.push_back((TBox *)b->Clone());
                }
            }
        }
    }

    // --- extract stat graph and syst boxes (ΔW canvas) ---
    TGraphErrors *g_dW_stat = nullptr;
    std::vector<TBox *> dW_systBoxes;
    {
        TIter next(c_dW_sub->GetListOfPrimitives());
        TObject *obj;
        while ((obj = next()))
        {
            if (obj->InheritsFrom(TGraphErrors::Class()))
            {
                g_dW_stat = (TGraphErrors *)obj->Clone("g_dW_stat_clone");
            }
            else if (obj->InheritsFrom(TBox::Class()))
            {
                TBox *b = (TBox *)obj;
                // keep only systematic boxes (exclude frame/legend)
                if (b->GetFillColor() == 1179 &&
                    !(b->GetX1() <= 0.5 && b->GetX2() >= 5.5))
                {
                    dW_systBoxes.push_back((TBox *)b->Clone());
                }
            }
        }
    }

    // --- read theory file ---
    std::ifstream fin(txtfile);
    if (!fin.is_open())
    {
        std::cerr << "Cannot open txt file: " << txtfile << std::endl;
        return;
    }

    std::vector<std::vector<double>> dMs, dWs;
    std::vector<double> xs, tauBs, as;

    std::string line;
    bool inBlock = false;
    std::vector<double> curM, curW;
    double curX = 0, curTau = 0, curA = 0;

    while (std::getline(fin, line))
    {
        if (line.empty())
            continue;
        if (line.find("x=") != std::string::npos)
        {
            if (inBlock)
            {
                xs.push_back(curX);
                tauBs.push_back(curTau);
                as.push_back(curA);
                dMs.push_back(curM);
                dWs.push_back(curW);
                curM.clear();
                curW.clear();
            }
            inBlock = true;
            std::istringstream ss(line);
            std::string tmp;
            ss >> tmp;
            curX = std::stod(tmp.substr(tmp.find("=") + 1));
            ss >> tmp;
            curTau = std::stod(tmp.substr(tmp.find("=") + 1));
            ss >> tmp;
            curA = std::stod(tmp.substr(tmp.find("=") + 1));
        }
        else if (line.find("Bin") != std::string::npos)
        {
            size_t posM = line.find("ΔM=");
            size_t posS = line.find("Δσ=");
            if (posM == std::string::npos || posS == std::string::npos)
                continue;
            double dM = std::stod(line.substr(posM + 4, line.find("GeV", posM) - (posM + 4)));
            double dW = std::stod(line.substr(posS + 5, line.find("GeV", posS) - (posS + 5)));
            curM.push_back(dM);
            curW.push_back(dW);
        }
    }
    if (inBlock)
    {
        xs.push_back(curX);
        tauBs.push_back(curTau);
        as.push_back(curA);
        dMs.push_back(curM);
        dWs.push_back(curW);
    }

    // --- make new canvases ---
    TCanvas *c_dM_new = new TCanvas("c_dM_new", "ΔM with theory", 800, 800);
    TCanvas *c_dW_new = new TCanvas("c_dW_new", "Δσ with theory", 800, 800);

    // Frame histograms (x = 5 bins)
    int n_cent = 5;
    TH1F *frameM = new TH1F("frameM", ";Centrality;#DeltaM (PbPb - pp)", n_cent, 0.5, n_cent + 0.5);
    frameM->SetMinimum(-0.6);
    frameM->SetMaximum(0.2);
    frameM->GetXaxis()->SetBinLabel(1, "0-10%");
    frameM->GetXaxis()->SetBinLabel(2, "10-20%");
    frameM->GetXaxis()->SetBinLabel(3, "20-30%");
    frameM->GetXaxis()->SetBinLabel(4, "30-100%");
    frameM->GetXaxis()->SetBinLabel(5, "0-100%");

    TH1F *frameW = new TH1F("frameW", ";Centrality;#Delta#sigma (PbPb - pp)", n_cent, 0.5, n_cent + 0.5);
    frameW->SetFillStyle(0); // don’t paint frame background
    frameW->SetMinimum(-0.2);
    frameW->SetMaximum(1.0);
    frameW->GetXaxis()->SetBinLabel(1, "0-10%");
    frameW->GetXaxis()->SetBinLabel(2, "10-20%");
    frameW->GetXaxis()->SetBinLabel(3, "20-30%");
    frameW->GetXaxis()->SetBinLabel(4, "30-100%");
    frameW->GetXaxis()->SetBinLabel(5, "0-100%");

    // --- draw ΔM data first ---
    c_dM_new->cd();
    frameM->Draw();

    // --- draw Δσ data first ---
    c_dW_new->cd();
    frameW->Draw();

    // --- add theory predictions ---
    int color = 2;
    for (size_t i = 0; i < dMs.size(); i++)
    {
        int n = dMs[i].size();
        std::vector<double> xvals(n);
        for (int j = 0; j < n; j++)
            xvals[j] = j + 1;

        TGraph *grM = new TGraph(n, xvals.data(), dMs[i].data());
        grM->SetLineColor(color);
        grM->SetMarkerColor(color);
        grM->SetMarkerStyle(24);
        grM->SetLineWidth(2);
        c_dM_new->cd();
        grM->Draw("PL SAME");

        TGraph *grW = new TGraph(n, xvals.data(), dWs[i].data());
        grW->SetLineColor(color);
        grW->SetMarkerColor(color);
        grW->SetMarkerStyle(25);
        grW->SetLineWidth(2);
        c_dW_new->cd();
        grW->Draw("PL SAME");

        color++;
        if (color == 5)
            color = 6;
    }
    c_dM_new->cd();

    for (auto *b : dM_systBoxes)
    {
        b->SetFillColorAlpha(kGray + 1, 0.35); // gray, 35% opaque
        b->SetLineColor(kGray + 2);            // darker outline
        b->SetLineWidth(2);
        b->Draw("same");
    }
    if (g_dM_stat)
        g_dM_stat->Draw("P SAME");

    c_dW_new->cd();

    for (auto *b : dW_systBoxes)
    {
        b->SetFillColorAlpha(kBlue - 9, 0.25); // light blue, 25% opaque
        b->SetLineColor(kBlue - 7);
        b->SetLineStyle(2);
        b->SetLineWidth(2);
        b->Draw("same");
    }
    if (g_dW_stat)
        g_dW_stat->Draw("P SAME");

    // --- save results ---
    c_dM_new->SaveAs("./final_plot_with_theory/dM_with_theory.pdf");
    c_dW_new->SaveAs("./final_plot_with_theory/dW_with_theory.pdf");
}