TH1D * xprof_HFm1_HFp1;
TH1D * xprof_trackmid1_HFm1;
TH1D * xprof_trackmid1_HFp1;
TH1D * xprof_trackm1_HFm1;
TH1D * xprof_trackm1_trackp1;
TH1D * xprof_trackm1_HFp1;
TH1D * xprof_trackp1_HFm1;
TH1D * xprof_trackp1_trackm1;
TH1D * xprof_trackp1_HFp1;

TH1D * xprof_HFm1_trackm2;
TH1D * xprof_HFm1_trackp2;
TH1D * xprof_HFp1_trackm2;
TH1D * xprof_HFp1_trackp2;
TH1D * xprof_HFm1_HFp2;
TH1D * xprof_HFp1_HFm2;
TH1D * xprof_trackmid1_HFm2;
TH1D * xprof_trackmid1_HFp2;
TH1D * xprof_trackm1_HFm2;
TH1D * xprof_trackm1_trackp2;
TH1D * xprof_trackm1_HFp2;
TH1D * xprof_trackp1_HFm2;
TH1D * xprof_trackp1_trackm2;
TH1D * xprof_trackp1_HFp2;

TH1D * xprof_HFm2_HFp2;
TH1D * xprof_trackmid2_HFm2;
TH1D * xprof_trackmid2_HFp2;
TH1D * xprof_trackm2_HFm2;
TH1D * xprof_trackm2_trackp2;
TH1D * xprof_trackm2_HFp2;
TH1D * xprof_trackp2_HFm2;
TH1D * xprof_trackp2_HFp2;

void EPXprofile()
{
    xprof_HFm1_HFp1 = (TH1D *) corr_HFm1_HFp1->ProfileX()->Clone("xprof_HFm1_HFp1");
    xprof_trackmid1_HFm1 = (TH1D *) corr_trackmid1_HFm1->ProfileX()->Clone("xprof_trackmid1_HFm1");
    xprof_trackmid1_HFp1 = (TH1D *) corr_trackmid1_HFp1->ProfileX()->Clone("xprof_trackmid1_HFp1");
    xprof_trackm1_HFm1 = (TH1D *) corr_trackm1_HFm1->ProfileX()->Clone("xprof_trackm1_HFm1");
    xprof_trackm1_trackp1 = (TH1D *) corr_trackm1_trackp1->ProfileX()->Clone("xprof_trackm1_trackp1");
    xprof_trackm1_HFp1 = (TH1D *) corr_trackm1_HFp1->ProfileX()->Clone("xprof_trackm1_HFp1");
    xprof_trackp1_HFm1 = (TH1D *) corr_trackp1_HFm1->ProfileX()->Clone("xprof_trackp1_HFm1");
    xprof_trackp1_trackm1 = (TH1D *) corr_trackp1_trackm1->ProfileX()->Clone("xprof_trackp1_trackm1");
    xprof_trackp1_HFp1 = (TH1D *) corr_trackp1_HFp1->ProfileX()->Clone("xprof_trackp1_HFp1");

    xprof_HFm1_trackm2 = (TH1D *) corr_HFm1_trackm2->ProfileX()->Clone("xprof_HFm1_trackm2");
    xprof_HFm1_trackp2 = (TH1D *) corr_HFm1_trackp2->ProfileX()->Clone("xprof_HFm1_trackp2");
    xprof_HFp1_trackm2 = (TH1D *) corr_HFp1_trackm2->ProfileX()->Clone("xprof_HFp1_trackm2");
    xprof_HFp1_trackp2 = (TH1D *) corr_HFp1_trackp2->ProfileX()->Clone("xprof_HFp1_trackp2");
    xprof_HFm1_HFp2 = (TH1D *) corr_HFm1_HFp2->ProfileX()->Clone("xprof_HFm1_HFp2");
    xprof_HFp1_HFm2 = (TH1D *) corr_HFp1_HFm2->ProfileX()->Clone("xprof_HFp1_HFm2");
    xprof_trackmid1_HFm2 = (TH1D *) corr_trackmid1_HFm2->ProfileX()->Clone("xprof_trackmid1_HFm2");
    xprof_trackmid1_HFp2 = (TH1D *) corr_trackmid1_HFp2->ProfileX()->Clone("xprof_trackmid1_HFp2");
    xprof_trackm1_HFm2 = (TH1D *) corr_trackm1_HFm2->ProfileX()->Clone("xprof_trackm1_HFm2");
    xprof_trackm1_trackp2 = (TH1D *) corr_trackm1_trackp2->ProfileX()->Clone("xprof_trackm1_trackp2");
    xprof_trackm1_HFp2 = (TH1D *) corr_trackm1_HFp2->ProfileX()->Clone("xprof_trackm1_HFp2");
    xprof_trackp1_HFm2 = (TH1D *) corr_trackp1_HFm2->ProfileX()->Clone("xprof_trackp1_HFm2");
    xprof_trackp1_trackm2 = (TH1D *) corr_trackp1_trackm2->ProfileX()->Clone("xprof_trackp1_trackm2");
    xprof_trackp1_HFp2 = (TH1D *) corr_trackp1_HFp2->ProfileX()->Clone("xprof_trackp1_HFp2");

    xprof_HFm2_HFp2 = (TH1D *) corr_HFm2_HFp2->ProfileX()->Clone("xprof_HFm2_HFp2");
    xprof_trackmid2_HFm2 = (TH1D *) corr_trackmid2_HFm2->ProfileX()->Clone("xprof_trackmid2_HFm2");
    xprof_trackmid2_HFp2 = (TH1D *) corr_trackmid2_HFp2->ProfileX()->Clone("xprof_trackmid2_HFp2");
    xprof_trackm2_HFm2 = (TH1D *) corr_trackm2_HFm2->ProfileX()->Clone("xprof_trackm2_HFm2");
    xprof_trackm2_trackp2 = (TH1D *) corr_trackm2_trackp2->ProfileX()->Clone("xprof_trackm2_trackp2");
    xprof_trackm2_HFp2 = (TH1D *) corr_trackm2_HFp2->ProfileX()->Clone("xprof_trackm2_HFp2");
    xprof_trackp2_HFm2 = (TH1D *) corr_trackp2_HFm2->ProfileX()->Clone("xprof_trackp2_HFm2");
    xprof_trackp2_HFp2 = (TH1D *) corr_trackp2_HFp2->ProfileX()->Clone("xprof_trackp2_HFp2");
}

void EPXprofPrint( TString tag )
{
    if (!fopen(Form("plots/%s/Xprofiles",tag.Data()),"r")) system(Form("mkdir plots/%s/Xprofiles",tag.Data()));
    if (!fopen(Form("plots/%s/Xprofiles/png",tag.Data()),"r")) system(Form("mkdir plots/%s/Xprofiles/png",tag.Data()));
    if (!fopen(Form("plots/%s/Xprofiles/pdf",tag.Data()),"r")) system(Form("mkdir plots/%s/Xprofiles/pdf",tag.Data()));

    TCanvas * cxprof_HFm1_HFp1 = new TCanvas("cxprof_HFm1_HFp1","cxprof_HFm1_HFp1",600,530);
    TPad * padxprof_HFm1_HFp1 = (TPad *) cxprof_HFm1_HFp1->cd();
    padxprof_HFm1_HFp1->SetRightMargin(rghtmrg);
    xprof_HFm1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    xprof_HFm1_HFp1->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFm1_HFp1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFm1_HFp1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFm1_HFp1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFm1_HFp1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFm1_HFp1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFm1_HFp1->Draw();
    if (print_plot) cxprof_HFm1_HFp1->Print(Form("plots/%s/Xprofiles/png/HFp1_vs_HFm1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFm1_HFp1->Print(Form("plots/%s/Xprofiles/pdf/HFp1_vs_HFm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFm1_HFp1->Close();


    TCanvas * cxprof_trackmid1_HFm1 = new TCanvas("cxprof_trackmid1_HFm1","cxprof_trackmid1_HFm1",600,530);
    TPad * padxprof_trackmid1_HFm1 = (TPad *) cxprof_trackmid1_HFm1->cd();
    padxprof_trackmid1_HFm1->SetRightMargin(rghtmrg);
    xprof_trackmid1_HFm1->SetYTitle("#Psi_{1}{HF-}");
    xprof_trackmid1_HFm1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackmid1_HFm1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackmid1_HFm1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackmid1_HFm1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackmid1_HFm1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackmid1_HFm1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackmid1_HFm1->Draw();
    if (print_plot) cxprof_trackmid1_HFm1->Print(Form("plots/%s/Xprofiles/png/HFm1_vs_trackmid1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackmid1_HFm1->Print(Form("plots/%s/Xprofiles/pdf/HFm1_vs_trackmid1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackmid1_HFm1->Close();


    TCanvas * cxprof_trackmid1_HFp1 = new TCanvas("cxprof_trackmid1_HFp1","cxprof_trackmid1_HFp1",600,530);
    TPad * padxprof_trackmid1_HFp1 = (TPad *) cxprof_trackmid1_HFp1->cd();
    padxprof_trackmid1_HFp1->SetRightMargin(rghtmrg);
    xprof_trackmid1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    xprof_trackmid1_HFp1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackmid1_HFp1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackmid1_HFp1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackmid1_HFp1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackmid1_HFp1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackmid1_HFp1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackmid1_HFp1->Draw();
    if (print_plot) cxprof_trackmid1_HFp1->Print(Form("plots/%s/Xprofiles/png/HFp1_vs_trackmid1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackmid1_HFp1->Print(Form("plots/%s/Xprofiles/pdf/HFp1_vs_trackmid1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackmid1_HFp1->Close();


    TCanvas * cxprof_trackm1_HFm1 = new TCanvas("cxprof_trackm1_HFm1","cxprof_trackm1_HFm1",600,530);
    TPad * padxprof_trackm1_HFm1 = (TPad *) cxprof_trackm1_HFm1->cd();
    padxprof_trackm1_HFm1->SetRightMargin(rghtmrg);
    xprof_trackm1_HFm1->SetYTitle("#Psi_{1}{HF-}");
    xprof_trackm1_HFm1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm1_HFm1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm1_HFm1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm1_HFm1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm1_HFm1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm1_HFm1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackm1_HFm1->Draw();
    if (print_plot) cxprof_trackm1_HFm1->Print(Form("plots/%s/Xprofiles/png/HFm1_vs_trackm1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm1_HFm1->Print(Form("plots/%s/Xprofiles/pdf/HFm1_vs_trackm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm1_HFm1->Close();


    TCanvas * cxprof_trackm1_trackp1 = new TCanvas("cxprof_trackm1_trackp1","cxprof_trackm1_trackp1",600,530);
    TPad * padxprof_trackm1_trackp1 = (TPad *) cxprof_trackm1_trackp1->cd();
    padxprof_trackm1_trackp1->SetRightMargin(rghtmrg);
    xprof_trackm1_trackp1->SetYTitle("#Psi_{1}{track+}");
    xprof_trackm1_trackp1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm1_trackp1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm1_trackp1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm1_trackp1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm1_trackp1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm1_trackp1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackm1_trackp1->Draw();
    if (print_plot) cxprof_trackm1_trackp1->Print(Form("plots/%s/Xprofiles/png/trackp1_vs_trackm1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm1_trackp1->Print(Form("plots/%s/Xprofiles/pdf/trackp1_vs_trackm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm1_trackp1->Close();


    TCanvas * cxprof_trackm1_HFp1 = new TCanvas("cxprof_trackm1_HFp1","cxprof_trackm1_HFp1",600,530);
    TPad * padxprof_trackm1_HFp1 = (TPad *) cxprof_trackm1_HFp1->cd();
    padxprof_trackm1_HFp1->SetRightMargin(rghtmrg);
    xprof_trackm1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    xprof_trackm1_HFp1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm1_HFp1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm1_HFp1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm1_HFp1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm1_HFp1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm1_HFp1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackm1_HFp1->Draw();
    if (print_plot) cxprof_trackm1_HFp1->Print(Form("plots/%s/Xprofiles/png/HFp1_vs_trackm1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm1_HFp1->Print(Form("plots/%s/Xprofiles/pdf/HFp1_vs_trackm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm1_HFp1->Close();


    TCanvas * cxprof_trackp1_HFm1 = new TCanvas("cxprof_trackp1_HFm1","cxprof_trackp1_HFm1",600,530);
    TPad * padxprof_trackp1_HFm1 = (TPad *) cxprof_trackp1_HFm1->cd();
    padxprof_trackp1_HFm1->SetRightMargin(rghtmrg);
    xprof_trackp1_HFm1->SetYTitle("#Psi_{1}{HF-}");
    xprof_trackp1_HFm1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp1_HFm1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp1_HFm1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp1_HFm1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp1_HFm1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp1_HFm1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackp1_HFm1->Draw();
    if (print_plot) cxprof_trackp1_HFm1->Print(Form("plots/%s/Xprofiles/png/HFm1_vs_trackp1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp1_HFm1->Print(Form("plots/%s/Xprofiles/pdf/HFm1_vs_trackp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp1_HFm1->Close();


    TCanvas * cxprof_trackp1_trackm1 = new TCanvas("cxprof_trackp1_trackm1","cxprof_trackp1_trackm1",600,530);
    TPad * padxprof_trackp1_trackm1 = (TPad *) cxprof_trackp1_trackm1->cd();
    padxprof_trackp1_trackm1->SetRightMargin(rghtmrg);
    xprof_trackp1_trackm1->SetYTitle("#Psi_{1}{track-}");
    xprof_trackp1_trackm1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp1_trackm1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp1_trackm1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp1_trackm1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp1_trackm1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp1_trackm1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackp1_trackm1->Draw();
    if (print_plot) cxprof_trackp1_trackm1->Print(Form("plots/%s/Xprofiles/png/trackm1_vs_trackp1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp1_trackm1->Print(Form("plots/%s/Xprofiles/pdf/trackm1_vs_trackp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp1_trackm1->Close();


    TCanvas * cxprof_trackp1_HFp1 = new TCanvas("cxprof_trackp1_HFp1","cxprof_trackp1_HFp1",600,530);
    TPad * padxprof_trackp1_HFp1 = (TPad *) cxprof_trackp1_HFp1->cd();
    padxprof_trackp1_HFp1->SetRightMargin(rghtmrg);
    xprof_trackp1_HFp1->SetYTitle("#Psi_{1}{HF+}");
    xprof_trackp1_HFp1->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp1_HFp1->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp1_HFp1->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp1_HFp1->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp1_HFp1->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp1_HFp1->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackp1_HFp1->Draw();
    if (print_plot) cxprof_trackp1_HFp1->Print(Form("plots/%s/Xprofiles/png/HFp1_vs_trackp1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp1_HFp1->Print(Form("plots/%s/Xprofiles/pdf/HFp1_vs_trackp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp1_HFp1->Close();


    TCanvas * cxprof_HFm1_trackm2 = new TCanvas("cxprof_HFm1_trackm2","cxprof_HFm1_trackm2",600,530);
    TPad * padxprof_HFm1_trackm2 = (TPad *) cxprof_HFm1_trackm2->cd();
    padxprof_HFm1_trackm2->SetRightMargin(rghtmrg);
    xprof_HFm1_trackm2->SetYTitle("#Psi_{2}{track-}");
    xprof_HFm1_trackm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFm1_trackm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFm1_trackm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFm1_trackm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFm1_trackm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFm1_trackm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFm1_trackm2->Draw();
    if (print_plot) cxprof_HFm1_trackm2->Print(Form("plots/%s/Xprofiles/png/trackm2_vs_HFm1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFm1_trackm2->Print(Form("plots/%s/Xprofiles/pdf/trackm2_vs_HFm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFm1_trackm2->Close();


    TCanvas * cxprof_HFm1_trackp2 = new TCanvas("cxprof_HFm1_trackp2","cxprof_HFm1_trackp2",600,530);
    TPad * padxprof_HFm1_trackp2 = (TPad *) cxprof_HFm1_trackp2->cd();
    padxprof_HFm1_trackp2->SetRightMargin(rghtmrg);
    xprof_HFm1_trackp2->SetYTitle("#Psi_{2}{track+}");
    xprof_HFm1_trackp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFm1_trackp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFm1_trackp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFm1_trackp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFm1_trackp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFm1_trackp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFm1_trackp2->Draw();
    if (print_plot) cxprof_HFm1_trackp2->Print(Form("plots/%s/Xprofiles/png/trackp2_vs_HFm1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFm1_trackp2->Print(Form("plots/%s/Xprofiles/pdf/trackp2_vs_HFm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFm1_trackp2->Close();


    TCanvas * cxprof_HFp1_trackm2 = new TCanvas("cxprof_HFp1_trackm2","cxprof_HFp1_trackm2",600,530);
    TPad * padxprof_HFp1_trackm2 = (TPad *) cxprof_HFp1_trackm2->cd();
    padxprof_HFp1_trackm2->SetRightMargin(rghtmrg);
    xprof_HFp1_trackm2->SetYTitle("#Psi_{2}{track-}");
    xprof_HFp1_trackm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFp1_trackm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFp1_trackm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFp1_trackm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFp1_trackm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFp1_trackm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFp1_trackm2->Draw();
    if (print_plot) cxprof_HFp1_trackm2->Print(Form("plots/%s/Xprofiles/png/trackm2_vs_HFp1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFp1_trackm2->Print(Form("plots/%s/Xprofiles/pdf/trackm2_vs_HFp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFp1_trackm2->Close();


    TCanvas * cxprof_HFp1_trackp2 = new TCanvas("cxprof_HFp1_trackp2","cxprof_HFp1_trackp2",600,530);
    TPad * padxprof_HFp1_trackp2 = (TPad *) cxprof_HFp1_trackp2->cd();
    padxprof_HFp1_trackp2->SetRightMargin(rghtmrg);
    xprof_HFp1_trackp2->SetYTitle("#Psi_{2}{track+}");
    xprof_HFp1_trackp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFp1_trackp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFp1_trackp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFp1_trackp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFp1_trackp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFp1_trackp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFp1_trackp2->Draw();
    if (print_plot) cxprof_HFp1_trackp2->Print(Form("plots/%s/Xprofiles/png/trackp2_vs_HFp1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFp1_trackp2->Print(Form("plots/%s/Xprofiles/pdf/trackp2_vs_HFp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFp1_trackp2->Close();


    TCanvas * cxprof_HFm1_HFp2 = new TCanvas("cxprof_HFm1_HFp2","cxprof_HFm1_HFp2",600,530);
    TPad * padxprof_HFm1_HFp2 = (TPad *) cxprof_HFm1_HFp2->cd();
    padxprof_HFm1_HFp2->SetRightMargin(rghtmrg);
    xprof_HFm1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_HFm1_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFm1_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFm1_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFm1_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFm1_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFm1_HFp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFm1_HFp2->Draw();
    if (print_plot) cxprof_HFm1_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_HFm1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFm1_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_HFm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFm1_HFp2->Close();


    TCanvas * cxprof_HFp1_HFm2 = new TCanvas("cxprof_HFp1_HFm2","cxprof_HFp1_HFm2",600,530);
    TPad * padxprof_HFp1_HFm2 = (TPad *) cxprof_HFp1_HFm2->cd();
    padxprof_HFp1_HFm2->SetRightMargin(rghtmrg);
    xprof_HFp1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_HFp1_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFp1_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFp1_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFp1_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFp1_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFp1_HFm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_HFp1_HFm2->Draw();
    if (print_plot) cxprof_HFp1_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_HFp1.png",tag.Data()),"png");
    if (print_plot) cxprof_HFp1_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_HFp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFp1_HFm2->Close();


    TCanvas * cxprof_trackmid1_HFm2 = new TCanvas("cxprof_trackmid1_HFm2","cxprof_trackmid1_HFm2",600,530);
    TPad * padxprof_trackmid1_HFm2 = (TPad *) cxprof_trackmid1_HFm2->cd();
    padxprof_trackmid1_HFm2->SetRightMargin(rghtmrg);
    xprof_trackmid1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_trackmid1_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackmid1_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackmid1_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackmid1_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackmid1_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackmid1_HFm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackmid1_HFm2->Draw();
    if (print_plot) cxprof_trackmid1_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_trackmid1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackmid1_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_trackmid1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackmid1_HFm2->Close();


    TCanvas * cxprof_trackmid1_HFp2 = new TCanvas("cxprof_trackmid1_HFp2","cxprof_trackmid1_HFp2",600,530);
    TPad * padxprof_trackmid1_HFp2 = (TPad *) cxprof_trackmid1_HFp2->cd();
    padxprof_trackmid1_HFp2->SetRightMargin(rghtmrg);
    xprof_trackmid1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_trackmid1_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackmid1_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackmid1_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackmid1_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackmid1_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackmid1_HFp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackmid1_HFp2->Draw();
    if (print_plot) cxprof_trackmid1_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_trackmid1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackmid1_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_trackmid1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackmid1_HFp2->Close();


    TCanvas * cxprof_trackm1_HFm2 = new TCanvas("cxprof_trackm1_HFm2","cxprof_trackm1_HFm2",600,530);
    TPad * padxprof_trackm1_HFm2 = (TPad *) cxprof_trackm1_HFm2->cd();
    padxprof_trackm1_HFm2->SetRightMargin(rghtmrg);
    xprof_trackm1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_trackm1_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm1_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm1_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm1_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm1_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm1_HFm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackm1_HFm2->Draw();
    if (print_plot) cxprof_trackm1_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_trackm1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm1_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_trackm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm1_HFm2->Close();


    TCanvas * cxprof_trackm1_trackp2 = new TCanvas("cxprof_trackm1_trackp2","cxprof_trackm1_trackp2",600,530);
    TPad * padxprof_trackm1_trackp2 = (TPad *) cxprof_trackm1_trackp2->cd();
    padxprof_trackm1_trackp2->SetRightMargin(rghtmrg);
    xprof_trackm1_trackp2->SetYTitle("#Psi_{2}{track+}");
    xprof_trackm1_trackp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm1_trackp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm1_trackp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm1_trackp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm1_trackp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm1_trackp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackm1_trackp2->Draw();
    if (print_plot) cxprof_trackm1_trackp2->Print(Form("plots/%s/Xprofiles/png/trackp2_vs_trackm1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm1_trackp2->Print(Form("plots/%s/Xprofiles/pdf/trackp2_vs_trackm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm1_trackp2->Close();


    TCanvas * cxprof_trackm1_HFp2 = new TCanvas("cxprof_trackm1_HFp2","cxprof_trackm1_HFp2",600,530);
    TPad * padxprof_trackm1_HFp2 = (TPad *) cxprof_trackm1_HFp2->cd();
    padxprof_trackm1_HFp2->SetRightMargin(rghtmrg);
    xprof_trackm1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_trackm1_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm1_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm1_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm1_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm1_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm1_HFp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackm1_HFp2->Draw();
    if (print_plot) cxprof_trackm1_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_trackm1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm1_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_trackm1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm1_HFp2->Close();


    TCanvas * cxprof_trackp1_HFm2 = new TCanvas("cxprof_trackp1_HFm2","cxprof_trackp1_HFm2",600,530);
    TPad * padxprof_trackp1_HFm2 = (TPad *) cxprof_trackp1_HFm2->cd();
    padxprof_trackp1_HFm2->SetRightMargin(rghtmrg);
    xprof_trackp1_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_trackp1_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp1_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp1_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp1_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp1_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp1_HFm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackp1_HFm2->Draw();
    if (print_plot) cxprof_trackp1_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_trackp1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp1_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_trackp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp1_HFm2->Close();


    TCanvas * cxprof_trackp1_trackm2 = new TCanvas("cxprof_trackp1_trackm2","cxprof_trackp1_trackm2",600,530);
    TPad * padxprof_trackp1_trackm2 = (TPad *) cxprof_trackp1_trackm2->cd();
    padxprof_trackp1_trackm2->SetRightMargin(rghtmrg);
    xprof_trackp1_trackm2->SetYTitle("#Psi_{2}{track-}");
    xprof_trackp1_trackm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp1_trackm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp1_trackm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp1_trackm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp1_trackm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp1_trackm2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackp1_trackm2->Draw();
    if (print_plot) cxprof_trackp1_trackm2->Print(Form("plots/%s/Xprofiles/png/trackm2_vs_trackp1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp1_trackm2->Print(Form("plots/%s/Xprofiles/pdf/trackm2_vs_trackp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp1_trackm2->Close();


    TCanvas * cxprof_trackp1_HFp2 = new TCanvas("cxprof_trackp1_HFp2","cxprof_trackp1_HFp2",600,530);
    TPad * padxprof_trackp1_HFp2 = (TPad *) cxprof_trackp1_HFp2->cd();
    padxprof_trackp1_HFp2->SetRightMargin(rghtmrg);
    xprof_trackp1_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_trackp1_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp1_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp1_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp1_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp1_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp1_HFp2->GetYaxis()->SetRangeUser(xprofymin,xprofymax);
    xprof_trackp1_HFp2->Draw();
    if (print_plot) cxprof_trackp1_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_trackp1.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp1_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_trackp1.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp1_HFp2->Close();


    TCanvas * cxprof_HFm2_HFp2 = new TCanvas("cxprof_HFm2_HFp2","cxprof_HFm2_HFp2",600,530);
    TPad * padxprof_HFm2_HFp2 = (TPad *) cxprof_HFm2_HFp2->cd();
    padxprof_HFm2_HFp2->SetRightMargin(rghtmrg);
    xprof_HFm2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_HFm2_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_HFm2_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_HFm2_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_HFm2_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_HFm2_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_HFm2_HFp2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_HFm2_HFp2->Draw();
    if (print_plot) cxprof_HFm2_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_HFm2.png",tag.Data()),"png");
    if (print_plot) cxprof_HFm2_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_HFm2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_HFm2_HFp2->Close();


    TCanvas * cxprof_trackmid2_HFm2 = new TCanvas("cxprof_trackmid2_HFm2","cxprof_trackmid2_HFm2",600,530);
    TPad * padxprof_trackmid2_HFm2 = (TPad *) cxprof_trackmid2_HFm2->cd();
    padxprof_trackmid2_HFm2->SetRightMargin(rghtmrg);
    xprof_trackmid2_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_trackmid2_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackmid2_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackmid2_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackmid2_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackmid2_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackmid2_HFm2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackmid2_HFm2->Draw();
    if (print_plot) cxprof_trackmid2_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_trackmid2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackmid2_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_trackmid2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackmid2_HFm2->Close();


    TCanvas * cxprof_trackmid2_HFp2 = new TCanvas("cxprof_trackmid2_HFp2","cxprof_trackmid2_HFp2",600,530);
    TPad * padxprof_trackmid2_HFp2 = (TPad *) cxprof_trackmid2_HFp2->cd();
    padxprof_trackmid2_HFp2->SetRightMargin(rghtmrg);
    xprof_trackmid2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_trackmid2_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackmid2_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackmid2_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackmid2_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackmid2_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackmid2_HFp2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackmid2_HFp2->Draw();
    if (print_plot) cxprof_trackmid2_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_trackmid2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackmid2_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_trackmid2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackmid2_HFp2->Close();


    TCanvas * cxprof_trackm2_HFm2 = new TCanvas("cxprof_trackm2_HFm2","cxprof_trackm2_HFm2",600,530);
    TPad * padxprof_trackm2_HFm2 = (TPad *) cxprof_trackm2_HFm2->cd();
    padxprof_trackm2_HFm2->SetRightMargin(rghtmrg);
    xprof_trackm2_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_trackm2_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm2_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm2_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm2_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm2_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm2_HFm2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackm2_HFm2->Draw();
    if (print_plot) cxprof_trackm2_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_trackm2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm2_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_trackm2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm2_HFm2->Close();


    TCanvas * cxprof_trackm2_trackp2 = new TCanvas("cxprof_trackm2_trackp2","cxprof_trackm2_trackp2",600,530);
    TPad * padxprof_trackm2_trackp2 = (TPad *) cxprof_trackm2_trackp2->cd();
    padxprof_trackm2_trackp2->SetRightMargin(rghtmrg);
    xprof_trackm2_trackp2->SetYTitle("#Psi_{2}{track+}");
    xprof_trackm2_trackp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm2_trackp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm2_trackp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm2_trackp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm2_trackp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm2_trackp2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackm2_trackp2->Draw();
    if (print_plot) cxprof_trackm2_trackp2->Print(Form("plots/%s/Xprofiles/png/trackp2_vs_trackm2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm2_trackp2->Print(Form("plots/%s/Xprofiles/pdf/trackp2_vs_trackm2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm2_trackp2->Close();


    TCanvas * cxprof_trackm2_HFp2 = new TCanvas("cxprof_trackm2_HFp2","cxprof_trackm2_HFp2",600,530);
    TPad * padxprof_trackm2_HFp2 = (TPad *) cxprof_trackm2_HFp2->cd();
    padxprof_trackm2_HFp2->SetRightMargin(rghtmrg);
    xprof_trackm2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_trackm2_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackm2_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackm2_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackm2_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackm2_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackm2_HFp2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackm2_HFp2->Draw();
    if (print_plot) cxprof_trackm2_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_trackm2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackm2_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_trackm2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackm2_HFp2->Close();


    TCanvas * cxprof_trackp2_HFm2 = new TCanvas("cxprof_trackp2_HFm2","cxprof_trackp2_HFm2",600,530);
    TPad * padxprof_trackp2_HFm2 = (TPad *) cxprof_trackp2_HFm2->cd();
    padxprof_trackp2_HFm2->SetRightMargin(rghtmrg);
    xprof_trackp2_HFm2->SetYTitle("#Psi_{2}{HF-}");
    xprof_trackp2_HFm2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp2_HFm2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp2_HFm2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp2_HFm2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp2_HFm2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp2_HFm2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackp2_HFm2->Draw();
    if (print_plot) cxprof_trackp2_HFm2->Print(Form("plots/%s/Xprofiles/png/HFm2_vs_trackp2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp2_HFm2->Print(Form("plots/%s/Xprofiles/pdf/HFm2_vs_trackp2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp2_HFm2->Close();


    TCanvas * cxprof_trackp2_HFp2 = new TCanvas("cxprof_trackp2_HFp2","cxprof_trackp2_HFp2",600,530);
    TPad * padxprof_trackp2_HFp2 = (TPad *) cxprof_trackp2_HFp2->cd();
    padxprof_trackp2_HFp2->SetRightMargin(rghtmrg);
    xprof_trackp2_HFp2->SetYTitle("#Psi_{2}{HF+}");
    xprof_trackp2_HFp2->GetYaxis()->CenterTitle(kTRUE);
    xprof_trackp2_HFp2->GetYaxis()->SetTitleSize(ytlsize1x1);
    xprof_trackp2_HFp2->GetYaxis()->SetTitleOffset(ytloffset1x1);
    xprof_trackp2_HFp2->GetYaxis()->SetLabelSize(ylbsize1x1);
    xprof_trackp2_HFp2->GetYaxis()->SetLabelOffset(ylboffset1x1);
    xprof_trackp2_HFp2->GetYaxis()->SetRangeUser(xprofymin2ndOrd,xprofymax2ndOrd);
    xprof_trackp2_HFp2->Draw();
    if (print_plot) cxprof_trackp2_HFp2->Print(Form("plots/%s/Xprofiles/png/HFp2_vs_trackp2.png",tag.Data()),"png");
    if (print_plot) cxprof_trackp2_HFp2->Print(Form("plots/%s/Xprofiles/pdf/HFp2_vs_trackp2.pdf",tag.Data()),"pdf");
    if (close_plot) cxprof_trackp2_HFp2->Close();

}
