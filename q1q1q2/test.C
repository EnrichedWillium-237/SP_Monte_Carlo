# include "TCanvas.h"
# include "TF1.h"
# include "TH1.h"
# include "TMath.h"
# include "TRandom3.h"

Double_t bounds(int ord, double ang) {
    while (ang >  TMath::Pi()/ord) ang-=TMath::TwoPi()/ord;
    while (ang < -TMath::Pi()/ord) ang+=TMath::TwoPi()/ord;
    return ang;
}

void test() {
    int nevents = 1e4;
    int mult = 1000;
    double v1in = 0.020;
    double v2in = 0.060;

    TH1D * hphi = new TH1D("phi", "", 100, -3.5, 3.5);
    TH1D * hPsi1 = new TH1D("Psi1", "", 100, -3.5, 3.5);
    TH1D * hPsi2 = new TH1D("Psi2", "", 100, -2.2, 2.2);
    TH1D * hq1q1q2 = new TH1D("q1q1q2", "", 100, -0.1, 0.2);
    TF1 * phidist = new TF1("phidist","1+2*[0]*cos(x)+2*[1]*cos(2*x)",-TMath::Pi(),TMath::Pi());
    phidist->SetParameters(v1in,v2in);
    TRandom3 * ran = new TRandom3(0);
    Double_t c1, c2, s1, s2, psi1, psi2;
    for (int k = 0; k<10; k++) {
        double q1q1q2 = 0;
        for (int i = 0; i<nevents; i++) {
            c1 = 0;
            c2 = 0;
            s1 = 0;
            s2 = 0;
            psi1 = 0;
            psi2 = 0;
            // for (int j = 0; j<mult; j++) {
            //     double phi = phidist->GetRandom();
            //     c1 += TMath::Cos(phi);
            //     s1 += TMath::Sin(phi);
            //     c2 += TMath::Cos(2*phi);
            //     s2 += TMath::Sin(2*phi);
            //     hphi->Fill(phi);
            // }
            // psi1 = TMath::ATan2(s1,c1);
            // psi2 = TMath::ATan2(s2,c2)/2.;
            // psi1 = ran->Uniform(-TMath::Pi(), TMath::Pi());
            psi1 = phidist->GetRandom(-TMath::Pi(), TMath::Pi());
            psi2 = phidist->GetRandom(-TMath::Pi()/2., TMath::Pi()/2.);
            hPsi1->Fill(psi1);
            hPsi2->Fill(psi2);
            q1q1q2 += TMath::Cos(2*(psi1 - psi2));
        }
        q1q1q2/=nevents;
        hq1q1q2->Fill(q1q1q2);
        cout<<q1q1q2<<endl;
    }
    cout<<"\n<cos(2(Psi1 - Psi2))>: "<<hq1q1q2->GetMean()<<" +/- "<<hq1q1q2->GetMeanError()<<endl;

    TCanvas * c0 = new TCanvas("c0", "c0", 600, 550);
    c0->cd();
    hPsi1->GetYaxis()->SetRangeUser(0, 20000);
    hPsi1->Draw();
    hPsi2->SetLineColor(kRed);
    hPsi2->Draw("same");

    TCanvas * c11 = new TCanvas("c11", "c11", 600, 500);
    c11->cd();
    hphi->Draw();

}
