# include "EPCorr.h"

void RunCase( Int_t indx, Int_t nevents, Int_t evtmult, Bool_t v1odd, Double_t setv1, Double_t setv2, Double_t setv3, Bool_t eta_weights, Bool_t pt_weights, Bool_t conserve_pT, Bool_t addHoles, Bool_t flatten, Bool_t recenter, Int_t iseed, TString mtag, Int_t Ntries )
{

    int cindx = indx;
    TRandom3 * ran = new TRandom3(iseed+cindx);
    MCEvent * event[nv1Ebins];

    for (int vbin = 0; vbin<nv1Ebins; vbin++) {
        hinit_philab[vbin] = new TH1D(Form("init_philab_v%i_%i",vbin,cindx), "", 150, -4, 4);
        hinit_eta[vbin] = new TH1D(Form("init_eta_v%i_%i",vbin,cindx), "", 150, -6, 6);
        hinit_pt[vbin] = new TH1D(Form("init_pt_v%i_%i",vbin,cindx), "", 150, 0, 8);
        hinit_phiPsiRP[vbin] = new TH1D(Form("init_phiPsiRP_v%i_%i",vbin,cindx), "", 150, -4, 4);
        event[vbin] = new MCEvent(0.0, setv2, setv3, 0.0, 0.0, 0.0);
    }
    hphilab     = new TH1D(Form("philab_%i",   cindx), "", 150, -4,  4);
    hetain      = new TH1D(Form("etain_%i",    cindx), "", 150, -6,  6);
    hptin       = new TH1D(Form("ptin_%i",     cindx), "", 150,  0,  8);
    hphiPsiRP   = new TH1D(Form("phiPsiRP_%i", cindx), "", 150, -4,  4);
    hptAve2     = new TH1D(Form("ptAve2_%i",   cindx), "", 150,  0, 64);
    hweights    = new TH1D(Form("weights_%i",  cindx), "", 150, -1,  1);
    hwpt        = new TH1D(Form("wpt_%i",      cindx), "", 150, -6,  6);
    hpt2D       = new TH2D(Form("pt2D_%i",     cindx), "", 101, -300, 300, 101, -300, 300);
    hpt2D->SetOption("colz");

    for (int nep = 0; nep<numEP; nep++) {
        hsub_philab[nep]    = new TH1D(Form("philab_ep%i_%i",  nep,cindx), "", 150, -4,  4);
        hsub_etain[nep]     = new TH1D(Form("etain_ep%i_%i",   nep,cindx), "", 150, -6,  6);
        hsub_ptin[nep]      = new TH1D(Form("ptin_ep%i_%i",    nep,cindx), "", 150,  0,  8);
        hsub_phiPsiRP[nep]  = new TH1D(Form("phiPsiRP_ep%i_%i",nep,cindx), "", 150, -4,  4);
        hsub_phiPsi1[nep]   = new TH1D(Form("phiPsi1_ep%i_%i", nep,cindx), "", 150, -3.5,  3.5);
        hsub_phiPsi2[nep]   = new TH1D(Form("phiPsi2_ep%i_%i", nep,cindx), "", 150, -2,  2);
        hsub_phiPsi3[nep]   = new TH1D(Form("phiPsi3_ep%i_%i", nep,cindx), "", 150, -1.5,  1.5);
        hsub_Psi1[nep]      = new TH1D(Form("Psi1_ep%i_%i",    nep,cindx), "", 150, -3.5,  3.5);
        hsub_Psi2[nep]      = new TH1D(Form("Psi2_ep%i_%i",    nep,cindx), "", 150, -2,  2);
        hsub_Psi3[nep]      = new TH1D(Form("Psi3_ep%i_%i",    nep,cindx), "", 150, -1.5,  1.5);
        hsub_Psi1lab[nep]   = new TH1D(Form("Psi1lab_ep%i_%i", nep,cindx), "", 150, -3.5,  3.5);
        hsub_Psi2lab[nep]   = new TH1D(Form("Psi2lab_ep%i_%i", nep,cindx), "", 150, -2,  2);
        hsub_Psi3lab[nep]   = new TH1D(Form("Psi3lab_ep%i_%i", nep,cindx), "", 150, -1.5,  1.5);
        hsub_ptAve[nep]     = new TH1D(Form("ptAve_ep%i_%i",   nep,cindx), "", 150,  0,  8);
        hsub_ptAve2[nep]    = new TH1D(Form("ptAve2_ep%i_%i",  nep,cindx), "", 150,  0, 64);
        hsub_weights[nep]   = new TH1D(Form("weights_ep%i_%i", nep,cindx), "", 150, -1,  1);
        hsub_wpt[nep]       = new TH1D(Form("wpt_ep%i_%i",     nep,cindx), "", 150, -6,  6);
        hsub_pt2D[nep]      = new TH2D(Form("pt2D_ep%i_%i",    nep,cindx), "", 101, -300, 300, 101, -300, 300);
        hsub_pt2D[nep]->SetOption("colz");

        flat1[nep] = new HiEvtPlaneFlatten();
        flat1[nep]->Init(9,1,1,1);
        flat2[nep] = new HiEvtPlaneFlatten();
        flat2[nep]->Init(9,1,1,2);
        flat3[nep] = new HiEvtPlaneFlatten();
        flat3[nep]->Init(9,1,1,3);
    }

    Double_t X1[numEP];
    Double_t Y1[numEP];
    Double_t X2[numEP];
    Double_t Y2[numEP];
    Double_t X3[numEP];
    Double_t Y3[numEP];

    Double_t c1[numEP];
    Double_t s1[numEP];
    Double_t c2[numEP];
    Double_t s2[numEP];
    Double_t c3[numEP];
    Double_t s3[numEP];

    Double_t rescor1_2;
    Double_t rescor1_3[2];
    Double_t rescor2HF_2;
    Double_t rescor2HF_3[2];
    Double_t rescor2Trk_3[2];
    Double_t rescor3HF_2;
    Double_t rescor3HF_3[2];
    Double_t rescor3Trk_3[2];
    Double_t rescorMix[2];
    Double_t rescorMix_3[2];
    Double_t epv1obs[2];
    Double_t epv2obs[2];
    Double_t epv3obs[2];
    Double_t epv1_2SE[2];
    Double_t epv1_3SE[2];
    Double_t epv2_2SE[2];
    Double_t epv2_3SE[2];
    Double_t epv3_2SE[2];
    Double_t epv3_3SE[2];
    Double_t epv1obs_mix[2];
    Double_t epv1_mix[2];
    Double_t epv1obs_mix_3[2];
    Double_t epv1_mix_3[2];

    Double_t spdenom1_2SE;
    Double_t spdenom1_3SE[2];
    Double_t spdenom2HF_2SE;
    Double_t spdenom2HF_3SE[2];
    Double_t spdenom2Trk_3SE[2];
    Double_t spdenom3HF_2SE;
    Double_t spdenom3HF_3SE[2];
    Double_t spdenom3Trk_3SE[2];
    Double_t spdenom_mix[2];
    Double_t spdenom_mix_3[2];
    Double_t spv1num[2];
    Double_t spv2num[2];
    Double_t spv3num[2];
    Double_t spv1_2SE[2];
    Double_t spv1_3SE[2];
    Double_t spv2_2SE[2];
    Double_t spv2_3SE[2];
    Double_t spv3_2SE[2];
    Double_t spv3_3SE[2];
    Double_t spv1num_mix[2];
    Double_t spv1_mix[2];
    Double_t spv1num_mix_3[2];
    Double_t spv1_mix_3[2];

    Double_t q1x[numEP];
    Double_t q1y[numEP];
    Double_t q2x[numEP];
    Double_t q2y[numEP];
    Double_t q3x[numEP];
    Double_t q3y[numEP];
    Double_t w1[numEP];
    Double_t w2[numEP];
    Double_t w3[numEP];
    Double_t qcnt[numEP];
    Double_t Psi1[numEP];
    Double_t Psi2[numEP];
    Double_t Psi3[numEP];
    Double_t subptave[numEP];
    Double_t subptave2[numEP];
    Int_t subcnt[numEP];

    comp Q1nA[2];
    comp Q1AB[2];
    comp Q1AC[2];
    comp Q1BC[2];
    comp Q2nA[2];
    comp Q2AB[2];
    comp Q2AC[2];
    comp Q2BC[2];
    comp Q3nA[2];
    comp Q3AB[2];
    comp Q3AC[2];
    comp Q3BC[2];
    comp QmixnAC[2];
    comp QmixABC[2];
    comp QmixnAC_3[2];
    comp QmixABC_3[2];

    comp Q1nAnorm[2];
    comp Q1ABnorm[2];
    comp Q1ACnorm[2];
    comp Q1BCnorm[2];
    comp Q2nAnorm[2];
    comp Q2ABnorm[2];
    comp Q2ACnorm[2];
    comp Q2BCnorm[2];
    comp Q3nAnorm[2];
    comp Q3ABnorm[2];
    comp Q3ACnorm[2];
    comp Q3BCnorm[2];
    comp QmixnACnorm[2];
    comp QmixABCnorm[2];
    comp QmixnACnorm_3[2];
    comp QmixABCnorm_3[2];

    double w1nAsum[2];
    double w1ABsum[2];
    double w1ACsum[2];
    double w1BCsum[2];
    double w2nAsum[2];
    double w2ABsum[2];
    double w2ACsum[2];
    double w2BCsum[2];
    double w3nAsum[2];
    double w3ABsum[2];
    double w3ACsum[2];
    double w3BCsum[2];
    double wmixnACsum[2];
    double wmixABCsum[2];
    double wmixnACsum_3[2];
    double wmixABCsum_3[2];

    //-- differential v1(pT)
    Double_t epv1obs_pt[2][nptbins];
    Double_t epv2obs_pt[2][nptbins];
    Double_t epv3obs_pt[2][nptbins];
    Double_t epv1_2SE_pt[2][nptbins];
    Double_t epv1_3SE_pt[2][nptbins];
    Double_t epv2_2SE_pt[2][nptbins];
    Double_t epv2_3SE_pt[2][nptbins];
    Double_t epv3_2SE_pt[2][nptbins];
    Double_t epv3_3SE_pt[2][nptbins];
    Double_t epv1obs_mix_pt[2][nptbins];
    Double_t epv1_mix_pt[2][nptbins];
    Double_t epv1obs_mix_3_pt[2][nptbins];
    Double_t epv1_mix_3_pt[2][nptbins];

    Double_t spv1num_pt[2][nptbins];
    Double_t spv2num_pt[2][nptbins];
    Double_t spv3num_pt[2][nptbins];
    Double_t spv1_2SE_pt[2][nptbins];
    Double_t spv1_3SE_pt[2][nptbins];
    Double_t spv2_2SE_pt[2][nptbins];
    Double_t spv2_3SE_pt[2][nptbins];
    Double_t spv3_2SE_pt[2][nptbins];
    Double_t spv3_3SE_pt[2][nptbins];
    Double_t spv1num_mix_pt[2][nptbins];
    Double_t spv1_mix_pt[2][nptbins];
    Double_t spv1num_mix_3_pt[2][nptbins];
    Double_t spv1_mix_3_pt[2][nptbins];

    Double_t q1x_pt[2][nptbins];
    Double_t q1y_pt[2][nptbins];
    Double_t q2x_pt[2][nptbins];
    Double_t q2y_pt[2][nptbins];
    Double_t q3x_pt[2][nptbins];
    Double_t q3y_pt[2][nptbins];
    Double_t w1_pt[2][nptbins];
    Double_t w2_pt[2][nptbins];
    Double_t w3_pt[2][nptbins];
    Double_t qcnt_pt[2][nptbins];
    Double_t subptave_pt[2][nptbins];
    Double_t subptave2_pt[2][nptbins];
    Int_t subcnt_pt[2][nptbins];

    comp Q1nA_pt[2][nptbins];
    comp Q2nA_pt[2][nptbins];
    comp Q3nA_pt[2][nptbins];
    comp QmixnAC_pt[2][nptbins];
    comp QmixnAC_3_pt[2][nptbins];

    comp Q1nAnorm_pt[2][nptbins];
    comp Q2nAnorm_pt[2][nptbins];
    comp Q3nAnorm_pt[2][nptbins];
    comp QmixnACnorm_pt[2][nptbins];
    comp QmixnACnorm_3_pt[2][nptbins];

    double w1nAsum_pt[2][nptbins];
    double w2nAsum_pt[2][nptbins];
    double w3nAsum_pt[2][nptbins];
    double wmixnACsum_pt[2][nptbins];
    double wmixnACsum_3_pt[2][nptbins];

    int ncount_pt[2][nptbins];

    //-- differential v1(eta)
    Double_t epv1obs_eta[2][netabins];
    Double_t epv1_2SE_eta[2][netabins];
    Double_t epv1_3SE_eta[2][netabins];
    Double_t epv1obs_mix_eta[2][netabins];
    Double_t epv1_mix_eta[2][netabins];
    Double_t epv1obs_mix_3_eta[2][netabins];
    Double_t epv1_mix_3_eta[2][netabins];

    Double_t spv1num_eta[2][netabins];
    Double_t spv1_2SE_eta[2][netabins];
    Double_t spv1_3SE_eta[2][netabins];
    Double_t spv1num_mix_eta[2][netabins];
    Double_t spv1_mix_eta[2][netabins];
    Double_t spv1num_mix_3_eta[2][netabins];
    Double_t spv1_mix_3_eta[2][netabins];

    Double_t q1x_eta[2][netabins];
    Double_t q1y_eta[2][netabins];
    Double_t w1_eta[2][netabins];
    Double_t qcnt_eta[2][netabins];
    Double_t subetaave_eta[2][netabins];
    Double_t subetaave2_eta[2][netabins];
    Int_t subcnt_eta[2][netabins];

    comp Q1nA_eta[2][netabins];
    comp QmixnAC_eta[2][netabins];
    comp QmixnAC_3_eta[2][netabins];

    comp Q1nAnorm_eta[2][netabins];
    comp QmixnACnorm_eta[2][netabins];
    comp QmixnACnorm_3_eta[2][netabins];

    double w1nAsum_eta[2][netabins];
    double wmixnACsum_eta[2][netabins];
    double wmixnACsum_3_eta[2][netabins];

    int ncount_eta[2][netabins];
    //--

    int nevts = 0;
    int ncount[2];

    for (int iside = 0; iside<2; iside++) {
        Q1nA[iside] = zero;
        Q1AB[iside] = zero;
        Q1AC[iside] = zero;
        Q1BC[iside] = zero;
        Q2nA[iside] = zero;
        Q2AB[iside] = zero;
        Q2AC[iside] = zero;
        Q2BC[iside] = zero;
        Q3nA[iside] = zero;
        Q3AB[iside] = zero;
        Q3AC[iside] = zero;
        Q3BC[iside] = zero;
        QmixnAC[iside] = zero;
        QmixABC[iside] = zero;
        QmixnAC_3[iside] = zero;
        QmixABC_3[iside] = zero;

        Q1nAnorm[iside] = zero;
        Q1ABnorm[iside] = zero;
        Q1ACnorm[iside] = zero;
        Q1BCnorm[iside] = zero;
        Q2nAnorm[iside] = zero;
        Q2ABnorm[iside] = zero;
        Q2ACnorm[iside] = zero;
        Q2BCnorm[iside] = zero;
        Q3nAnorm[iside] = zero;
        Q3ABnorm[iside] = zero;
        Q3ACnorm[iside] = zero;
        Q3BCnorm[iside] = zero;
        QmixnACnorm[iside] = zero;
        QmixABCnorm[iside] = zero;
        QmixnACnorm_3[iside] = zero;
        QmixABCnorm_3[iside] = zero;

        w1nAsum[iside] = 0;
        w1ABsum[iside] = 0;
        w1ACsum[iside] = 0;
        w1BCsum[iside] = 0;
        w2nAsum[iside] = 0;
        w2ABsum[iside] = 0;
        w2ACsum[iside] = 0;
        w2BCsum[iside] = 0;
        w3nAsum[iside] = 0;
        w3ABsum[iside] = 0;
        w3ACsum[iside] = 0;
        w3BCsum[iside] = 0;
        wmixnACsum[iside] = 0;
        wmixABCsum[iside] = 0;
        wmixnACsum_3[iside] = 0;
        wmixABCsum_3[iside] = 0;

        ncount[iside] = 0;

        for (int pbin = 0; pbin<nptbins; pbin++) {
            Q1nA_pt[iside][pbin] = zero;
            Q2nA_pt[iside][pbin] = zero;
            Q3nA_pt[iside][pbin] = zero;
            QmixnAC_pt[iside][pbin] = zero;
            QmixnAC_3_pt[iside][pbin] = zero;

            Q1nAnorm_pt[iside][pbin] = zero;
            Q2nAnorm_pt[iside][pbin] = zero;
            Q3nAnorm_pt[iside][pbin] = zero;
            QmixnACnorm_pt[iside][pbin] = zero;
            QmixnACnorm_3_pt[iside][pbin] = zero;

            w1nAsum_pt[iside][pbin] = 0;
            w2nAsum_pt[iside][pbin] = 0;
            w3nAsum_pt[iside][pbin] = 0;
            wmixnACsum_pt[iside][pbin] = 0;
            wmixnACsum_3_pt[iside][pbin] = 0;

            ncount_pt[iside][pbin] = 0;
        }
        for (int ebin = 0; ebin<netabins; ebin++) {
            Q1nA_eta[iside][ebin] = zero;
            QmixnAC_eta[iside][ebin] = zero;
            QmixnAC_3_eta[iside][ebin] = zero;

            Q1nAnorm_eta[iside][ebin] = zero;
            QmixnACnorm_eta[iside][ebin] = zero;
            QmixnACnorm_3_eta[iside][ebin] = zero;

            w1nAsum_eta[iside][ebin] = 0;
            wmixnACsum_eta[iside][ebin] = 0;
            wmixnACsum_3_eta[iside][ebin] = 0;

            ncount_eta[iside][ebin] = 0;
        }
    }

    Double_t phiTracks_init[nv1Ebins][multMax];
    Double_t etaTracks_init[nv1Ebins][multMax];
    Double_t ptTracks_init[nv1Ebins][multMax];

    Double_t phiTracks_merge[multMax];
    Double_t etaTracks_merge[multMax];
    Double_t ptTracks_merge[multMax];

    Double_t phiTracks_sub[numEP][multMax];
    Double_t etaTracks_sub[numEP][multMax];
    Double_t ptTracks_sub[numEP][multMax];


    //-- get recentering values
    for (int nep = 0; nep<numEP; nep++) {
        X1[nep] = 0;
        Y1[nep] = 0;
        X2[nep] = 0;
        Y2[nep] = 0;
        X3[nep] = 0;
        Y3[nep] = 0;
    }
    if (recenter) {
        cout << "Getting recentering values... " << endl;
        for (int vbin = 0; vbin<nv1Ebins; vbin++) {
            double v1tmp;
            if (v1odd) v1tmp = v1oddin[vbin];
            else v1tmp = setv1;
            int multtmp = v1Emult[vbin];
            event[vbin]->SetSeed(iseed+cindx);
            event[vbin]->Setv1(v1tmp);
            event[vbin]->SetMult(multtmp);
            event[vbin]->SetEventParms();
        }
        // merge throws into single event
        for (Int_t ievent = 0; ievent<nevents; ievent++) {
            event[0]->SetPsiRandom();
            Double_t PsiRP = event[0]->GetPsi();
            for (int vbin = 0; vbin<nv1Ebins; vbin++) {
                event[vbin]->SetPsi(PsiRP);
            }
            double ptx = 0;
            double pty = 0;
            Int_t jcnt = 0;
            for (int vbin = 0; vbin<nv1Ebins; vbin++) {
                double vmin = v1Ebins[vbin];
                double vmax = v1Ebins[vbin+1];
                event[vbin]->SetPsi(PsiRP);
                event[vbin]->GetThrowPhi(phiTracks_init[vbin]);
                event[vbin]->GetEtaRandom(etaTracks_init[vbin],vmin,vmax);
                event[vbin]->GetPtDist(ptTracks_init[vbin]);
                // event[vbin]->GetPtRandom(ptTracks_init[vbin]);
                for (Int_t j = 0; j<event[vbin]->GetMult(); j++) {
                    double ph = bounds(1,phiTracks_init[vbin][j]);
                    double eta = etaTracks_init[vbin][j];
                    double pt = ptTracks_init[vbin][j];
                    phiTracks_merge[jcnt] = ph;
                    etaTracks_merge[jcnt] = eta;
                    ptTracks_merge[jcnt] = pt;

                    ptx += pt*TMath::Cos(ph);
                    pty += pt*TMath::Sin(ph);
                    jcnt++;
                }
            }

            // conserve momentum
            if (conserve_pT) {
                double ptx_shift = 0;
                double pty_shift = 0;
                double ptave_shift = 0;
                double ptave2_shift = 0;
                for (Int_t j = 0; j<evtmult; j++) {
                    double ph = bounds(1,phiTracks_merge[j]);
                    double pt = ptTracks_merge[j];

                    double xrecoil = pt*TMath::Cos(ph) - ptx/evtmult;
                    double yrecoil = pt*TMath::Sin(ph) - pty/evtmult;
                    double ph_shift = TMath::ATan2(yrecoil,xrecoil);
                    double pt_shift = sqrt(pow(xrecoil,2) + pow(yrecoil,2));

                    phiTracks_merge[j] = ph_shift;
                    ptTracks_merge[j] = pt_shift;
                    ptx_shift+=xrecoil;
                    pty_shift+=yrecoil;
                    ptave_shift+=ptave_shift;
                    ptave2_shift+=pow(ptave_shift,2);
                }
                ptave_shift/=evtmult;
                ptave2_shift/=evtmult;
            }

            // divide into subevents
            for (int nep = 0; nep<numEP; nep++) {
                c1[nep] = 0;
                s1[nep] = 0;
                c2[nep] = 0;
                s2[nep] = 0;
                c3[nep] = 0;
                s3[nep] = 0;
            }
            for (Int_t j = 0; j<evtmult; j++) {
                double ph = bounds(1,phiTracks_merge[j]);
                double eta = etaTracks_merge[j];
                double pt = ptTracks_merge[j];
                for (int nep = 0; nep<numEP-1; nep++) {
                    if (addHoles && eta>=etaMinHole[nep] && eta<etaMaxHole[nep] && ph>=phiMinHole[nep] && ph<phiMaxHole[nep]) pt = 0;
                    if (eta>=ecutmin[nep] && eta<ecutmax[nep] && pt>=pcutmin[nep] && pt<pcutmax[nep]) {
                        c1[nep] += TMath::Cos(ph);
                        s1[nep] += TMath::Sin(ph);
                        c2[nep] += TMath::Cos(2*ph);
                        s2[nep] += TMath::Sin(2*ph);
                        c3[nep] += TMath::Cos(3*ph);
                        s3[nep] += TMath::Sin(3*ph);
                    }
                }
            }
            for (int nep = 0; nep<numEP-1; nep++) {
                X1[nep] += c1[nep];
                Y1[nep] += s1[nep];
                X2[nep] += c2[nep];
                Y2[nep] += s2[nep];
                X3[nep] += c3[nep];
                Y3[nep] += s3[nep];
            }
        }
        for (int nep = 0; nep<numEP-1; nep++) {
            X1[nep]/=nevents;
            Y1[nep]/=nevents;
            X2[nep]/=nevents;
            Y2[nep]/=nevents;
            X3[nep]/=nevents;
            Y3[nep]/=nevents;
            cout << "offsets: "
            << EPName[nep].data() << "\t" << X1[nep] << "\t" << Y1[nep] << "\t" << X2[nep] << "\t" << Y2[nep] <<  "\t" << X3[nep] << "\t" << Y3[nep] << endl;
        }
        cout << " ...recentering complete \n" << endl;
    }


    //-- event plane flattening
    if (flatten) {
        cout << "Generating flattening parameters... " << endl;
        for (int vbin = 0; vbin<nv1Ebins; vbin++) {
            double v1tmp;
            if (v1odd) v1tmp = v1oddin[vbin];
            else v1tmp = setv1;
            int multtmp = v1Emult[vbin];
            event[vbin]->SetSeed(iseed+cindx);
            event[vbin]->Setv1(v1tmp);
            event[vbin]->SetMult(multtmp);
            event[vbin]->SetEventParms();
        }
        // merge throws into single event
        for (Int_t ievent = 0; ievent<nevents; ievent++) {
            event[0]->SetPsiRandom();
            Double_t PsiRP = event[0]->GetPsi();
            for (int vbin = 0; vbin<nv1Ebins; vbin++) {
                event[vbin]->SetPsi(PsiRP);
            }
            double ptx = 0;
            double pty = 0;
            Int_t jcnt = 0;
            for (int vbin = 0; vbin<nv1Ebins; vbin++) {
                double vmin = v1Ebins[vbin];
                double vmax = v1Ebins[vbin+1];
                event[vbin]->SetPsi(PsiRP);
                event[vbin]->GetThrowPhi(phiTracks_init[vbin]);
                event[vbin]->GetEtaRandom(etaTracks_init[vbin],vmin,vmax);
                event[vbin]->GetPtDist(ptTracks_init[vbin]);
                // event[vbin]->GetPtRandom(ptTracks_init[vbin]);
                for (Int_t j = 0; j<event[vbin]->GetMult(); j++) {
                    double ph = bounds(1,phiTracks_init[vbin][j]);
                    double eta = etaTracks_init[vbin][j];
                    double pt = ptTracks_init[vbin][j];
                    phiTracks_merge[jcnt] = ph;
                    etaTracks_merge[jcnt] = eta;
                    ptTracks_merge[jcnt] = pt;

                    ptx += pt*TMath::Cos(ph);
                    pty += pt*TMath::Sin(ph);
                    jcnt++;
                }
            }

            // conserve momentum
            if (conserve_pT) {
                double ptx_shift = 0;
                double pty_shift = 0;
                double ptave_shift = 0;
                double ptave2_shift = 0;
                for (Int_t j = 0; j<evtmult; j++) {
                    double ph = bounds(1,phiTracks_merge[j]);
                    double pt = ptTracks_merge[j];

                    double xrecoil = pt*TMath::Cos(ph) - ptx/evtmult;
                    double yrecoil = pt*TMath::Sin(ph) - pty/evtmult;
                    double ph_shift = TMath::ATan2(yrecoil,xrecoil);
                    double pt_shift = sqrt(pow(xrecoil,2) + pow(yrecoil,2));

                    phiTracks_merge[j] = ph_shift;
                    ptTracks_merge[j] = pt_shift;
                    ptx_shift+=xrecoil;
                    pty_shift+=yrecoil;
                    ptave_shift+=ptave_shift;
                    ptave2_shift+=pow(ptave_shift,2);
                }
                ptave_shift/=evtmult;
                ptave2_shift/=evtmult;
            }

            // divide into subevents
            for (int nep = 0; nep<numEP; nep++) {
                c1[nep] = 0;
                s1[nep] = 0;
                c2[nep] = 0;
                s2[nep] = 0;
                c3[nep] = 0;
                s3[nep] = 0;
            }
            for (Int_t j = 0; j<evtmult; j++) {
                double ph = bounds(1,phiTracks_merge[j]);
                double eta = etaTracks_merge[j];
                double pt = ptTracks_merge[j];
                for (int nep = 0; nep<numEP-1; nep++) {
                    if (addHoles && eta>=etaMinHole[nep] && eta<etaMaxHole[nep] && ph>=phiMinHole[nep] && ph<phiMaxHole[nep]) pt = 0;
                    if (eta>=ecutmin[nep] && eta<ecutmax[nep] && pt>=pcutmin[nep] && pt<pcutmax[nep]) {
                        c1[nep] += TMath::Cos(ph);
                        s1[nep] += TMath::Sin(ph);
                        c2[nep] += TMath::Cos(2*ph);
                        s2[nep] += TMath::Sin(2*ph);
                        c3[nep] += TMath::Cos(2*ph);
                        s3[nep] += TMath::Sin(2*ph);
                    }
                }
            }
            for (int nep = 0; nep<numEP-1; nep++) {
                Psi1[nep] = TMath::ATan2(s1[nep]-Y1[nep],c1[nep]-X1[nep]);
                Psi2[nep] = TMath::ATan2(s2[nep]-Y2[nep],c2[nep]-X2[nep]);
                Psi3[nep] = TMath::ATan2(s3[nep]-Y3[nep],c3[nep]-X3[nep]);
                flat1[nep]->Fill(Psi1[nep],0,0);
                flat2[nep]->Fill(Psi2[nep],0,0);
                flat3[nep]->Fill(Psi3[nep],0,0);
            }
        }
        for (int nep = 0; nep<numEP; nep++) {
            for (int j = 0; j<flat1[nep]->GetHBins(); j++) {
                if (flat1[nep]->GetCnt(j)>0) {
                    flat1[nep]->SetXDB(j,flat1[nep]->GetX(j)/flat1[nep]->GetCnt(j));
                    flat1[nep]->SetXDB(j,flat1[nep]->GetY(j)/flat1[nep]->GetCnt(j));
                }
            }
            for (int j = 0; j<flat2[nep]->GetHBins(); j++) {
                if (flat2[nep]->GetCnt(j)>0) {
                    flat2[nep]->SetXDB(j,flat2[nep]->GetX(j)/flat2[nep]->GetCnt(j));
                    flat2[nep]->SetXDB(j,flat2[nep]->GetY(j)/flat2[nep]->GetCnt(j));
                }
            }
            for (int j = 0; j<flat3[nep]->GetHBins(); j++) {
                if (flat3[nep]->GetCnt(j)>0) {
                    flat3[nep]->SetXDB(j,flat3[nep]->GetX(j)/flat3[nep]->GetCnt(j));
                    flat3[nep]->SetXDB(j,flat3[nep]->GetY(j)/flat3[nep]->GetCnt(j));
                }
            }
        }
        cout << " ...flattening done \n" << endl;
    }


    //-- main event loop

    cout << "Entering primary event loop" << endl;
    for (int vbin = 0; vbin<nv1Ebins; vbin++) {
        double v1tmp;
        if (v1odd) v1tmp = v1oddin[vbin];
        else v1tmp = setv1;
        int multtmp = v1Emult[vbin];
        event[vbin]->SetSeed(iseed+cindx);
        event[vbin]->Setv1(v1tmp);
        event[vbin]->SetMult(multtmp);
        event[vbin]->SetEventParms();
        hinit_v1in->SetBinContent(vbin+1,v1tmp);
        hinit_v1in->SetBinError(vbin+1,0.00001);
    }

    //-- make initial throws and merge into single event
    for (Int_t ievent = 0; ievent<nevents; ievent++) {
        if (fmod(double(ievent+1), nevents/20.)==0) cout << " event: " << ievent+1
            << "/" << nevents << "\trun: " << cindx+1 << "/" << Ntries <<  endl;

        event[0]->SetPsiRandom();
        Double_t PsiRP = event[0]->GetPsi();
//        Double_t PsiRP = ran->Uniform(-TMath::Pi(), TMath::Pi());
        double ptx = 0;
        double pty = 0;
        double ptave = 0;
        double ptave2 = 0;
        Int_t jcnt = 0;
        for (int vbin = 0; vbin<nv1Ebins; vbin++) {
            double vmin = v1Ebins[vbin];
            double vmax = v1Ebins[vbin+1];
            event[vbin]->SetPsi(PsiRP);
            event[vbin]->GetThrowPhi(phiTracks_init[vbin]);
            event[vbin]->GetEtaRandom(etaTracks_init[vbin],vmin,vmax);
            event[vbin]->GetPtDist(ptTracks_init[vbin]);
            // event[vbin]->GetPtRandom(ptTracks_init[vbin]);
            for (Int_t j = 0; j<event[vbin]->GetMult(); j++) {
                double ph = bounds(1,phiTracks_init[vbin][j]);
                double eta = etaTracks_init[vbin][j];
                double pt = ptTracks_init[vbin][j];
                phiTracks_merge[jcnt] = ph;
                etaTracks_merge[jcnt] = eta;
                ptTracks_merge[jcnt] = pt;

                ptx += pt*TMath::Cos(ph);
                pty += pt*TMath::Sin(ph);
                ptave += pt;
                ptave2 += pt*pt;

                hinit_philab[vbin]->Fill(bounds(1,ph));
                hinit_eta[vbin]->Fill(eta);
                hinit_pt[vbin]->Fill(pt);
                hinit_phiPsiRP[vbin]->Fill(bounds(1,ph - PsiRP));
                hetain->Fill(eta);
                if (!conserve_pT) {
                    hphilab->Fill(bounds(1,ph));
                    hptin->Fill(pt);
                    hphiPsiRP->Fill(bounds(1,ph - PsiRP));
                    hptAve2->Fill(pow(pt,2));
                }

                jcnt++;
            }
        }
        if (!conserve_pT) hpt2D->Fill(ptx,pty);


        //-- apply momentum conserving shift
        if (conserve_pT) {
            double ptx_shift = 0;
            double pty_shift = 0;
            double ptave_shift = 0;
            double ptave2_shift = 0;
            for (Int_t j = 0; j<evtmult; j++) {
                double ph = bounds(1,phiTracks_merge[j]);
                double pt = ptTracks_merge[j];

                double xrecoil = pt*TMath::Cos(ph) - ptx/evtmult;
                double yrecoil = pt*TMath::Sin(ph) - pty/evtmult;
                double ph_shift = TMath::ATan2(yrecoil,xrecoil);
                double pt_shift = sqrt(pow(xrecoil,2) + pow(yrecoil,2));

                phiTracks_merge[j] = ph_shift;
                ptTracks_merge[j] = pt_shift;
                ptx_shift+=xrecoil;
                pty_shift+=yrecoil;
                ptave_shift+=ptave_shift;
                ptave2_shift+=pow(ptave_shift,2);

                hphilab->Fill(bounds(1,ph_shift));
                hptin->Fill(pt_shift);
                hphiPsiRP->Fill(bounds(1,ph_shift - PsiRP));
                hptAve2->Fill(pow(pt_shift,2));
            }
            ptave_shift/=evtmult;
            ptave2_shift/=evtmult;
            hpt2D->Fill(ptx_shift,pty_shift);
        }


        //-- sort particles into subevents
        for (int nep = 0; nep<numEP; nep++) {
            subptave[nep] = 0;
            subptave2[nep] = 0;
            subcnt[nep] = 0;
        }
        for (Int_t j = 0; j<evtmult; j++) {
            double ph = bounds(1,phiTracks_merge[j]);
            double eta = etaTracks_merge[j];
            double pt = ptTracks_merge[j];
            // HF-, track-, track+, HF+, and trackmid
            bool isinsc = 0;
            for (int nep = 0; nep<numEP-1; nep++) {
                if (addHoles && eta>=etaMinHole[nep] && eta<etaMaxHole[nep] && ph>=phiMinHole[nep] && ph<phiMaxHole[nep]) pt = 0;
                if (eta>=ecutmin[nep] && eta<ecutmax[nep] && pt>=pcutmin[nep] && pt<pcutmax[nep]) {
                    isinsc = 1;
                    phiTracks_sub[nep][subcnt[nep]] = ph;
                    etaTracks_sub[nep][subcnt[nep]] = eta;
                    ptTracks_sub[nep][subcnt[nep]] = pt;
                    subptave[nep] += pt;
                    subptave2[nep] += pow(pt,2);
                    ++subcnt[nep];
                }
            }
            // particles not accepted
            if (!isinsc) {
                phiTracks_sub[numEP-1][subcnt[numEP-1]] = ph;
                etaTracks_sub[numEP-1][subcnt[numEP-1]] = eta;
                ptTracks_sub[numEP-1][subcnt[numEP-1]] = pt;
                subptave[numEP-1] += pt;
                subptave2[numEP-1] += pow(pt,2);
                ++subcnt[numEP-1];
            }
        }

        //-- calculate event plane angles
        for (int nep = 0; nep<numEP; nep++) {
            subptave[nep]/=subcnt[nep];
            subptave2[nep]/=subcnt[nep];
            hsub_ptAve[nep]->Fill(subptave[nep]);
            hsub_ptAve2[nep]->Fill(subptave2[nep]);

            q1x[nep] = 0;
            q1y[nep] = 0;
            q2x[nep] = 0;
            q2y[nep] = 0;
            q3x[nep] = 0;
            q3y[nep] = 0;
            Psi1[nep] = 0;
            Psi2[nep] = 0;
            Psi3[nep] = 0;
            w1[nep] = 0;
            w2[nep] = 0;
            w3[nep] = 0;
            qcnt[nep] = 0;
            double ptxsub = 0;
            double ptysub = 0;

            for (Int_t j = 0; j<subcnt[nep]; j++) {
                double ph = bounds(1,phiTracks_sub[nep][j]);
                double eta = etaTracks_sub[nep][j];
                double pt = ptTracks_sub[nep][j];
                double ww1 = 1;
                double ww2 = ww1;
                double ww3 = ww1;
                if (pt_weights) ww1 = pt - (subptave2[nep]/subptave[nep]);
                if (eta_weights && eta<0) ww1*=-1;
                ptxsub += pt*TMath::Cos(ph);
                ptysub += pt*TMath::Sin(ph);

                hsub_philab[nep]->Fill(bounds(1,ph));
                hsub_etain[nep]->Fill(eta);
                hsub_ptin[nep]->Fill(pt);
                hsub_phiPsiRP[nep]->Fill(bounds(1,ph - PsiRP));
                hweights->Fill(ww1);
                hwpt->Fill(ww1*pt);

                q1x[nep] += ww1*TMath::Cos(ph);
                q1y[nep] += ww1*TMath::Sin(ph);
                q2x[nep] += ww2*TMath::Cos(2*ph);
                q2y[nep] += ww2*TMath::Sin(2*ph);
                q3x[nep] += ww3*TMath::Cos(3*ph);
                q3y[nep] += ww3*TMath::Sin(3*ph);
                w1[nep]  += ww1;
                w2[nep]  += ww2;
                w3[nep]  += ww3;
                qcnt[nep]++;
            }
            Psi1[nep] = TMath::ATan2(q1y[nep]-Y1[nep],q1x[nep]-X1[nep]);
            Psi2[nep] = TMath::ATan2(q2y[nep]-Y2[nep],q2x[nep]-X2[nep])/2.;
            Psi3[nep] = TMath::ATan2(q3y[nep]-Y3[nep],q3x[nep]-X3[nep])/3.;
            if (flatten) {
                flat1[nep]->Fill(Psi1[nep],0,0);
                flat2[nep]->Fill(Psi2[nep],0,0);
                flat3[nep]->Fill(Psi3[nep],0,0);
                Double_t ps1 = Psi1[nep];
                Double_t ps2 = Psi2[nep];
                Double_t ps3 = Psi3[nep];
                Double_t ps1flat = flat1[nep]->GetFlatPsi(ps1,0,0);
                Double_t ps2flat = flat2[nep]->GetFlatPsi(ps2,0,0);
                Double_t ps3flat = flat3[nep]->GetFlatPsi(ps3,0,0);
                Psi1[nep] = ps1flat;
                Psi2[nep] = ps2flat;
                Psi3[nep] = ps3flat;
            }

            EPCorr(Psi1[0], Psi1[1], Psi1[2], Psi1[3], Psi1[4], Psi2[0], Psi2[1], Psi2[2], Psi2[3], Psi2[4]);

            hsub_Psi1[nep]->Fill(bounds(1,Psi1[nep] - PsiRP));
            hsub_Psi2[nep]->Fill(bounds(2,Psi2[nep] - PsiRP));
            hsub_Psi3[nep]->Fill(bounds(3,Psi3[nep] - PsiRP));
            hsub_Psi1lab[nep]->Fill(bounds(1,Psi1[nep]));
            hsub_Psi2lab[nep]->Fill(bounds(2,Psi2[nep]));
            hsub_Psi3lab[nep]->Fill(bounds(3,Psi3[nep]));
            hsub_pt2D[nep]->Fill(ptxsub,ptysub);
            for (Int_t j = 0; j<subcnt[nep]; j++) {
                double ph = phiTracks_sub[nep][j];
                hsub_phiPsi1[nep]->Fill(bounds(1,ph - Psi1[nep]));
                hsub_phiPsi2[nep]->Fill(bounds(2,ph - Psi2[nep]));
                hsub_phiPsi3[nep]->Fill(bounds(3,ph - Psi3[nep]));
            }
        }

        //-- calculate vn
        int ep1A, ep1B, ep1C;
        int ep2A, ep2B, ep2C;
        int ep3A, ep3B, ep3C;
        int epPOI;
        for (int iside = 0; iside<2; iside++) {
            if (iside == 0) {
                ep1A = 0; // Psi1A: HFm1
                ep1B = 3; // Psi1B: HFp1
                ep1C = 1; // Psi1C: trackm1
                ep2A = 0;
                ep2B = 3;
                ep2C = 1;
                ep3A = 0;
                ep3B = 3;
                ep3C = 1;
                epPOI = 2;
            }
            // iside = 1: A = HFp (3), B = HFm (0), C = trackp (2), POI = trackm (1)
            if (iside == 1) {
                ep1A = 3;
                ep1B = 0;
                ep1C = 2;
                ep2A = 3;
                ep2B = 0;
                ep2C = 2;
                ep3A = 3;
                ep3B = 0;
                ep3C = 2;
                epPOI = 1;
            }

            w1nAsum[iside] += w1[ep1A] * (double)qcnt[epPOI];
            w1ABsum[iside] += w1[ep1A] * w1[ep1B];
            w1ACsum[iside] += w1[ep1A] * w1[ep1C];
            w1BCsum[iside] += w1[ep1B] * w1[ep1C];
            w2nAsum[iside] += w2[ep2A] * (double)qcnt[epPOI];
            w2ABsum[iside] += w2[ep2A] * w2[ep2B];
            w2ACsum[iside] += w2[ep2A] * w2[ep2C];
            w2BCsum[iside] += w2[ep2B] * w2[ep2C];

            w3nAsum[iside] += w3[ep3A] * (double)qcnt[epPOI];
            w3ABsum[iside] += w3[ep3A] * w3[ep3B];
            w3ACsum[iside] += w3[ep3A] * w3[ep3C];
            w3BCsum[iside] += w3[ep3B] * w3[ep3C];

            wmixnACsum[iside] += w1[ep1A] * w2[ep2C] * (double)qcnt[epPOI];
            wmixABCsum[iside] += w1[ep1A] * w1[ep1B] * w2[ep2C];
            ////
            wmixnACsum_3[iside] += w2[ep2A] * w3[ep3C] * (double)qcnt[epPOI];
            wmixABCsum_3[iside] += w1[ep1B] * w2[ep2A] * w3[ep3C];
            ////

            comp Q1n(q1x[epPOI], q1y[epPOI]);
            comp Q1A(q1x[ep1A], q1y[ep1A]);
            comp Q1B(q1x[ep1B], q1y[ep1B]);
            comp Q1C(q1x[ep1C], q1y[ep1C]);

            comp Q2n(q2x[epPOI], q2y[epPOI]);
            comp Q2A(q2x[ep2A], q2y[ep2A]);
            comp Q2B(q2x[ep2B], q2y[ep2B]);
            comp Q2C(q2x[ep2C], q2y[ep2C]);

            comp Q3n(q3x[epPOI], q3y[epPOI]);
            comp Q3A(q3x[ep3A], q3y[ep3A]);
            comp Q3B(q3x[ep3B], q3y[ep3B]);
            comp Q3C(q3x[ep3C], q3y[ep3C]);

            comp nA1 = Q1n * std::conj(Q1A);
            comp AB1 = Q1A * std::conj(Q1B);
            comp AC1 = Q1A * std::conj(Q1C);
            comp BC1 = Q1B * std::conj(Q1C);
            comp nA2 = Q2n * std::conj(Q2A);
            comp AB2 = Q2A * std::conj(Q2B);
            comp AC2 = Q2A * std::conj(Q2C);
            comp BC2 = Q2B * std::conj(Q2C);
            comp nA3 = Q3n * std::conj(Q3A);
            comp AB3 = Q3A * std::conj(Q3B);
            comp AC3 = Q3A * std::conj(Q3C);
            comp BC3 = Q3B * std::conj(Q3C);
            comp mixnAC = Q1n * Q1A * std::conj(Q2C);
            comp mixABC = Q1A * Q1B * std::conj(Q2C);
            ////
            comp mixnAC_3 = Q1n * Q2A * std::conj(Q3C);
            comp mixABC_3 = Q1B * Q2A * std::conj(Q3C);
            ////

            Q1nA[iside] += nA1;
            Q1AB[iside] += AB1;
            Q1AC[iside] += AC1;
            Q1BC[iside] += BC1;
            Q2nA[iside] += nA2;
            Q2AB[iside] += AB2;
            Q2AC[iside] += AC2;
            Q2BC[iside] += BC2;
            Q3nA[iside] += nA3;
            Q3AB[iside] += AB3;
            Q3AC[iside] += AC3;
            Q3BC[iside] += BC3;
            QmixnAC[iside] += mixnAC;
            QmixABC[iside] += mixABC;
            QmixnAC_3[iside] += mixnAC_3;
            QmixABC_3[iside] += mixABC_3;

            comp nA1norm = Q1n * (std::conj(Q1A)/std::abs(Q1A));
            comp AB1norm = (Q1A/std::abs(Q1A)) * (std::conj(Q1B)/std::abs(Q1B));
            comp AC1norm = (Q1A/std::abs(Q1A)) * (std::conj(Q1C)/std::abs(Q1C));
            comp BC1norm = (Q1B/std::abs(Q1B)) * (std::conj(Q1C)/std::abs(Q1C));
            comp nA2norm = Q2n * (std::conj(Q2A)/std::abs(Q2A));
            comp AB2norm = (Q2A/std::abs(Q2A)) * (std::conj(Q2B)/std::abs(Q2B));
            comp AC2norm = (Q2A/std::abs(Q2A)) * (std::conj(Q2C)/std::abs(Q2C));
            comp BC2norm = (Q2B/std::abs(Q2B)) * (std::conj(Q2C)/std::abs(Q2C));
            comp nA3norm = Q3n * (std::conj(Q3A)/std::abs(Q3A));
            comp AB3norm = (Q3A/std::abs(Q3A)) * (std::conj(Q3B)/std::abs(Q3B));
            comp AC3norm = (Q3A/std::abs(Q3A)) * (std::conj(Q3C)/std::abs(Q3C));
            comp BC3norm = (Q3B/std::abs(Q3B)) * (std::conj(Q3C)/std::abs(Q3C));
            comp mixnACnorm = Q1n * Q1A/std::abs(Q1A) * (std::conj(Q2C)/std::abs(Q2C));
            comp mixABCnorm = (Q1A/std::abs(Q1A)) * (Q1B/std::abs(Q1B)) * (std::conj(Q2C)/std::abs(Q2C));
            ////
            comp mixnACnorm_3 = Q1n * Q2A/std::abs(Q2A) * (std::conj(Q3C)/std::abs(Q3C));
            comp mixABCnorm_3 = (Q1B/std::abs(Q1B)) * (Q2A/std::abs(Q2A)) * (std::conj(Q3C)/std::abs(Q3C));
            ////

            Q1nAnorm[iside] += nA1norm;
            Q1ABnorm[iside] += AB1norm;
            Q1ACnorm[iside] += AC1norm;
            Q1BCnorm[iside] += BC1norm;
            Q2nAnorm[iside] += nA2norm;
            Q2ABnorm[iside] += AB2norm;
            Q2ACnorm[iside] += AC2norm;
            Q2BCnorm[iside] += BC2norm;
            Q3nAnorm[iside] += nA3norm;
            Q3ABnorm[iside] += AB3norm;
            Q3ACnorm[iside] += AC3norm;
            Q3BCnorm[iside] += BC3norm;
            QmixnACnorm[iside] += mixnACnorm;
            QmixABCnorm[iside] += mixABCnorm;
            QmixnACnorm_3[iside] += mixnACnorm_3;
            QmixABCnorm_3[iside] += mixABCnorm_3;

            ncount[iside]+=(double)qcnt[epPOI];


            //-- differential vn
            for (int pbin = 0; pbin<nptbins; pbin++) {
                q1x_pt[iside][pbin] = 0;
                q1y_pt[iside][pbin] = 0;
                q2x_pt[iside][pbin] = 0;
                q2y_pt[iside][pbin] = 0;
                q3x_pt[iside][pbin] = 0;
                q3y_pt[iside][pbin] = 0;
                w1_pt[iside][pbin] = 0;
                w2_pt[iside][pbin] = 0;
                w3_pt[iside][pbin] = 0;
                qcnt_pt[iside][pbin] = 0;
            }
            for (int ebin = 0; ebin<netabins; ebin++) {
                q1x_eta[iside][ebin] = 0;
                q1y_eta[iside][ebin] = 0;
                w1_eta[iside][ebin] = 0;
                qcnt_eta[iside][ebin] = 0;
            }
            for (int j = 0; j<evtmult; j++) {
                double ph = bounds(1,phiTracks_merge[j]);
                double eta = etaTracks_merge[j];
                double pt = ptTracks_merge[j];
                double ww1 = 1;
                double ww2 = ww1;
                double ww3 = ww1;
                if (pt_weights) ww1 = pt - (subptave2[1]+subptave2[2])/(subptave[1]+subptave[2]);
                if (eta_weights && iside == 0) ww1*=-1;
                for (int pbin = 0; pbin<nptbins; pbin++) {
                    if (pt>=ptbins[pbin] && pt<ptbins[pbin+1] && fabs(eta)<2.4) {
                        q1x_pt[iside][pbin] += ww1*TMath::Cos(ph);
                        q1y_pt[iside][pbin] += ww1*TMath::Sin(ph);
                        q2x_pt[iside][pbin] += ww2*TMath::Cos(2*ph);
                        q2y_pt[iside][pbin] += ww2*TMath::Sin(2*ph);
                        q3x_pt[iside][pbin] += ww3*TMath::Cos(3*ph);
                        q3y_pt[iside][pbin] += ww3*TMath::Sin(3*ph);
                        w1_pt[iside][pbin] += ww1;
                        w2_pt[iside][pbin] += ww2;
                        w3_pt[iside][pbin] += ww3;
                        qcnt_pt[iside][pbin]++;
                    }
                }
                for (int ebin = 0; ebin<netabins; ebin++) {
                    if (eta>=etabins[ebin] && eta<etabins[ebin+1] && pt>0.3 && pt<3.0) {
                        q1x_eta[iside][ebin] += ww1*TMath::Cos(ph);
                        q1y_eta[iside][ebin] += ww1*TMath::Sin(ph);
                        w1_eta[iside][ebin] += ww1;
                        qcnt_eta[iside][ebin]++;
                    }
                }
            }
            for (int pbin = 0; pbin<nptbins; pbin++) {
                w1nAsum_pt[iside][pbin] += w1[ep1A] * (double)qcnt_pt[iside][pbin];
                w2nAsum_pt[iside][pbin] += w2[ep2A] * (double)qcnt_pt[iside][pbin];
                w3nAsum_pt[iside][pbin] += w3[ep3A] * (double)qcnt_pt[iside][pbin];
                wmixnACsum_pt[iside][pbin] += w1[ep1A] * w2[ep2C] * (double)qcnt_pt[iside][pbin];
                ////
                wmixnACsum_3_pt[iside][pbin] += w2[ep2A] * w3[ep3C] * (double)qcnt_pt[iside][pbin];
                ////

                comp Q1n_pt(q1x_pt[iside][pbin], q1y_pt[iside][pbin]);
                comp Q2n_pt(q2x_pt[iside][pbin], q2y_pt[iside][pbin]);
                comp Q3n_pt(q2x_pt[iside][pbin], q3y_pt[iside][pbin]);

                comp nA1_pt = Q1n_pt * std::conj(Q1A);
                comp nA2_pt = Q2n_pt * std::conj(Q2A);
                comp nA3_pt = Q3n_pt * std::conj(Q3A);
                comp mixnAC_pt = Q1n_pt * Q1A * std::conj(Q2C);
                ////
                comp mixnAC_3_pt = Q1n_pt * Q2A * std::conj(Q3C);
                ////

                Q1nA_pt[iside][pbin] += nA1_pt;
                Q2nA_pt[iside][pbin] += nA2_pt;
                Q3nA_pt[iside][pbin] += nA3_pt;
                QmixnAC_pt[iside][pbin] += mixnAC_pt;
                QmixnAC_3_pt[iside][pbin] += mixnAC_3_pt;

                comp nA1norm_pt = Q1n_pt * (std::conj(Q1A)/std::abs(Q1A));
                comp nA2norm_pt = Q2n_pt * (std::conj(Q2A)/std::abs(Q2A));
                comp nA3norm_pt = Q3n_pt * (std::conj(Q3A)/std::abs(Q3A));
                comp mixnACnorm_pt = Q1n_pt * Q1A/std::abs(Q1A) * (std::conj(Q2C)/std::abs(Q2C));
                ////
                comp mixnACnorm_3_pt = Q1n_pt * Q2A/std::abs(Q2A) * (std::conj(Q3C)/std::abs(Q3C));
                ////

                Q1nAnorm_pt[iside][pbin] += nA1norm_pt;
                Q2nAnorm_pt[iside][pbin] += nA2norm_pt;
                Q3nAnorm_pt[iside][pbin] += nA3norm_pt;
                QmixnACnorm_pt[iside][pbin] += mixnACnorm_pt;
                QmixnACnorm_3_pt[iside][pbin] += mixnACnorm_3_pt;

                ncount_pt[iside][pbin]+=(double)qcnt_pt[iside][pbin];
                hqcnt_pt[iside][pbin]->Fill(qcnt_pt[iside][pbin]);
            }

            for (int ebin = 0; ebin<netabins; ebin++) {
                w1nAsum_eta[iside][ebin] += w1[ep1A] * (double)qcnt_eta[iside][ebin];
                wmixnACsum_eta[iside][ebin] += w1[ep1A] * w2[ep2C] * (double)qcnt_eta[iside][ebin];
                ////
                wmixnACsum_3_eta[iside][ebin] += w2[ep2A] * w3[ep3C] * (double)qcnt_eta[iside][ebin];
                ////

                comp Q1n_eta(q1x_eta[iside][ebin], q1y_eta[iside][ebin]);

                comp nA1_eta = Q1n_eta * std::conj(Q1A);
                comp mixnAC_eta = Q1n_eta * Q1A * std::conj(Q2C);
                ////
                comp mixnAC_3_eta = Q1n_eta * Q2A * std::conj(Q3C);
                ////

                Q1nA_eta[iside][ebin] += nA1_eta;
                QmixnAC_eta[iside][ebin] += mixnAC_eta;
                QmixnAC_3_eta[iside][ebin] += mixnAC_3_eta;

                comp nA1norm_eta = Q1n_eta * (std::conj(Q1A)/std::abs(Q1A));
                comp mixnACnorm_eta = Q1n_eta * Q1A/std::abs(Q1A) * (std::conj(Q2C)/std::abs(Q2C));
                ////
                comp mixnACnorm_3_eta = Q1n_eta * Q2A/std::abs(Q2A) * (std::conj(Q3C)/std::abs(Q3C));
                ////

                Q1nAnorm_eta[iside][ebin] += nA1norm_eta;
                QmixnACnorm_eta[iside][ebin] += mixnACnorm_eta;
                QmixnACnorm_3_eta[iside][ebin] += mixnACnorm_3_eta;

                ncount_eta[iside][ebin]+=(double)qcnt_eta[iside][ebin];
                hqcnt_eta[iside][ebin]->Fill(qcnt_eta[iside][ebin]);
            }
        }
        ++nevts;
    }

    rescor1_2 = sqrt( fabs(Q1ABnorm[0].real())/(double)nevts );
    rescor1_2SE->Fill(rescor1_2);

    rescor2HF_2 = sqrt( fabs(Q2ABnorm[0].real())/(double)nevts );
    rescor2HF_2SE->Fill(rescor2HF_2);

    rescor3HF_2 = sqrt( fabs(Q3ABnorm[0].real())/(double)nevts );
    rescor3HF_2SE->Fill(rescor3HF_2);

    spdenom1_2SE = sqrt( fabs(Q1AB[0].real()/w1ABsum[0]) );
    SPdenom1_2SE->Fill(spdenom1_2SE);

    spdenom2HF_2SE = sqrt( fabs(Q2AB[0].real()/w2ABsum[0]) );
    SPdenom2HF_2SE->Fill(spdenom2HF_2SE);

    spdenom3HF_2SE = sqrt( fabs(Q3AB[0].real()/w3ABsum[0]) );
    SPdenom3HF_2SE->Fill(spdenom3HF_2SE);


    for (int iside = 0; iside<2; iside++) {
        double nEP1evts = (double)nevts;
        double nEP1count = (double)ncount[iside];
        if (eta_weights && iside == 0) {nEP1evts*=-1; nEP1count*=-1;}
        // This preserves the eta-weighting on the particles of interest
        // This is already accounted for in the scalar-product weights

        double na1ep = Q1nAnorm[iside].real()/nEP1count;
        double ab1ep = Q1ABnorm[iside].real()/nEP1evts;
        double ac1ep = Q1ACnorm[iside].real()/nEP1evts;
        double bc1ep = Q1BCnorm[iside].real()/nEP1evts;

        double nacep = QmixnACnorm[iside].real()/nEP1count;
        double abcep = QmixABCnorm[iside].real()/nEP1evts;
        double nacep_3 = QmixnACnorm_3[iside].real()/nEP1count;
        double abcep_3 = QmixABCnorm_3[iside].real()/nEP1evts;

        double na2ep = Q2nAnorm[iside].real()/(double)ncount[iside];
        double ab2ep = Q2ABnorm[iside].real()/(double)nevts;
        double ac2ep = Q2ACnorm[iside].real()/(double)nevts;
        double bc2ep = Q2BCnorm[iside].real()/(double)nevts;

        double na3ep = Q3nAnorm[iside].real()/(double)ncount[iside];
        double ab3ep = Q3ABnorm[iside].real()/(double)nevts;
        double ac3ep = Q3ACnorm[iside].real()/(double)nevts;
        double bc3ep = Q3BCnorm[iside].real()/(double)nevts;


        double na1sp = Q1nA[iside].real()/w1nAsum[iside];
        double ab1sp = Q1AB[iside].real()/w1ABsum[iside];
        double ac1sp = Q1AC[iside].real()/w1ACsum[iside];
        double bc1sp = Q1BC[iside].real()/w1BCsum[iside];

        double nacsp = QmixnAC[iside].real()/wmixnACsum[iside];
        double abcsp = QmixABC[iside].real()/wmixABCsum[iside];
        double nacsp_3 = QmixnAC_3[iside].real()/wmixnACsum_3[iside];
        double abcsp_3 = QmixABC_3[iside].real()/wmixABCsum_3[iside];

        double na2sp = Q2nA[iside].real()/w2nAsum[iside];
        double ab2sp = Q2AB[iside].real()/w2ABsum[iside];
        double ac2sp = Q2AC[iside].real()/w2ACsum[iside];
        double bc2sp = Q2BC[iside].real()/w2BCsum[iside];

        double na3sp = Q3nA[iside].real()/w3nAsum[iside];
        double ab3sp = Q3AB[iside].real()/w3ABsum[iside];
        double ac3sp = Q3AC[iside].real()/w3ACsum[iside];
        double bc3sp = Q3BC[iside].real()/w3BCsum[iside];

        // Q-vector averages
        hQ1AB_final[iside]->Fill(ab1sp);
        hQ1AC_final[iside]->Fill(ac1sp);
        hQ1BC_final[iside]->Fill(bc1sp);
        hQ2AB_final[iside]->Fill(ab2sp);
        hQ2AC_final[iside]->Fill(ac2sp);
        hQ2BC_final[iside]->Fill(bc2sp);
        hQ3AB_final[iside]->Fill(ab3sp);
        hQ3AC_final[iside]->Fill(ac3sp);
        hQ3BC_final[iside]->Fill(bc3sp);
        hQmixABC_final[iside]->Fill(abcsp);
        hQmixABC_3_final[iside]->Fill(abcsp_3);

        // normalized Q-vector averages
        hQ1ABnorm_final[iside]->Fill(ab1ep);
        hQ1ACnorm_final[iside]->Fill(ac1ep);
        hQ1BCnorm_final[iside]->Fill(bc1ep);
        hQ2ABnorm_final[iside]->Fill(ab2ep);
        hQ2ACnorm_final[iside]->Fill(ac2ep);
        hQ2BCnorm_final[iside]->Fill(bc2ep);
        hQ3ABnorm_final[iside]->Fill(ab3ep);
        hQ3ACnorm_final[iside]->Fill(ac3ep);
        hQ3BCnorm_final[iside]->Fill(bc3ep);
        hQmixABCnorm_final[iside]->Fill(abcep);
        hQmixABC_3norm_final[iside]->Fill(abcep_3);


        rescor1_3[iside] = sqrt( fabs(ab1ep*ac1ep/bc1ep) );
        rescor1_3SE[iside]->Fill(rescor1_3[iside]);

        rescor2HF_3[iside] = sqrt( fabs(ab2ep*ac2ep/bc2ep) );
        rescor2HF_3SE[iside]->Fill(rescor2HF_3[iside]);

        rescor2Trk_3[iside] = sqrt( fabs(bc2ep*ac2ep/ab2ep) );
        rescor2Trk_3SE[iside]->Fill(rescor2Trk_3[iside]);

        rescor3HF_3[iside] = sqrt( fabs(ab3ep*ac3ep/bc3ep) );
        rescor3HF_3SE[iside]->Fill(rescor3HF_3[iside]);

        rescor3Trk_3[iside] = sqrt( fabs(bc3ep*ac3ep/ab3ep) );
        rescor3Trk_3SE[iside]->Fill(rescor3Trk_3[iside]);

        //rescorMix[iside] = sqrt( fabs(abcep * rescor2Trk_3[iside]) );
        rescorMix[iside] = rescor1_3[iside] * rescor2Trk_3[iside];
        rescor1_mix[iside]->Fill(rescorMix[iside]);

        ////
        // rescorMix_3[iside] = sqrt( fabs(rescor1_3[iside] * rescor2Trk_3[iside] * rescor3Trk_3[iside]) );
        rescorMix_3[iside] = rescor2HF_3[iside] * rescor3Trk_3[iside];
        rescor1_mix3[iside]->Fill(rescorMix_3[iside]);
        ////

        epv1obs[iside] = na1ep;
        EPv1obs[iside]->Fill(epv1obs[iside]);

        epv2obs[iside] = na2ep;
        EPv2obs[iside]->Fill(epv2obs[iside]);

        epv3obs[iside] = na3ep;
        EPv3obs[iside]->Fill(epv3obs[iside]);

        epv1_2SE[iside] = na1ep/rescor1_2;
        EPv1_2SE[iside]->Fill(epv1_2SE[iside]);

        epv1_3SE[iside] = na1ep/rescor1_3[iside];
        EPv1_3SE[iside]->Fill(epv1_3SE[iside]);

        epv2_2SE[iside] = na2ep/rescor2HF_2;
        EPv2_2SE[iside]->Fill(epv2_2SE[iside]);

        epv2_3SE[iside] = na2ep/rescor2HF_3[iside];
        EPv2_3SE[iside]->Fill(epv2_3SE[iside]);

        epv3_2SE[iside] = na3ep/rescor3HF_2;
        EPv3_2SE[iside]->Fill(epv3_2SE[iside]);

        epv3_3SE[iside] = na3ep/rescor3HF_3[iside];
        EPv3_3SE[iside]->Fill(epv3_3SE[iside]);

        epv1obs_mix[iside] = nacep;
        EPv1obsMix[iside]->Fill(epv1obs_mix[iside]);

        epv1obs_mix_3[iside] = nacep_3;
        EPv1obsMix_3[iside]->Fill(epv1obs_mix_3[iside]);

        epv1_mix[iside] = nacep/rescorMix[iside];
        EPv1_mix[iside]->Fill(epv1_mix[iside]);

        epv1_mix_3[iside] = nacep_3/rescorMix_3[iside];
        EPv1_mix_3[iside]->Fill(epv1_mix_3[iside]);

        spdenom1_3SE[iside] = sqrt( fabs(ab1sp*ac1sp/bc1sp) );
        SPdenom1_3SE[iside]->Fill(spdenom1_3SE[iside]);

        spdenom2HF_3SE[iside] = sqrt( fabs(ab2sp*ac2sp/bc2sp) );
        SPdenom2HF_3SE[iside]->Fill(spdenom2HF_3SE[iside]);

        spdenom2Trk_3SE[iside] = sqrt( fabs(bc2sp*ac2sp/ab2sp) );
        SPdenom2Trk_3SE[iside]->Fill(spdenom2Trk_3SE[iside]);

        spdenom3HF_3SE[iside] = sqrt( fabs(ab3sp*ac3sp/bc3sp) );
        SPdenom3HF_3SE[iside]->Fill(spdenom3HF_3SE[iside]);

        spdenom3Trk_3SE[iside] = sqrt( fabs(bc3sp*ac3sp/ab3sp) );
        SPdenom3Trk_3SE[iside]->Fill(spdenom3Trk_3SE[iside]);

        //spdenom_mix[iside] = sqrt( fabs(abcsp) * spdenom2Trk_3SE[iside] );
        spdenom_mix[iside] = spdenom1_3SE[iside] * spdenom2Trk_3SE[iside];
        SPdenom1_mix[iside]->Fill(spdenom_mix[iside]);

        ////
        //spdenom_mix_3[iside] = sqrt( fabs(abcsp_3 * spdenom2Trk_3SE[iside] * spdenom3Trk_3SE[iside]) );
        // spdenom_mix_3[iside] = sqrt( fabs(spdenom1_3SE[iside] * spdenom2Trk_3SE[iside] * spdenom3Trk_3SE[iside]) );
        spdenom_mix_3[iside] =  spdenom2HF_3SE[iside] * spdenom3Trk_3SE[iside];
        SPdenom1_mix_3[iside]->Fill(spdenom_mix_3[iside]);
        ////

        spv1num[iside] = na1sp;
        SPv1num[iside]->Fill(na1sp);

        spv1num_mix[iside] = nacsp;
        SPv1numMix[iside]->Fill(spv1num_mix[iside]);

        spv1num_mix_3[iside] = nacsp_3;
        SPv1numMix_3[iside]->Fill(spv1num_mix_3[iside]);

        spv2num[iside] = na2sp;
        SPv2num[iside]->Fill(na2sp);

        spv3num[iside] = na3sp;
        SPv3num[iside]->Fill(na3sp);

        spv1_2SE[iside] = na1sp/spdenom1_2SE;
        SPv1_2SE[iside]->Fill(spv1_2SE[iside]);

        spv1_3SE[iside] = na1sp/spdenom1_3SE[iside];
        SPv1_3SE[iside]->Fill(spv1_3SE[iside]);

        spv1_mix[iside] = nacsp/spdenom_mix[iside];
        SPv1_mix[iside]->Fill(spv1_mix[iside]);

        spv1_mix_3[iside] = nacsp_3/spdenom_mix_3[iside];
        SPv1_mix_3[iside]->Fill(spv1_mix_3[iside]);

        spv2_2SE[iside] = na2sp/spdenom2HF_2SE;
        SPv2_2SE[iside]->Fill(spv2_2SE[iside]);

        spv2_3SE[iside] = na2sp/spdenom2HF_3SE[iside];
        SPv2_3SE[iside]->Fill(spv2_3SE[iside]);

        spv3_2SE[iside] = na3sp/spdenom3HF_2SE;
        SPv3_2SE[iside]->Fill(spv3_2SE[iside]);

        spv3_3SE[iside] = na3sp/spdenom3HF_3SE[iside];
        SPv3_3SE[iside]->Fill(spv3_3SE[iside]);

        for (int pbin = 0; pbin<nptbins; pbin++) {
            double nEP1count_pt = (double)ncount_pt[iside][pbin];
            if (eta_weights && iside == 0) nEP1count_pt*=-1;
            // This preserves the eta-weighting on the particles of interest
            // This is already accounted for in the scalar-product weights

            double na1ep_pt = Q1nAnorm_pt[iside][pbin].real()/nEP1count_pt;
            double nacep_pt = QmixnACnorm_pt[iside][pbin].real()/nEP1count_pt;
            double nacep_3_pt = QmixnACnorm_3_pt[iside][pbin].real()/nEP1count_pt;
            double na2ep_pt = Q2nAnorm_pt[iside][pbin].real()/(double)ncount_pt[iside][pbin];
            double na3ep_pt = Q3nAnorm_pt[iside][pbin].real()/(double)ncount_pt[iside][pbin];

            double na1sp_pt = Q1nA_pt[iside][pbin].real()/w1nAsum_pt[iside][pbin];
            double nacsp_pt = QmixnAC_pt[iside][pbin].real()/wmixnACsum_pt[iside][pbin];
            double nacsp_3_pt = QmixnAC_3_pt[iside][pbin].real()/wmixnACsum_3_pt[iside][pbin];
            double na2sp_pt = Q2nA_pt[iside][pbin].real()/w2nAsum_pt[iside][pbin];
            double na3sp_pt = Q3nA_pt[iside][pbin].real()/w3nAsum_pt[iside][pbin];

            epv1obs_pt[iside][pbin] = na1ep_pt;
            EPv1obs_pt[iside][pbin]->Fill(epv1obs_pt[iside][pbin]);

            epv2obs_pt[iside][pbin] = na2ep_pt;
            EPv2obs_pt[iside][pbin]->Fill(epv2obs_pt[iside][pbin]);

            epv3obs_pt[iside][pbin] = na3ep_pt;
            EPv3obs_pt[iside][pbin]->Fill(epv3obs_pt[iside][pbin]);

            epv1_2SE_pt[iside][pbin] = na1ep_pt/rescor1_2;
            EPv1_2SE_pt[iside][pbin]->Fill(epv1_2SE_pt[iside][pbin]);

            epv1_3SE_pt[iside][pbin] = na1ep_pt/rescor1_3[iside];
            EPv1_3SE_pt[iside][pbin]->Fill(epv1_3SE_pt[iside][pbin]);

            epv2_2SE_pt[iside][pbin] = na2ep_pt/rescor2HF_2;
            EPv2_2SE_pt[iside][pbin]->Fill(epv2_2SE_pt[iside][pbin]);

            epv2_3SE_pt[iside][pbin] = na2ep_pt/rescor2HF_3[iside];
            EPv2_3SE_pt[iside][pbin]->Fill(epv2_3SE_pt[iside][pbin]);

            epv3_2SE_pt[iside][pbin] = na3ep_pt/rescor3HF_2;
            EPv3_2SE_pt[iside][pbin]->Fill(epv3_2SE_pt[iside][pbin]);

            epv3_3SE_pt[iside][pbin] = na3ep_pt/rescor3HF_3[iside];
            EPv3_3SE_pt[iside][pbin]->Fill(epv3_3SE_pt[iside][pbin]);

            epv1obs_mix_pt[iside][pbin] = nacep_pt;
            EPv1obsMix_pt[iside][pbin]->Fill(epv1obs_mix_pt[iside][pbin]);

            epv1obs_mix_3_pt[iside][pbin] = nacep_3_pt;
            EPv1obsMix_3_pt[iside][pbin]->Fill(epv1obs_mix_3_pt[iside][pbin]);

            epv1_mix_pt[iside][pbin] = nacep_pt/rescorMix[iside];
            EPv1_mix_pt[iside][pbin]->Fill(epv1_mix_pt[iside][pbin]);

            epv1_mix_3_pt[iside][pbin] = nacep_3_pt/rescorMix_3[iside];
            EPv1_mix_3_pt[iside][pbin]->Fill(epv1_mix_3_pt[iside][pbin]);

            spv1num_pt[iside][pbin] = na1sp_pt;
            SPv1num_pt[iside][pbin]->Fill(na1sp_pt);

            spv1num_mix_pt[iside][pbin] = nacsp_pt;
            SPv1numMix_pt[iside][pbin]->Fill(spv1num_mix_pt[iside][pbin]);

            spv1num_mix_3_pt[iside][pbin] = nacsp_3_pt;
            SPv1numMix_3_pt[iside][pbin]->Fill(spv1num_mix_3_pt[iside][pbin]);

            spv2num_pt[iside][pbin] = na2sp_pt;
            SPv2num_pt[iside][pbin]->Fill(na2sp_pt);

            spv3num_pt[iside][pbin] = na3sp_pt;
            SPv3num_pt[iside][pbin]->Fill(na3sp_pt);

            spv1_2SE_pt[iside][pbin] = na1sp_pt/spdenom1_2SE;
            SPv1_2SE_pt[iside][pbin]->Fill(spv1_2SE_pt[iside][pbin]);

            spv1_3SE_pt[iside][pbin] = na1sp_pt/spdenom1_3SE[iside];
            SPv1_3SE_pt[iside][pbin]->Fill(spv1_3SE_pt[iside][pbin]);

            spv1_mix_pt[iside][pbin] = nacsp_pt/spdenom_mix[iside];
            SPv1_mix_pt[iside][pbin]->Fill(spv1_mix_pt[iside][pbin]);

            spv1_mix_3_pt[iside][pbin] = nacsp_3_pt/spdenom_mix_3[iside];
            if (iside == 0) spv1_mix_3_pt[iside][pbin] *= -1; // preserve reflection symmetry
            SPv1_mix_3_pt[iside][pbin]->Fill(spv1_mix_3_pt[iside][pbin]);

            spv2_2SE_pt[iside][pbin] = na2sp_pt/spdenom2HF_2SE;
            SPv2_2SE_pt[iside][pbin]->Fill(spv2_2SE_pt[iside][pbin]);

            spv2_3SE_pt[iside][pbin] = na2sp_pt/spdenom2HF_3SE[iside];
            SPv2_3SE_pt[iside][pbin]->Fill(spv2_3SE_pt[iside][pbin]);

            spv3_2SE_pt[iside][pbin] = na3sp_pt/spdenom3HF_2SE;
            SPv3_2SE_pt[iside][pbin]->Fill(spv3_2SE_pt[iside][pbin]);

            spv3_3SE_pt[iside][pbin] = na3sp_pt/spdenom3HF_3SE[iside];
            SPv3_3SE_pt[iside][pbin]->Fill(spv3_3SE_pt[iside][pbin]);

            hncnt_pt[iside][pbin]->Fill(ncount_pt[iside][pbin]);
        }
        for (int ebin = 0; ebin<netabins; ebin++) {
            double nEP1count_eta = (double)ncount_eta[iside][ebin];
            if (eta_weights && iside == 0) nEP1count_eta*=-1;
            // This preserves the eta-weighting on the particles of interest
            // This is already accounted for in the scalar-product weights

            double na1ep_eta = Q1nAnorm_eta[iside][ebin].real()/nEP1count_eta;
            double nacep_eta = QmixnACnorm_eta[iside][ebin].real()/nEP1count_eta;
            double nacep_3_eta = QmixnACnorm_3_eta[iside][ebin].real()/nEP1count_eta;

            double na1sp_eta = Q1nA_eta[iside][ebin].real()/w1nAsum_eta[iside][ebin];
            double nacsp_eta = QmixnAC_eta[iside][ebin].real()/wmixnACsum_eta[iside][ebin];
            double nacsp_3_eta = QmixnAC_3_eta[iside][ebin].real()/wmixnACsum_3_eta[iside][ebin];

            epv1obs_eta[iside][ebin] = na1ep_eta;
            EPv1obs_eta[iside][ebin]->Fill(epv1obs_eta[iside][ebin]);

            epv1_2SE_eta[iside][ebin] = na1ep_eta/rescor1_2;
            EPv1_2SE_eta[iside][ebin]->Fill(epv1_2SE_eta[iside][ebin]);

            epv1_3SE_eta[iside][ebin] = na1ep_eta/rescor1_3[iside];
            EPv1_3SE_eta[iside][ebin]->Fill(epv1_3SE_eta[iside][ebin]);

            epv1obs_mix_eta[iside][ebin] = nacep_eta;
            EPv1obsMix_eta[iside][ebin]->Fill(epv1obs_mix_eta[iside][ebin]);

            epv1obs_mix_3_eta[iside][ebin] = nacep_eta;
            EPv1obsMix_3_eta[iside][ebin]->Fill(epv1obs_mix_3_eta[iside][ebin]);

            epv1_mix_eta[iside][ebin] = nacep_eta/rescorMix[iside];
            EPv1_mix_eta[iside][ebin]->Fill(epv1_mix_eta[iside][ebin]);

            epv1_mix_3_eta[iside][ebin] = nacep_3_eta/rescorMix_3[iside];
            EPv1_mix_3_eta[iside][ebin]->Fill(epv1_mix_3_eta[iside][ebin]);

            spv1num_eta[iside][ebin] = na1sp_eta;
            SPv1num_eta[iside][ebin]->Fill(na1sp_eta);

            spv1num_mix_eta[iside][ebin] = nacsp_eta;
            SPv1numMix_eta[iside][ebin]->Fill(spv1num_mix_eta[iside][ebin]);

            spv1num_mix_3_eta[iside][ebin] = nacsp_3_eta;
            SPv1numMix_3_eta[iside][ebin]->Fill(spv1num_mix_3_eta[iside][ebin]);

            spv1_2SE_eta[iside][ebin] = na1sp_eta/spdenom1_2SE;
            SPv1_2SE_eta[iside][ebin]->Fill(spv1_2SE_eta[iside][ebin]);

            spv1_3SE_eta[iside][ebin] = na1sp_eta/spdenom1_3SE[iside];
            SPv1_3SE_eta[iside][ebin]->Fill(spv1_3SE_eta[iside][ebin]);

            spv1_mix_eta[iside][ebin] = nacsp_eta/spdenom_mix[iside];
            SPv1_mix_eta[iside][ebin]->Fill(spv1_mix_eta[iside][ebin]);

            spv1_mix_3_eta[iside][ebin] = nacsp_3_eta/spdenom_mix_3[iside];
            if (iside == 0) spv1_mix_3_eta[iside][ebin] *= -1; // preserve reflection symmetry
            SPv1_mix_3_eta[iside][ebin]->Fill(spv1_mix_3_eta[iside][ebin]);

            hncnt_eta[iside][ebin]->Fill(ncount_eta[iside][ebin]);

            // Q-vector averages
            hQ1nA[iside][ebin]->Fill(na1sp_eta);
            hQmixnAC[iside][ebin]->Fill(nacsp_eta);
            hQmixnAC_3[iside][ebin]->Fill(nacsp_3_eta);

            // normalized Q-vector averages
            hQ1nAnorm[iside][ebin]->Fill(na1ep_eta);
            hQmixnACnorm[iside][ebin]->Fill(nacep_eta);
            hQmixnAC_3norm[iside][ebin]->Fill(nacep_3_eta);
        }
    }

    ofstream fout;
    if (!fopen("logs","r")) system(Form("mkdir logs"));
    TString tag = Form("logs_%s",mtag.Data());
    TString foutname = "logs/"+tag+".dat";
    fout.open(foutname.Data());
    fout << "====================" << endl;
    fout << "nevents:    " << nevents << endl;
    if (v1odd) fout << "Rapidity-odd v1 input " << endl;
    else fout << "Rapidity-even v1:  " << setv1 << endl;
    fout << "Input v2:     " << setv2 << endl;
    fout << "Input v3:     " << setv3 << endl;
    if (eta_weights) fout << "Using eta-dependent weights..." <<endl;
    if (pt_weights) fout << "Using pt-dependent weights..." <<endl;
    if (conserve_pT) fout << "Conserving transverse momentum event by event... " << endl;
    if (addHoles) fout << "Holes in detector acceptance... " << endl;
    fout << "iseed:    " << iseed << endl;
    for (int nep = 0; nep<numEP-1; nep++) {
        fout << "offsets: "
        << EPName[nep].data() << "\t" << X1[nep] << "\t" << Y1[nep] << "\t" << X2[nep] << "\t" << Y2[nep] << endl;
    }

    fout<<"   -----             "<<endl;
    fout<<"nevts:               "<<nevts<<endl;
    fout<<"evtmult:             "<<evtmult<<endl;
    fout<<"ncount (HFm):        "<<ncount[0]<<endl;
    fout<<"ncount (HFp):        "<<ncount[1]<<endl;
    fout<<"   -----             "<<endl;
    fout<<"HFm mult:            "<<subcnt[0]<<endl;
    fout<<"trackm mult:         "<<subcnt[1]<<endl;
    fout<<"trackmid mult:       "<<subcnt[4]<<endl;
    fout<<"trackp mult:         "<<subcnt[2]<<endl;
    fout<<"HFp mult:            "<<subcnt[3]<<endl;
    fout<<"particles rejected:  "<<subcnt[5]<<endl;
    fout<<"   -----             "<<endl;
    fout<<"rescor1_2:           "<<Form("%.6f",rescor1_2)       <<endl;
    fout<<"rescor1_3     (HFm): "<<Form("%.6f",rescor1_3[0])    <<"\t (HFp): "<<Form("%.6f",rescor1_3[1])<<endl;
    fout<<"rescor2HF_2:           "<<Form("%.6f",rescor2HF_2)       <<endl;
    fout<<"rescor2HF_3     (HFm): "<<Form("%.6f",rescor2HF_3[0])    <<"\t (HFp): "<<Form("%.6f",rescor2HF_3[1])<<endl;
    fout<<"rescor3HF_2:           "<<Form("%.6f",rescor3HF_2)       <<endl;
    fout<<"rescor3HF_3     (HFm): "<<Form("%.6f",rescor3HF_3[0])    <<"\t (HFp): "<<Form("%.6f",rescor3HF_3[1])<<endl;
    fout<<"rescorMix     (HFm): "<<Form("%.6f",rescorMix[0])    <<"\t (HFp): "<<Form("%.6f",rescorMix[1])<<endl;
    fout<<"rescorMix_3   (HFm): "<<Form("%.6f",rescorMix_3[0])  <<"\t (HFp): "<<Form("%.6f",rescorMix_3[1])<<endl;
    fout<<"epv1obs       (HFm): "<<Form("%.6f",epv1obs[0])      <<"\t (HFp): "<<Form("%.6f",epv1obs[1])<<endl;
    fout<<"epv1obs_mix   (HFm): "<<Form("%.6f",epv1obs_mix[0])  <<"\t (HFp): "<<Form("%.6f",epv1obs_mix[1])<<endl;
    fout<<"epv1obs_mix_3 (HFm): "<<Form("%.6f",epv1obs_mix_3[0])<<"\t (HFp): "<<Form("%.6f",epv1obs_mix_3[1])<<endl;
    fout<<"epv2obs       (HFm): "<<Form("%.6f",epv2obs[0])      <<"\t (HFp): "<<Form("%.6f",epv2obs[1])<<endl;
    fout<<"epv3obs       (HFm): "<<Form("%.6f",epv3obs[0])      <<"\t (HFp): "<<Form("%.6f",epv3obs[1])<<endl;
    fout<<"   -----             "<<endl;
    fout<<"spdenom1_2SE:        "<<Form("%.6f",spdenom1_2SE)    <<endl;
    fout<<"spdenom1_3SE  (HFm): "<<Form("%.6f",spdenom1_3SE[0])<<"\t (HFp): "<<Form("%.6f",spdenom1_3SE[1])<<endl;
    fout<<"SPdenom2HF_2SE:        "<<Form("%.6f",spdenom2HF_2SE)    <<endl;
    fout<<"SPdenom2HF_3SE  (HFm): "<<Form("%.6f",spdenom2HF_3SE[0])<<"\t (HFp): "<<Form("%.6f",spdenom2HF_3SE[1])<<endl;
    fout<<"SPdenom3HF_2SE:        "<<Form("%.6f",spdenom3HF_2SE)    <<endl;
    fout<<"SPdenom3HF_3SE  (HFm): "<<Form("%.6f",spdenom3HF_3SE[0])<<"\t (HFp): "<<Form("%.6f",spdenom3HF_3SE[1])<<endl;
    fout<<"spdenom_mix   (HFm): "<<Form("%.6f",spdenom_mix[0])  <<"\t (HFp): "<<Form("%.6f",spdenom_mix[1])<<endl;
    fout<<"spdenom_mix_3 (HFm): "<<Form("%.6f",spdenom_mix_3[0])<<"\t (HFp): "<<Form("%.6f",spdenom_mix_3[1])<<endl;
    fout<<"spv1num       (HFm): "<<Form("%.6f",spv1num[0])      <<"\t (HFp): "<<Form("%.6f",spv1num[1])<<endl;
    fout<<"spv2num       (HFm): "<<Form("%.6f",spv2num[0])      <<"\t (HFp): "<<Form("%.6f",spv2num[1])<<endl;
    fout<<"spv3num       (HFm): "<<Form("%.6f",spv3num[0])      <<"\t (HFp): "<<Form("%.6f",spv3num[1])<<endl;
    fout<<"spv1num_mix   (HFm): "<<Form("%.6f",spv1num_mix[0])  <<"\t (HFp): "<<Form("%.6f",spv1num_mix[1])<<endl;
    fout<<"spv1num_mix_3 (HFm): "<<Form("%.6f",spv1num_mix_3[0])<<"\t (HFp): "<<Form("%.6f",spv1num_mix_3[1])<<endl;
    fout<<"   -----             "<<endl;
    fout<<"epv1_2SE      (HFm): "<<Form("%.6f",epv1_2SE[0])     <<"\t (HFp): "<<Form("%.6f",epv1_2SE[1])<<endl;
    fout<<"epv1_3SE      (HFm): "<<Form("%.6f",epv1_3SE[0])     <<"\t (HFp): "<<Form("%.6f",epv1_3SE[1])<<endl;
    fout<<"epv1_mix      (HFm): "<<Form("%.6f",epv1_mix[0])     <<"\t (HFp): "<<Form("%.6f",epv1_mix[1])<<endl;
    fout<<"epv1_mix_3    (HFm): "<<Form("%.6f",epv1_mix_3[0])   <<"\t (HFp): "<<Form("%.6f",epv1_mix_3[1])<<endl;
    fout<<"   -----             "<<endl;
    fout<<"spv1_2SE      (HFm): "<<Form("%.6f",spv1_2SE[0])    <<"\t (HFp): "<<Form("%.6f",spv1_2SE[1])<<endl;
    fout<<"spv1_3SE      (HFm): "<<Form("%.6f",spv1_3SE[0])    <<"\t (HFp): "<<Form("%.6f",spv1_3SE[1])<<endl;
    fout<<"spv1_mix      (HFm): "<<Form("%.6f",spv1_mix[0])    <<"\t (HFp): "<<Form("%.6f",spv1_mix[1])<<endl;
    fout<<"spv1_mix_3    (HFm): "<<Form("%.6f",spv1_mix_3[0])  <<"\t (HFp): "<<Form("%.6f",spv1_mix_3[1])<<endl;
    fout<<"   -----             "<<endl;
    fout<<"epv2_2SE      (HFm): "<<Form("%.6f",epv2_2SE[0])    <<"\t (HFp): "<<Form("%.6f",epv2_2SE[1])<<endl;
    fout<<"epv2_3SE      (HFm): "<<Form("%.6f",epv2_3SE[0])    <<"\t (HFp): "<<Form("%.6f",epv2_3SE[1])<<endl;
    fout<<"spv2_2SE      (HFm): "<<Form("%.6f",spv2_2SE[0])    <<"\t (HFp): "<<Form("%.6f",spv2_2SE[1])<<endl;
    fout<<"spv2_3SE      (HFm): "<<Form("%.6f",spv2_3SE[0])    <<"\t (HFp): "<<Form("%.6f",spv2_3SE[1])<<endl;
    fout<<"   -----             "<<endl;
    fout<<"epv3_2SE      (HFm): "<<Form("%.6f",epv3_2SE[0])    <<"\t (HFp): "<<Form("%.6f",epv3_2SE[1])<<endl;
    fout<<"epv3_3SE      (HFm): "<<Form("%.6f",epv3_3SE[0])    <<"\t (HFp): "<<Form("%.6f",epv3_3SE[1])<<endl;
    fout<<"spv3_2SE      (HFm): "<<Form("%.6f",spv3_2SE[0])    <<"\t (HFp): "<<Form("%.6f",spv3_2SE[1])<<endl;
    fout<<"spv3_3SE      (HFm): "<<Form("%.6f",spv3_3SE[0])    <<"\t (HFp): "<<Form("%.6f",spv3_3SE[1])<<endl;
    fout<<"   -----             \n"<<endl;
    fout.close();

    cout<<"   -----             "<<endl;
    cout<<"nevts:               "<<nevts<<endl;
    cout<<"evtmult:             "<<evtmult<<endl;
    cout<<"ncount (HFm):        "<<ncount[0]<<endl;
    cout<<"ncount (HFp):        "<<ncount[1]<<endl;
    cout<<"   -----             "<<endl;
    cout<<"HFm mult:            "<<subcnt[0]<<endl;
    cout<<"trackm mult:         "<<subcnt[1]<<endl;
    cout<<"trackmid mult:       "<<subcnt[4]<<endl;
    cout<<"trackp mult:         "<<subcnt[2]<<endl;
    cout<<"HFp mult:            "<<subcnt[3]<<endl;
    cout<<"particles rejected:  "<<subcnt[5]<<endl;
    cout<<"   -----             "<<endl;
    cout<<"rescor1_2:           "<<Form("%.6f",rescor1_2)       <<endl;
    cout<<"rescor1_3     (HFm): "<<Form("%.6f",rescor1_3[0])    <<"\t (HFp): "<<Form("%.6f",rescor1_3[1])<<endl;
    cout<<"rescor2HF_2:           "<<Form("%.6f",rescor2HF_2)       <<endl;
    cout<<"rescor2HF_3     (HFm): "<<Form("%.6f",rescor2HF_3[0])    <<"\t (HFp): "<<Form("%.6f",rescor2HF_3[1])<<endl;
    cout<<"rescor3HF_2:           "<<Form("%.6f",rescor3HF_2)       <<endl;
    cout<<"rescor3HF_3     (HFm): "<<Form("%.6f",rescor3HF_3[0])    <<"\t (HFp): "<<Form("%.6f",rescor3HF_3[1])<<endl;
    cout<<"rescorMix     (HFm): "<<Form("%.6f",rescorMix[0])    <<"\t (HFp): "<<Form("%.6f",rescorMix[1])<<endl;
    cout<<"rescorMix_3   (HFm): "<<Form("%.6f",rescorMix_3[0])  <<"\t (HFp): "<<Form("%.6f",rescorMix_3[1])<<endl;
    cout<<"epv1obs       (HFm): "<<Form("%.6f",epv1obs[0])      <<"\t (HFp): "<<Form("%.6f",epv1obs[1])<<endl;
    cout<<"epv1obs_mix   (HFm): "<<Form("%.6f",epv1obs_mix[0])  <<"\t (HFp): "<<Form("%.6f",epv1obs_mix[1])<<endl;
    cout<<"epv1obs_mix_3 (HFm): "<<Form("%.6f",epv1obs_mix_3[0])<<"\t (HFp): "<<Form("%.6f",epv1obs_mix_3[1])<<endl;
    cout<<"epv2obs       (HFm): "<<Form("%.6f",epv2obs[0])      <<"\t (HFp): "<<Form("%.6f",epv2obs[1])<<endl;
    cout<<"epv3obs       (HFm): "<<Form("%.6f",epv3obs[0])      <<"\t (HFp): "<<Form("%.6f",epv3obs[1])<<endl;
    cout<<"   -----             "<<endl;
    cout<<"spdenom1_2SE:        "<<Form("%.6f",spdenom1_2SE)    <<endl;
    cout<<"spdenom1_3SE  (HFm): "<<Form("%.6f",spdenom1_3SE[0])<<"\t (HFp): "<<Form("%.6f",spdenom1_3SE[1])<<endl;
    cout<<"SPdenom2HF_2SE:        "<<Form("%.6f",spdenom2HF_2SE)    <<endl;
    cout<<"SPdenom2HF_3SE  (HFm): "<<Form("%.6f",spdenom2HF_3SE[0])<<"\t (HFp): "<<Form("%.6f",spdenom2HF_3SE[1])<<endl;
    cout<<"SPdenom3HF_2SE:        "<<Form("%.6f",spdenom3HF_2SE)    <<endl;
    cout<<"SPdenom3HF_3SE  (HFm): "<<Form("%.6f",spdenom3HF_3SE[0])<<"\t (HFp): "<<Form("%.6f",spdenom3HF_3SE[1])<<endl;
    cout<<"spdenom_mix   (HFm): "<<Form("%.6f",spdenom_mix[0])  <<"\t (HFp): "<<Form("%.6f",spdenom_mix[1])<<endl;
    cout<<"spdenom_mix_3 (HFm): "<<Form("%.6f",spdenom_mix_3[0])<<"\t (HFp): "<<Form("%.6f",spdenom_mix_3[1])<<endl;
    cout<<"spv1num       (HFm): "<<Form("%.6f",spv1num[0])      <<"\t (HFp): "<<Form("%.6f",spv1num[1])<<endl;
    cout<<"spv2num       (HFm): "<<Form("%.6f",spv2num[0])      <<"\t (HFp): "<<Form("%.6f",spv2num[1])<<endl;
    cout<<"spv3num       (HFm): "<<Form("%.6f",spv3num[0])      <<"\t (HFp): "<<Form("%.6f",spv3num[1])<<endl;
    cout<<"spv1num_mix   (HFm): "<<Form("%.6f",spv1num_mix[0])  <<"\t (HFp): "<<Form("%.6f",spv1num_mix[1])<<endl;
    cout<<"spv1num_mix_3 (HFm): "<<Form("%.6f",spv1num_mix_3[0])<<"\t (HFp): "<<Form("%.6f",spv1num_mix_3[1])<<endl;
    cout<<"   -----             "<<endl;
    cout<<"epv1_2SE      (HFm): "<<Form("%.6f",epv1_2SE[0])     <<"\t (HFp): "<<Form("%.6f",epv1_2SE[1])<<endl;
    cout<<"epv1_3SE      (HFm): "<<Form("%.6f",epv1_3SE[0])     <<"\t (HFp): "<<Form("%.6f",epv1_3SE[1])<<endl;
    cout<<"epv1_mix      (HFm): "<<Form("%.6f",epv1_mix[0])     <<"\t (HFp): "<<Form("%.6f",epv1_mix[1])<<endl;
    cout<<"epv1_mix_3    (HFm): "<<Form("%.6f",epv1_mix_3[0])   <<"\t (HFp): "<<Form("%.6f",epv1_mix_3[1])<<endl;
    cout<<"   -----             "<<endl;
    cout<<"spv1_2SE      (HFm): "<<Form("%.6f",spv1_2SE[0])    <<"\t (HFp): "<<Form("%.6f",spv1_2SE[1])<<endl;
    cout<<"spv1_3SE      (HFm): "<<Form("%.6f",spv1_3SE[0])    <<"\t (HFp): "<<Form("%.6f",spv1_3SE[1])<<endl;
    cout<<"spv1_mix      (HFm): "<<Form("%.6f",spv1_mix[0])    <<"\t (HFp): "<<Form("%.6f",spv1_mix[1])<<endl;
    cout<<"spv1_mix_3    (HFm): "<<Form("%.6f",spv1_mix_3[0])  <<"\t (HFp): "<<Form("%.6f",spv1_mix_3[1])<<endl;
    cout<<"   -----             "<<endl;
    cout<<"epv2_2SE      (HFm): "<<Form("%.6f",epv2_2SE[0])    <<"\t (HFp): "<<Form("%.6f",epv2_2SE[1])<<endl;
    cout<<"epv2_3SE      (HFm): "<<Form("%.6f",epv2_3SE[0])    <<"\t (HFp): "<<Form("%.6f",epv2_3SE[1])<<endl;
    cout<<"spv2_2SE      (HFm): "<<Form("%.6f",spv2_2SE[0])    <<"\t (HFp): "<<Form("%.6f",spv2_2SE[1])<<endl;
    cout<<"spv2_3SE      (HFm): "<<Form("%.6f",spv2_3SE[0])    <<"\t (HFp): "<<Form("%.6f",spv2_3SE[1])<<endl;
    cout<<"   -----             "<<endl;
    cout<<"epv3_2SE      (HFm): "<<Form("%.6f",epv3_2SE[0])    <<"\t (HFp): "<<Form("%.6f",epv3_2SE[1])<<endl;
    cout<<"epv3_3SE      (HFm): "<<Form("%.6f",epv3_3SE[0])    <<"\t (HFp): "<<Form("%.6f",epv3_3SE[1])<<endl;
    cout<<"spv3_2SE      (HFm): "<<Form("%.6f",spv3_2SE[0])    <<"\t (HFp): "<<Form("%.6f",spv3_2SE[1])<<endl;
    cout<<"spv3_3SE      (HFm): "<<Form("%.6f",spv3_3SE[0])    <<"\t (HFp): "<<Form("%.6f",spv3_3SE[1])<<endl;
    cout<<"   -----             \n"<<endl;

}
