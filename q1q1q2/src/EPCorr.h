void EPCorr(Double_t HFm1, Double_t trackm1, Double_t trackp1, Double_t HFp1, Double_t trackmid1,
            Double_t HFm2, Double_t trackm2, Double_t trackp2, Double_t HFp2, Double_t trackmid2)
{
    corr_HFm1_HFp1->Fill(HFm1, HFp1);
    corr_trackmid1_HFm1->Fill(trackmid1, HFm1);
    corr_trackmid1_HFp1->Fill(trackmid1, HFp1);
    corr_trackm1_HFm1->Fill(trackm1, HFm1);
    corr_trackm1_trackp1->Fill(trackm1, trackp1);
    corr_trackm1_HFp1->Fill(trackm1, HFp1);
    corr_trackp1_HFm1->Fill(trackp1, HFm1);
    corr_trackp1_trackm1->Fill(trackp1, trackm1);
    corr_trackp1_HFp1->Fill(trackp1, HFp1);
    
    corr_HFm1_trackm2->Fill(HFm1, trackm2);
    corr_HFm1_trackp2->Fill(HFm1, trackp2);
    corr_HFp1_trackm2->Fill(HFp1, trackm2);
    corr_HFp1_trackp2->Fill(HFp1, trackp2);
    corr_HFm1_HFp2->Fill(HFm1, HFp2);
    corr_HFp1_HFm2->Fill(HFp1, HFm2);
    corr_trackmid1_HFm2->Fill(trackmid1, HFm2);
    corr_trackmid1_HFp2->Fill(trackmid1, HFp2);
    corr_trackm1_HFm2->Fill(trackm1, HFm2);
    corr_trackm1_trackp2->Fill(trackm1, trackp2);
    corr_trackm1_HFp2->Fill(trackm1, HFp2);
    corr_trackp1_HFm2->Fill(trackp1, HFm2);
    corr_trackp1_trackm2->Fill(trackp1, trackm2);
    corr_trackp1_HFp2->Fill(trackp1, HFp2);
    
    corr_HFm2_HFp2->Fill(HFm2, HFp2);
    corr_trackmid2_HFm2->Fill(trackmid2, HFm2);
    corr_trackmid2_HFp2->Fill(trackmid2, HFp2);
    corr_trackm2_HFm2->Fill(trackm2, HFm2);
    corr_trackm2_trackp2->Fill(trackm2, trackp2);
    corr_trackm2_HFp2->Fill(trackm2, HFp2);
    corr_trackp2_HFm2->Fill(trackp2, HFm2);
    corr_trackp2_HFp2->Fill(trackp2, HFp2);
    
    //-- cosines of correlations
    
    cos_HFm1_HFp1->Fill(TMath::Cos(HFm1 - HFp1));
    cos_trackmid1_HFm1->Fill(TMath::Cos(trackmid1 - HFm1));
    cos_trackmid1_HFp1->Fill(TMath::Cos(trackmid1 - HFp1));
    cos_trackm1_HFm1->Fill(TMath::Cos(trackm1 - HFm1));
    cos_trackm1_trackp1->Fill(TMath::Cos(trackm1 - trackp1));
    cos_trackm1_HFp1->Fill(TMath::Cos(trackm1 - HFp1));
    cos_trackp1_HFm1->Fill(TMath::Cos(trackp1 - HFm1));
    cos_trackp1_trackm1->Fill(TMath::Cos(trackm1 - trackp1));
    cos_trackp1_HFp1->Fill(TMath::Cos(trackp1 - HFp1));
    
    cos_HFm1_trackm2->Fill(TMath::Cos(2*(HFm1 - trackm2)));
    cos_HFm1_trackp2->Fill(TMath::Cos(2*(HFm1 - trackp2)));
    cos_HFp1_trackm2->Fill(TMath::Cos(2*(HFp1 - trackm2)));
    cos_HFp1_trackp2->Fill(TMath::Cos(2*(HFp1 - trackp2)));
    cos_HFm1_HFp2->Fill(TMath::Cos(2*(HFm1 - HFp2)));
    cos_HFp1_HFm2->Fill(TMath::Cos(2*(HFp1 - HFm2)));
    cos_trackmid1_HFm2->Fill(TMath::Cos(2*(trackmid1 - HFm2)));
    cos_trackmid1_HFp2->Fill(TMath::Cos(2*(trackmid1 - HFp2)));
    cos_trackm1_HFm2->Fill(TMath::Cos(2*(trackm1 - HFm2)));
    cos_trackm1_trackp2->Fill(TMath::Cos(2*(trackm1 - trackp2)));
    cos_trackm1_HFp2->Fill(TMath::Cos(2*(trackm1 - HFp2)));
    cos_trackp1_HFm2->Fill(TMath::Cos(2*(trackp1 - HFm2)));
    cos_trackp1_trackm2->Fill(TMath::Cos(2*(trackp1 - trackm2)));
    cos_trackp1_HFp2->Fill(TMath::Cos(2*(trackp1 - HFp2)));
    
    cos_HFm2_HFp2->Fill(TMath::Cos(HFm2 - HFp2));
    cos_trackmid2_HFm2->Fill(TMath::Cos(trackmid2 - HFm2));
    cos_trackmid2_HFp2->Fill(TMath::Cos(trackmid2 - HFp2));
    cos_trackm2_HFm2->Fill(TMath::Cos(trackm2 - HFm2));
    cos_trackm2_trackp2->Fill(TMath::Cos(trackm2 - trackp2));
    cos_trackm2_HFp2->Fill(TMath::Cos(trackm2 - HFp2));
    cos_trackp2_HFm2->Fill(TMath::Cos(trackp2 - HFm2));
    cos_trackp2_HFp2->Fill(TMath::Cos(trackp2 - HFp2));

}

