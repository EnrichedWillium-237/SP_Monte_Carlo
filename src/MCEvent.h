#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

class MCEvent
{
public:
    MCEvent(Double_t tv1,Double_t tv2,Double_t tv3,Double_t tv4,Double_t tv5,Double_t tv6);
    void SetEventParms();
    void SetMult(Int_t value);
    void SetPsiRandom();
    void SetPsi(Double_t value){Psi = value;}
    void GetStaticEta(Double_t value, Double_t etaTrackArray[]);
    void GetEtaRandom(Double_t etaTrackArray[], Double_t etamin, Double_t etamax);
    void GetStaticPt(Double_t value, Double_t ptTrackArray[]);
    void GetPtRandom(Double_t ptTrackArray[]);
    void Setv1(Double_t value);
    void Setv2(Double_t value);
    void Setv3(Double_t value);
    void Setv4(Double_t value);
    void Setv5(Double_t value);
    void Setv6(Double_t value);
    Int_t GetMult();
    void GetThrowPhi(Double_t phiTrackArray[]);
    void GetThrowPhiHole(Double_t phiTrackArray[], Double_t minh, Double_t maxh);
    Double_t GetPsi();
    Bool_t IsPhiHole();
    void SetSeed(Int_t seed);
    void SetOffset(Double_t val){valueOff = val;}

private:
    Double_t v1;
    Double_t v2;
    Double_t v3;
    Double_t v4;
    Double_t v5;
    Double_t v6;
    Double_t Psi;
    Double_t valueRad;
    Double_t valueOff;
    Int_t mult;
    Double_t cent;
    TRandom * ran;
    TF1 * phidist;
    TF1 * ptdist;
    TF1 * dNdetadist;
    Bool_t setphiHole;
};


MCEvent::MCEvent(Double_t tv1,Double_t tv2,Double_t tv3,Double_t tv4,Double_t tv5,Double_t tv6) {
    setphiHole = kFALSE;
    valueRad = 1.0;
    valueOff = 0.0;
    ran = new TRandom3(0);
    v1 = tv1;
    v2 = tv2;
    v3 = tv3;
    v4 = tv4;
    v5 = tv5;
    v6 = tv6;
    phidist = new TF1("phidist","1+2*[0]*cos(x)+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x)+2*[5]*cos(6*x)",-TMath::Pi(),TMath::Pi());
//    ptdist = new TF1("ptdist","x*((([0]-1)*([0]-2))/([1]*[1]))*[2]*TMath::Power(1+x/[1],-[0])",0.0,5.0);
    ptdist = new TF1("ptdist","x*((([0]-1)*([0]-2))/([1]*[1]))*[2]*TMath::Power(1+x/[1],-[0])",0.0,8.0);
    ptdist->SetParameters(6.06271, 1.08372, 586.474);

    //    etadist = new TF1("etadist","( [0] + [1]*x - [2]*x*x - [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x )/[7]",-5,5);
    //    etadist->SetParameters(1.01669e+06, 6124.96, 169458., 1323.27, 1191.52, 71.6021, 565.621, 1.0e+06);
    dNdetadist = new TF1("dNdetadist","[0]*TMath::Exp((-1*x*x)/([1]*[1])) - [2]*TMath::Exp((-1*x*x)/([3]*[3]))",-5.6,5.6);
    //    dNdetadist->SetParameters(2088, 5.23, 470, 1.63); // cent 0-5%
    //    dNdetadist->SetParameters(1704, 5.30, 387, 1.65); // cent 5-10%
    //    dNdetadist->SetParameters(1896, 5.26, 429, 1.64);  // cent 0-10%
    //    dNdetadist->SetParameters(1263, 5.39, 281, 1.60); // cent 10-20%
    dNdetadist->SetParameters(873, 5.42, 205, 1.69);  // cent 20-30%
    //    dNdetadist->SetParameters(567, 5.54, 143, 1.76);  // cent 30-40%
    //    dNdetadist->SetParameters(357, 5.51, 97, 1.88);  // cent 40-50%
}

void MCEvent::SetMult(Int_t value) {mult=value;}
void MCEvent::SetPsiRandom() {
    Psi = ran->Uniform(-TMath::Pi(),TMath::Pi());
}

void MCEvent::SetEventParms() {
    phidist->SetParameters(v1,v2,v3,v4,v5,v6);
}

Int_t MCEvent::GetMult() {return mult;}

void MCEvent::GetStaticEta(Double_t value, Double_t etaTrackArray[]) {
    for(Int_t i = 0; i<mult; i++) {
        etaTrackArray[i] = value;
    }
}

void MCEvent::GetEtaRandom(Double_t etaTrackArray[], Double_t emin, Double_t emax) {
    for(Int_t i = 0; i<mult; i++) {
        etaTrackArray[i] = dNdetadist->GetRandom(emin,emax);
    }
}

void MCEvent::GetStaticPt(Double_t value, Double_t ptTrackArray[]) {
    for(Int_t i = 0; i<mult; i++) {
        ptTrackArray[i] = value;
    }
}

void MCEvent::GetPtRandom(Double_t ptTrackArray[]) {
    for (Int_t i = 0; i<mult; i++) {
        ptTrackArray[i] = ran->Uniform(0.1,7.8);
//        ptTrackArray[i] = ptdist->GetRandom();
    }
}

Double_t MCEvent::GetPsi() {
    return Psi;
}

void MCEvent::GetThrowPhi(Double_t phiTrackArray[]) {
    setphiHole = kFALSE;
    for (Int_t i = 0; i<mult; i++) {
        Double_t phi = phidist->GetRandom()+Psi;
//        Double_t phi = phidist->GetRandom();
        Double_t phiOffset = phi;
        Double_t R = valueRad;
        Double_t d = valueOff;
        if (valueOff!=0) {
            if (phi>=0) {
                phiOffset = TMath::ACos((d - R*TMath::Cos(TMath::Pi()-phi))/TMath::Power((d*d + R*R - 2*d*R*TMath::Cos(TMath::Pi()-phi)),.5));
            } else {
                phiOffset = -TMath::ACos((d - R*TMath::Cos(TMath::Pi()-phi))/TMath::Power((d*d + R*R - 2*d*R*TMath::Cos(TMath::Pi()-phi)),.5));
            }
        }
        phiTrackArray[i] = phiOffset;
    }
}

void MCEvent::GetThrowPhiHole(Double_t phiTrackArray[], Double_t minh, Double_t maxh) {
    setphiHole = kTRUE;
    for(Int_t i = 0; i<mult; i++) {
        Double_t ph = phidist->GetRandom()+Psi;
        while(ph>minh && ph<maxh) {
            ph=phidist->GetRandom()+Psi;
        }
        phiTrackArray[i] = ph;
    }
}

void MCEvent::Setv1(Double_t value) { v1 = value; }
void MCEvent::Setv2(Double_t value) { v2 = value; }
void MCEvent::Setv3(Double_t value) { v3 = value; }
void MCEvent::Setv4(Double_t value) { v4 = value; }
void MCEvent::Setv5(Double_t value) { v5 = value; }
void MCEvent::Setv6(Double_t value) { v6 = value; }
Bool_t MCEvent::IsPhiHole() {return setphiHole;}

void MCEvent::SetSeed(Int_t seed) {
    ran->SetSeed(seed);
    gRandom->SetSeed(seed);
}
