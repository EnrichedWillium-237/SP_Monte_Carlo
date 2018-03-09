# include "TCanvas.h"
# include "TF1.h"
# include <iostream>

using namespace std;

static const double v1Ebins[] = {
        -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0,
        -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,
         2.4,  2.8,  3.2,  3.6,  4.0,  4.4,  4.8,  5.2,  5.6};

void ALICEcent() {
    TF1 * f0to5 = new TF1("f0to5","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f0to5->SetParameters(1671.1, -8.825, 12.275, 1.9336, -4.2109, -0.0778, 0.0955);
    f0to5->Draw();
    cout<<"0 to 5%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f0to5->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;

    TF1 * f5to10 = new TF1("f5to10","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f5to10->SetParameters(1372.6, -5.1456, 4.0703, -1.9022, -3.6521, 0.0643, 0.095);
    f5to10->Draw("same");
    cout<<"5 to 10%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f5to10->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;

    TF1 * f10to20 = new TF1("f10to20","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f10to20->SetParameters(1013, -4.6706, 10.84, 1.0302, -2.7517, -0.0418, 0.0621);
    f10to20->Draw("same");
    cout<<"10 to 20%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f10to20->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;

    TF1 * f30to40 = new TF1("f30to40","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f30to40->SetParameters(431.84, -0.6658, 12.898, 0.11, -2.0092, -0.0012, 0.0494);
    f30to40->Draw("same");
    cout<<"30 to 40%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f30to40->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;

    TF1 * f40to50 = new TF1("f40to50","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f40to50->SetParameters(263.69, - 0.3419, 8.4487, 0.0707, -1.2448, -0.0027, 0.0306);
    f40to50->Draw("same");
    cout<<"40 to 50%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f40to50->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;

    TF1 * f50to60 = new TF1("f50to60","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f50to60->SetParameters(149.31, -0.1938, 5.4365, 0.0445, -0.7445, -0.0022, 0.0182);
    f50to60->Draw("same");
    cout<<"50 to 60%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f50to60->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;

    TF1 * f60to70 = new TF1("f60to70","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",-5.7,5.7);
    f60to70->SetParameters(76.21, -0.1249, 3.1879, 0.0458, -0.4012, -0.0037, 0.01);
    f60to70->Draw("same");
    cout<<"60 to 70%"<<endl;
    for (int i = 0; i<28; i++) {
        cout<<(int)f60to70->Integral(v1Ebins[i], v1Ebins[i+1])<<", ";
    }
    cout<<"\n"<<endl;


}
