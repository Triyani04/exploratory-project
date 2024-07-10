#include<bits/stdc++.h>
using namespace std;
int main()
{
    // Given parameters
    double phi1,phi2;
    double K1,K2;
    double T=100;   // in Degree celcius
    double x1=0.958;
    double x2=1-x1;
    double R=8.314;
    double  P1_sat=exp(16.5938-(3644.3/(T+239.76)))*1000;
    double  P2_sat=exp(16.2620-(3799.89/(T+226.35)))*1000;
    // cout << "P1_sat " << P1_sat << endl;
    // cout << "P2_sat " << P2_sat << endl;
    double a12=107.38*4.186;
    double a21=469.55*4.186;
    double V1=40.73;
    double V2=18.07;
    // Values Given
    double Tc1=512.6;
    double Tc2=647.1;
    double Pc1=80.97*pow(10,5);
    double Pc2=220.55*pow(10,5);
    double w1=0.564;
    double w2=0.345;
    double Zc1=0.224;
    double Zc2=0.229;
    double S_old,S_new;
    double Vc1=(Zc1*R*Tc1)/Pc1;
    double Vc2=(Zc2*R*Tc2)/Pc2;
    cout << "Vc1 "  << Vc1 << endl;
    cout << "Vc2 "  <<Vc2 << endl;
    //Calculate Wilson coefficient
    double K12 = (V2 / V1) * exp(-a12 / (R * (T + 273.15)));
    double K21 = (V1 / V2) * exp(-a21 / (R * (T + 273.15)));
    cout << "K12 " << K12 << endl;
    cout << "K21 "<<  K21 << endl;
   // Calculating the activity coefficients
    double Y1 = exp(-log(x1 + x2 * K12) + x2 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
    double Y2 = exp(-log(x2 + x1 * K21) - x1 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
     cout << "Y1" << Y1 << endl;
     cout << "Y2" << Y2 << endl;
    // Calculate B11,B22,B12
    double Tr1=(T+273.15)/Tc1;
    double Tr2=(T+273.15)/Tc2;
     cout << "Tr1 " << Tr1 << endl;
     cout << "Tr2 " << Tr2 << endl;
    // For component 1
    double B_01=0.083-0.422/(pow(Tr1,1.6));
    double B_11=0.139-0.172/(pow(Tr1,4.2));
    double B11=(R*Tc1)*(B_01+w1*B_11)/(Pc1);
    //For component 2
    double B_02=0.083-0.422/(pow(Tr2,1.6));
    double B_12=0.139-0.172/(pow(Tr2,4.2));
    double B22=(R*Tc2)*(B_02+w2*B_12)/(Pc2);
    cout<<"B11 "<<B11<<endl;
    cout<<"B22 "<<B22<<endl;
    // To calculate B12 through mix
    double Tc12=pow((Tc1*Tc2),0.5);
    double Tr12=(T+273.15)/Tc12;
    cout << "Tr12 " << Tr12 << endl;
    cout << "Tc12 "  << Tc12 << endl;
    double w12=(w1+w2)/2;
    double Zc12=(Zc1+Zc2)/2;
    double Vc=(cbrt(Vc1)+cbrt(Vc2));
    cout<<"Vc "<<Vc<<endl;
    double Vc12=pow(Vc/2,3);//Power change
    cout << "Vc1 "  << Vc1 << endl;
    cout << "Vc2 "  << Vc2 << endl;
    cout << "Vc12 " << Vc12 << endl;
    double Pc12=(Zc12*R*Tc12)/Vc12;
    cout << "Pc1 " << Pc1 << endl;
    cout << "Pc2 " << Pc2 << endl;
    cout << "Pc12 " << Pc12 << endl;
    double B_012=0.083-0.422/(pow(Tr12,1.6));
    double B_112=0.139-0.172/(pow(Tr12,4.2));
    double B12=(R*Tc12)*(B_012+w12*B_112)/(Pc12);
    cout << "B12 " << B12 << endl;
    // Calculate Del12
    double del12=2*B12-B11-B22;
     cout << "del12 " << del12 << endl;
    // To calculate Fugacity
    double P_old=((x1*Y1*P1_sat)+(x2*Y2*P2_sat));
    double y1=(x1*Y1*P1_sat)/P_old;
    double y2=(x2*Y2*P2_sat)/P_old;
    cout << y1 << endl;
    cout << y2 << endl;
     double P_new=P_old;
     cout << P_new << endl;
     int c=0;
    do
    { 
          cout << "Outer" << endl;
        P_old=P_new;
       cout << "P_old " << P_old << endl;
    double F1=P1_sat*exp(B11*P1_sat/(R*(T+273.15)))*exp(V1*pow(10,-6)(P_old-P1_sat)/((R*(T+273.15))));
    double F2=P2_sat*exp(B22*P2_sat/(R*(T+273.15)))*exp(V2*pow(10,-6)(P_old-P2_sat)/((R*(T+273.15))));
     phi1=exp(P_old/(R*(T+273.15))*(B11+y2*y2*del12));         
        phi2=exp(P_old/(R*(T+273.15))*(B22+y1*y1*del12));
   double  k1=(Y1*F1)/(phi1*P_old);
         double k2=(Y2*F2)/(phi2*P_old); 
            S_new=k1*x1+k2*x2;
             y1=(k1*x1)/S_new;
         y2=(k2*x2)/S_new;
    do
    {   S_old=S_new;
          cout << "Inner" << endl;
        phi1=exp(P_old/(R*(T+273.15))*(B11+y2*y2*del12));
        phi2=exp(P_old/(R*(T+273.15))*(B22+y1*y1*del12));
      k1=(Y1*F1)/(phi1*P_old);
         k2=(Y2*F2)/(phi2*P_old);
        cout<<"phi="<<phi2<<endl;
         cout<<"k2="<<k1<<endl;
        S_new=k1*x1+k2*x2;
        cout<<S_old<<endl;
        y1=(k1*x1)/S_new;
        y2=(k2*x2)/S_new;
     cout<<"y1= "<<y1<<endl;
        cout<<"y2= "<<y2<<endl;
        

         cout<<S_new<<endl;
    }
    while(fabs(S_new-S_old)>0.001);
    P_new=((x1*Y1*F1)/phi1)+((x2*Y2*F2)/phi2);
    y1=(x1*Y1*F1)/(phi1*P_new);
    y2=(x2*Y2*F2)/(phi2*P_new);
    cout << "P_new " << P_new << endl;
 
}
while(fabs(P_new-P_old)>1);
cout << "P_new is " << P_new/1000 << endl;
 cout << "y1 is" << y1 << endl;
}