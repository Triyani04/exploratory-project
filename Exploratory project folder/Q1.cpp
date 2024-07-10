#include <bits/stdc++.h>
using namespace std;
int main()
{
    double temp = 100; // Given temperature is in Degree Celcius
    double x1 = 0.5;   // Mole fraction
    double x2 = 1 - x1;
    double r = 1.987;                                           // cal/mol K
    double P1_sat = exp(16.678 - (3640.2 / (temp + 219.61)));   // Calculating p1_sat
    double P2_sat = exp(16.2887 - (3816.44 / (temp + 227.02))); // Calculating p2_sat
    // Calculating Wilson's Coefficients
    double a12 = 437.98; // in cal/mol
    double a21 = 1238;   // in cal/mol
    double V1 = 76.92;   // cm^3/mol
    double V2 = 18.07;   // cm^3/mol
    double K12 = (V2 / V1) * exp(-a12 / (r * (temp + 273.15)));
    double K21 = (V1 / V2) * exp(-a21 / (r * (temp + 273.15)));
    // Calculating the activity coefficients
    double Y1 = exp(-log(x1 + x2 * K12) + x2 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
    cout << Y1 << endl;
    double Y2 = exp(-log(x2 + x1 * K21) - x1 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
    cout << Y2 << endl;
    // Calculating Final Pressure
    //  P=x1*Y1*P1_sat+x2*Y2*P2_sat   (using formula)
    double P = x1 * Y1 * P1_sat + x2 * Y2 * P2_sat;
    cout << "Final Pressure is : " << P << " Kpa";
}