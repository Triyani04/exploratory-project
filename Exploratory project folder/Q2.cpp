#include <bits/stdc++.h>
using namespace std;
int main()
{
    double T_old; int c=0;
    double P = 101.325; // Given temperature is in Degree Celcius
    double x1 = 0.5;    // Mole fraction
    double x2 = 1 - x1;
    double r = 1.987;    // cal/mol K
    double a12 = 437.98; // in cal/mol
    double a21 = 1238;   // in cal/mol
    double V1 = 76.92;   // cm^3/mol
    double V2 = 18.07;   // cm^3/mol
    double t1_sat = 3640.2 / (16.6788 - log(101.325)) - 219.61;
    double t2_sat = 3816.44 / (16.2887 - log(101.325)) - 227.02;
    double T = x1 * t1_sat + x2 * t2_sat;
    // Calculating Wilson's Coefficients
    double K12 = (V2 / V1) * exp(-a12 / (r * (T + 273.15)));
    double K21 = (V1 / V2) * exp(-a21 / (r * (T + 273.15)));
    // Calculating the activity coefficients
    double Y1 = exp(-log(x1 + x2 * K12) + x2 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
    double Y2 = exp(-log(x2 + x1 * K21) - x1 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
    double P1_sat = exp(16.678 - (3640.2 / (T + 219.61)));   // Calculating p1_sat
    double P2_sat = exp(16.2887 - (3816.44 / (T + 227.02))); // Calculating p2_sat
    double P1_upsat = P / (x1 * Y1 + (x2 * Y2 * P2_sat / P1_sat));
    double T_new = 3640.2 / (16.6788 - log(P1_upsat)) - 219.61;
    do
    {
        T_old = T_new; c++;
        P1_sat = exp(16.678 - (3640.2 / (T_new + 219.61)));   // Calculating p1_sat Updated
        P2_sat = exp(16.2887 - (3816.44 / (T_new + 227.02))); // Calculating p2_sat Updated
        K12 = (V2 / V1) * exp(-a12 / (r * (T_new + 273.15)));
        K21 = (V1 / V2) * exp(-a21 / (r * (T_new + 273.15)));
        Y1 = exp(-log(x1 + x2 * K12) + x2 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
        Y2 = exp(-log(x2 + x1 * K21) - x1 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
        P1_upsat = P / (x1 * Y1 + (x2 * Y2 * P2_sat / P1_sat));
        T_new = 3640.2 / (16.6788 - log(P1_upsat)) - 219.61;
    } while (fabs(T_new - T_old) > 0.001);
    cout << T_new << " " << c;
}