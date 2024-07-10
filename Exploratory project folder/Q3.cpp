#include<bits/stdc++.h>
using namespace std;
int main()
{
    double t=100,y1=0.5,y2=0.5;
    double P1sat,P2sat;
    double Y1old,Y2old;
    double tnew;
    double P;
    double x1,x2;
    double Y1new=1,Y2new=1;
     P1sat=exp(16.678-(3640.2/(t+219.61)));
    P2sat=exp(16.2887-(3816.44/(t+227.02)));
    double Pnew=1/((y1/(Y1new*P1sat))+(y2/(Y2new*P2sat)));
                //  double t1_sat = 3640.2 / (16.6788 - log(P)) - 219.61;
        //  double t2_sat = 3816.44 / (16.2887 - log(P)) - 227.02;
          // tnew = x1*Y1old * t1_sat + x2*Y2old * t2_sat;
       
    double Pold;
    double a12 = 437.98; // in cal/mol
    double a21 = 1238;   // in cal/mol
    double V1 = 76.92;   // cm^3/mol
    double V2 = 18.07;   // cm^3/mol
    double r = 1.987;
    do
    {
      Pold=Pnew;
      
         cout<<"Outer"<<endl;
        // told=tnew;
    //     P1sat=exp(16.678-(3640.2/(t+219.61)));
    // P2sat=exp(16.2887-(3816.44/(told+227.02)));
 
        // P=1/((y1/(Y1new*P1sat))+(y2/(Y2new*P2sat)));
        
        do
        {
           cout<<"Inner"<<endl;

            Y1old=Y1new;
            Y2old=Y2new;

            x1=(y1*Pold)/(Y1old*P1sat);
            x2=(y2*Pold)/(Y2old*P2sat);
            x1=x1/(x1+x2);
            x2=x2/(x1+x2);
        //     double t1_sat = 3640.2 / (16.6788 - log(P)) - 219.61;
        //  double t2_sat = 3816.44 / (16.2887 - log(P)) - 227.02;
          // tnew = y1 * t1_sat + y2 * t2_sat;
             
    // Calculating Wilson's Coefficients
    double K12 = (V2 / V1) * exp(-a12 / (r * (t + 273.15)));
    double K21 = (V1 / V2) * exp(-a21 / (r * (t + 273.15)));
    // Calculating the activity coefficients
     Y1new = exp(-log(x1 + x2 * K12) + x2 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
     Y2new = exp(-log(x2 + x1 * K21) - x1 * (K12 / (x1 + x2 * K12) - (K21 / (x2 + x1 * K21))));
           
        } while (fabs(Y1new-Y1old)>0.01 && fabs(Y2new- Y2old)>0.01);
            
    //         P1sat=exp(16.678-(3640.2/(tnew+219.61)));
    // P2sat=exp(16.2887-(3816.44/(tnew+227.02)));
    // tnew = 3640.2 / (16.6788 - log(P1sat)) - 219.61;
    cout<<"P1sat= "<<P1sat<<endl;
      cout<<"P2sat= "<<P2sat<<endl;
  Pnew=1/((y1/(Y1new*P1sat))+(y2/(Y2new*P2sat)));
    
        }
     while (fabs(Pnew-Pold)>0.1);
    //  cout <<  1/((y1/(Y1new*P1sat))+(y2/(Y2new*P2sat))) << endl;
    cout << "P_OLD IS" << Pold << endl;
     cout << "P_NEW IS" << Pnew << endl;

     
}