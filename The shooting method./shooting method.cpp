#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include<numeric>
#include<valarray>
#include<fstream>
using namespace std;
double max(double x1, double x2)
{
    if (abs(x1)>abs(x2))
    {
        return (x1-x2)/abs(x1);
    }
    else
    {
        return (x1-x2)/abs(x2);
    }
}
double Numerov(double E)
{
    double a = 5.0,V_0 = 5.0;
    double h = 0.00001;
    double delta = 0.1; // the shooting angel
    double K1,K2,K3; // temp variables
    double r_max=1.0;
    double r_m = 2.0/a;
    double delta_deri=0;
    double r_min = -1.0;
    int N = 2.0/h+1, N1 = (int) ((r_m-r_min)/h)+1, N2 = (int) ((r_max-r_m)/h)+1; // the iteration step;
    double *phi1;
    double *phi2;
    double *phi;
    double gamma = ((2*pow(a,2)*V_0));// the auxiliary variable
    phi1 = new double[N1]; //From r_min to r_m, the needed steps
    phi2 = new double[N2];//From r_max to r_m, the needed steps
    phi = new double[N];
    phi1[0] = 0;
    phi2[N2-1] = 0;// correspond to r_max
    double *K_iter;
    K_iter = new double[N];
    double *x_in;
    x_in = new double[N];
    double *r_in;
    r_in = new double[N];
    double *V_in;
    V_in = new double[N];
    int i;
    for (i=0; i<N; i++)
    {
        x_in[i] = (-1+h*i)*a;
    }
//    for (i=0; i<N; i++)
//    {
//        r_in[i] = x_in[i]/a;
//    }
    for (i=0; i<N; i++)
    {
        V_in[i] = (-pow(x_in[i],2)+pow(x_in[i],4)/20.0)/V_0;
    }
    for (i=0; i<N; i++)
    {
        K_iter[i] = gamma*(E/V_0-V_in[i]);
    }
    phi1[1] = phi1[0]+delta*h+(1.0/6.0)*pow(h,3)*(-K_iter[0])*delta;
    phi2[N2-2] = -delta*h-(1.0/6.0)*pow(h,3)*(-K_iter[N-1])*delta;
    for (i=1; i<N1-1; i++)
    {
        K1 = 1+pow(h,2)*K_iter[i+1]/12.0;
        K2 = 1-pow(h,2)*K_iter[i]*5.0/12.0;
        K3 = 1+ pow(h,2)*K_iter[i-1]/12.0;
        phi1[i+1] = phi1[i]*(2*K2/K1) - phi1[i-1]*K3/K1;
        phi[i] = phi1[i];
    }
    phi[0] = phi1[0];
    phi[N1-1] = phi1[N1-1];
    for (i=N2-2; i>0; i--)
    {
        K1 = 1+pow((h),2)*K_iter[i+N-N2+1]/12.0;
        K2 = 1-pow((h),2)*K_iter[i+N-N2]*5.0/12.0;
        K3 = 1+ pow((h),2)*K_iter[i+N-N2-1]/12.0;
        phi2[i-1] = -phi2[i+1]*(K1/K3)+2*phi2[i]*(K2/K3);
    }
    double lambda = phi1[N1-1]/phi2[0];
    for (i=0; i<N2; i++)
    {
        phi[N1+i] = lambda*phi2[i];
    }
    delta_deri = max((phi1[N1-1]-phi1[N1-2])/h , lambda*(phi2[1]-phi2[0])/h);
    printf("The data is %f\n",delta_deri);
    fstream output_stream;
    output_stream.open("/Users/phi-0.564230.txt",ios::out | ios::app);
    for (i=0;i<N;i++)
    {
        output_stream << phi[i] << endl;
    }

    delete[]phi1;
    delete[]phi2;
    delete[]K_iter;
    delete[]V_in;
    delete[]x_in;
    return delta_deri;
}


int main(int argc, const char * argv[])
{
    int i,k;
    double En;
    double h_s = pow(10,-6);
    int N_s = 0.02/h_s;
    double *deltaE;
    deltaE = new double[N_s];
    Numerov(-0.564230);
    //try to find the accurate energy
//    for (i=0; i<N_s; i++)
//    {
//        En = -0.685+i*h_s;
//        deltaE[i] = Numerov(En);
////        printf("this is delta %f\n",deltaE[i]);
//        if (fabs(deltaE[i])<pow(10,-5))
//        {
//            printf("The data%f\n",En);
//        }
//    }
//    fstream output_stream;
//    output_stream.open("/Users/lichan/Desktop/研究生一年级/计算物理/Chapter8 Shooting method/deltaE.txt",ios::out | ios::app);
//    for (k=0;k<N_s;k++)
//    {
//        output_stream << deltaE[k] << endl;
//    }
    delete[]deltaE;
    return 0;
}
