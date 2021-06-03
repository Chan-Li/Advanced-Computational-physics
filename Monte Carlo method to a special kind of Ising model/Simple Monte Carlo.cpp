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
# define N 100
# define M 100000
//This is the random number generator-double
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / ((double)RAND_MAX + 1.0);
    return fMin + f * (fMax - fMin);
}
double Initial()
{
    double ran;
    double temp;
    ran = fRand(-1.0,1.0);
    if (ran>0)
    {
            temp = 0.5;
    }
    else
    {
        temp = -0.5;
            
    }
    return temp;
}
// Initialization
double D_Engergy1(int Flap,double *S)
{
    int sum1=0;
    if (Flap==0)
    {
        sum1=2*S[Flap]*(S[N-1]+S[1]);
    }
    else if (Flap==N-1)
    {
        sum1=2*S[Flap]*(S[N-2]+S[0]);
    }
    else{
        sum1 = 2*S[Flap]*(S[Flap-1]+S[Flap+1]);
    }
    return sum1;
}

double D_Engergy2(int Flap,double *S)
{
    int sum2=0;
    if (Flap==0 || Flap ==1)
    {
        sum2=2*S[Flap]*(S[N+Flap-2]+S[Flap+2]);
    }
    else if (Flap==N-1 || Flap == N-2)
    {
        sum2=2*S[Flap]*(S[Flap+2-N]+S[Flap-2]);
    }
    else{
        sum2 = 2*S[Flap]*(S[Flap-2]+S[Flap+2]);
    }
    return sum2;
}
double magnetization(double *S)
{
    int i;
    double mag=0;
    for (i=0; i<N; i++)
    {
        mag = mag+S[i]*pow(-1, (i+1));
        //mag = mag+S[i];
    }
    return mag;
}


double magnetization2(double *S)
{
    int i;
    double mag1=0,mag2=0,m3=0;
    for (i=0; i<N; i++)
    {
        mag1 = mag1+S[i]*pow(-1, (int)((i)/2));
    }
    for (i=0; i<N-1; i++)
    {
        mag2 = mag2+S[i+1]*pow(-1, (int)((i)/2));
    }
    mag2+=S[0]*(-1);
    m3 = fabs(mag1)+fabs(mag2);
    return m3;
}
double Hamilton(double *S,double J_1,double J_2,double beta)
{
    int i,j;
    double E1=0,E2=0,E_total=0;
    for(i=0;i<N-1;i++)
    {
        E1 += -J_1*(S[i]*S[i+1]);
    }
    E1+=-J_1*S[N-1]*S[0];
    for(j=0;j<N-2;j++)
    {
        E2 += -J_2*S[j]*S[j+2];
    }
    E2 +=-J_2*(S[N-2]*S[0]+S[N-1]+S[1]);
    E_total = E1+E2;
    return E_total;
}
double C_v(double beta,double H1,double H2)
{
    double cv;
    cv = beta*beta*(H1-H2)/N;
    return cv;
}

double Chi_(double beta,double M1,double M2)
{
    double Chi;
    Chi = (beta/N)*(M1-M2);
    return Chi;
}

double *Update(double S[],double J_1,double J_2,double beta,int num)
{
    double delta_E=0;
    double ran;
    int i,flap,j;
    double *S_now;
    S_now = new double[N];
    for (j=0;j<N;j++)
    {
        S_now[j] = S[j];
    }
    //S_now = S;
    for (i=0;i<num;i++)
    {
        delta_E=0;
        flap = (int)(fRand(0,N));
        S_now[flap] = -1*S_now[flap];
        delta_E = beta*J_1*D_Engergy1(flap,S_now)+beta*J_2*D_Engergy2(flap,S_now);
        if (fabs(delta_E)>fabs(beta*J_1*4+beta*J_2*4))
        {
            printf("Something is wrong!%f\n",delta_E);
        }
        else
        {
            ran = fRand(0.0,1.0);
            if (exp(delta_E)>ran)
            {
                S_now = S_now;
            }
            else
            {
                S_now[flap] = -1*S_now[flap];//翻转回去
            }
        }
        
    }
    return S_now;
}



int main()
{
    int i,j,m,k;
    double magnet=0,magnet2=0,M1=0,M2=0,chii,chii2,M3=0,M4=0;
    double Tot_E=0,H1=0,H2=0,CV,temp;
    double J_1 = -1.0,beta,J_2;
    double *S_Intial;
    S_Intial = new double[N];
    double *S_new;
    S_new=new double(N);
    double *M_data;
    M_data = new double(360);
    for (m=0; m<60; m++)//温度
    {
        Tot_E = 0;
        magnet=0;
        magnet2=0;
        //J_2 = (J_1)*(m)*0.01;
        J_2 =-0.4;
        temp=(m+1)*fabs(J_1)*0.1;
        //temp=0.05;
        beta = 1/temp;//定义温度
        CV=0;
        H2=0;
        H1=0;
        M1=0;
        M2=0;
        M3=0;
        M4=0;
        chii=0;
        chii2=0;
        srand((unsigned)time(NULL)+123 * (m + 12345));
        for (i=0;i<N;i++)
        {
            S_Intial[i] = Initial();
        }
        S_new = Update(S_Intial,J_1,J_2,beta,100000);
        for (j=0; j<M; j++)
        {
            S_new = Update(S_new,J_1,J_2,beta,5000);
//            if (j==M-1)
//            {
//                for (i=0; i<N; i++)
//                {
//                    printf("This is an example%f\n",S_new[i]);
//                }
//            }
            H1+= pow((Hamilton(S_new,J_1,J_2,beta)),2);
            H2+= Hamilton(S_new,J_1,J_2,beta);
            Tot_E += Hamilton(S_new,J_1,J_2,beta)/(fabs(J_1)*N);
            magnet = magnet+fabs(magnetization(S_new)/N);
            magnet2 = magnet2+fabs(magnetization2(S_new)/N);
            M1+=pow(magnetization(S_new),2);
            M2+=fabs(magnetization(S_new));
            M3+=pow(magnetization2(S_new),2);
            M4+=fabs(magnetization2(S_new));
        }
        CV = C_v(beta,H1/M, pow((H2/M),2));
        chii = Chi_(beta,M1/M,pow((M2/M),2));
        chii2 = Chi_(beta,M3/M,pow((M4/M),2));
        M_data[6*m] = (magnet/M);
        M_data[6*m+1] = (magnet2/M);
        M_data[6*m+2] = (Tot_E/M);
        M_data[6*m+3] = (CV);
        M_data[6*m+4] = (chii);
        M_data[6*m+5] = (chii2);
        printf("Iteration %d\n",m);
        printf("T %f\n",temp);
        printf("This is the M1 %f\n",(M_data[6*m]));
        printf("This is the M2 %f\n",(M_data[6*m+1]));
        printf("This is the E %f\n",(M_data[6*m+2]));
        printf("This is the CV %f\n",(M_data[6*m+3]));
        printf("This is the chii1 %f\n",(M_data[6*m+4]));
        printf("This is the chii2 %f\n",(M_data[6*m+5]));
    }
    fstream output_stream;
    output_stream.open("/Users/lichan/Desktop/研究生一年级/计算物理/Chapter2 Monte Carlo/entropyJ204TC2.txt",ios::out | ios::app);
    for (k=0;k<360;k++)
    {
        output_stream << M_data[k] << endl;
    }
    delete[]S_new;
    delete[]S_Intial;
    //delete M_data;
    return 0;
}



