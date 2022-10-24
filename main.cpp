#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;
double q(double x);
double p(double x);
double f(double x);
double u(double x);
double a(double x);
double b(double x);
double c(double x);
double NAlpha(double x,double previousalpha);
double NBeta(double x,double previousalpha,double previousbeta);

const double e=0.0;
const double e1=-sqrt(3);
const double e2=2*sqrt(3);
const double d=4.0;
const double d1=sqrt(2);
const double d2=-sqrt(2);
const int N=1000;
const double sb=1.0;
const double sa=0.0;
//для вычесления
const double h=(sb-sa)/N;
const double r1=-d2/(h*d1-d2);
const double s1=(h*d)/(d1*h-d2);
const double r2=e2/(e1*h+e2);
const double s2=e*h/(e1*h+e2);
//cd 'D:\KPFU C++\CryMethod'
//plot "CryMethod.dat" title "graphics" linetype 6 linecolor 4 with linespoints
// вариант 23
int main(){

    vector<double> X(N+1); 
    vector<double> U(N+1);
    for(int i=0;i<N;i++){
        X[i]=sa+h*i;
        U[i]=u(X[i]);
    }
    X[N]=sb;
    U[N]=u(X[N]);
    vector<double>Y(N+1);
    double* Alpha=new double[N+1];
    double* Betta=new double[N+1];
    Alpha[1]=r1;
    Betta[1]=s1;
    for(int i=2;i<=N;i++){
        Alpha[i]=NAlpha(X[i],Alpha[i-1]);
        Betta[i]=NBeta(X[i],Alpha[i-1],Betta[i-1]);
    }
    Y[N]=(Betta[N]*r2+s2)/(1-Alpha[N]*r2);
    for(int i=N;i>0;i--){
        Y[i-1]=Alpha[i]*Y[i]+Betta[i];
    }
    delete[] Alpha;
    delete[] Betta;
    double max=0;
    for(int i=0;i<=N;i++){
        if(max<abs(U[i]-Y[i])){
            max=abs(U[i]-Y[i]);
        }
    }
    cout<<max;
/*    RGBABitmapImageReference *imageRef= CreateRGBABitmapImageReference();
    DrawScatterPlot(imageRef,600,400,&X,&Y);
    vector<double> *pngData=ConvertToPNG(imageRef->image);
    WriteToFile(pngData,"plot.png");
    DeleteImage(imageRef->image);*/
    fstream fout;
    fout.open("CryMethod.dat",ios::out);
    if(fout){
        for(int i=0;i<=N;i++){
            fout<<X[i]<<"    ";
            fout<<abs(U[i]-Y[i])<<endl;
        }
        fout<<X[N];
        fout<<Y[N];
    }
    cout<<sa;
    return 0;
}
double NAlpha(double x,double previousalpha){
    return c(x)/(b(x)-previousalpha*a(x));
}
double NBeta(double x,double previousalpha,double previousbeta){
    return (previousbeta*a(x)-f(x))/(b(x)-previousalpha*a(x));
}
double c(double x){
    return (1/(h*h)+p(x)/(2*h));
}
double b(double x){
    return (2/(h*h)-q(x));
}
double a(double x){
    return (1/(h*h)-p(x)/(2*h));
}
double p(double x){
    return 2.0/(x+2);
}
double q(double x){
    return -3.0/(x+2)/(x+2);
}
double f(double x){
    return 3/(sqrt(x+2));
}
double u(double x){
    return 4.0*(x+2)*sqrt(x+2);
}
