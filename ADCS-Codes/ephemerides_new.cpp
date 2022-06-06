#include<bits/stdc++.h>
using namespace std;

double getJulianDate(double year,double month, double date, double hour, double min, double sec){
    long A,B,C,E;
    double JD;
    if(month==1 || month==2)
    {
        year=year-1;
        month=month+12;
    }
    A=(int)(year/100.0);
    B=2-A+(int)A/4;
    if(year<0)
    {
        C=(int)(365.25*year-0.75);
    }
    C=(int)(365.25*year);
    E=(int)(30.6001*(month+1));
    JD=B+C+date+E+1720994.5;
    JD=JD+ (hour)/24.0 + month/1440.0 + sec/86400.0;
    //cout<<"hi"<<endl;
    return JD;
}

double jdut2gmst(double JD)
{
    double H;
    double diff=fabs(floor(JD)-JD);
    if(diff<=0.5)
    {
        H=(diff+0.5)*24;
    }
    else
    {
        H=(diff-0.5)*24;
    }
    double JD0=JD-H/24;
    double D=JD-2451545.0;
    double D0=JD0-2451545.0;
    double T=D/36525;
    double GMST= 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T*T;
    long quotient=(long)GMST/24;
    //cout<<quotient<<endl;
    double GMST_0_24=fabs(24*quotient-GMST)*15.0;
    // cout<<GMST_0_24<<endl;
    return GMST_0_24;
}