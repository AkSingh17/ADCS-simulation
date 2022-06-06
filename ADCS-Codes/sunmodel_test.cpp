#include <bits/stdc++.h>
#include "programs_new.cpp"
using namespace std;
#define pi 3.14159265359

double** sun_model_purana(double jdut, double** r)
{
    int i;
    double Msun, Vsun, e1, Resmag;
    double** sun = getzeromatrix(3,1);
    double** rss = getzeromatrix(3,1);
    jdut -= 2451545;

    Msun = 357.528 + (0.9856003 * jdut);
    Vsun = 280.460 + (0.9856474 * jdut) + (1.915 * sin(Msun * pi / 180)) + (0.020 * sin(((2 * Msun * pi) / 180)));
    e1 = 23.439 - (0.0000004 * jdut);

    sun[0][0] = cos((Vsun * pi / 180));
    sun[1][0] = cos((e1 * pi / 180)) * sin((Vsun * pi / 180));
    sun[2][0] = sin((e1 * pi / 180)) * sin((Vsun * pi / 180));

    Resmag = (1.00014 - 0.01671 * cos(Msun * pi / 180) - 0.00014 * cos(2 * Msun * pi / 180)) * 149597870.700;

    for (i = 0; i < 3; i++)
        rss[i][0] = Resmag * sun[i][0] - r[i][0];

    double** Res = getzeromatrix(3,1);
    for (i = 0; i < 3; i++)
        Res[i][0] = Resmag * sun[i][0];

    double value = DotProduct(Unit_vector(Res,3,1), Unit_vector(r,3,1));
    double dat_1 = 1 - (value * value);
    double dat_2 = sqrt(dat_1) * norm(r);

    double eclipse;
    if (dat_2 < 6378.137 && value < 0)
        eclipse = 0;
    else
        eclipse = 1;

    double** send = getzeromatrix(4,1);;
    send[0][0] = rss[0][0];
    send[1][0] = rss[1][0];
    send[2][0] = rss[2][0];
    send[3][0] = eclipse;

    return send;
}

double **sun_model_naya(double julian_date,double **position){
    
    int i;
    double Msun, Vsun, e1, dist_sun_earth, solar_dist, sat_dist;
    double **sun, **rss, **res, **send;
    sun = getzeromatrix(3,1);
    rss = getzeromatrix(3,1);
    res = getzeromatrix(3,1);
    send = getzeromatrix(4,1);
    julian_date -= 2451545;

    //Calculating Mean Anomaly of the sun in radians.
    Msun = 357.528 + (0.9856003 * julian_date);

    //To ensure Msun <360deg
    if(Msun>360){
        int integer_part=Msun/360;
        double frac_part = Msun - 360*(integer_part);
        Msun = frac_part;
    }

    //Calculating the Longitude of the ecliptic in degrees
    Vsun = 280.460 + (0.9856474 * julian_date) + (1.915 * sin(Msun * pi / 180)) + (0.0199 * sin(((2 * Msun * pi) / 180)));
    
    //To ensure Vsun <360deg
    // if(Vsun>360){
    //     int integer_part=Vsun/360;
    //     double frac_part = Vsun - 360*(integer_part);
    //     Vsun = frac_part;
    // }
    
    //Obliquity of the ecliptic
    e1 = 23.439 - (0.0000004 * julian_date);
    
    //Earth to sun unit vector
    sun[0][0] = cos((Vsun * pi / 180));
    sun[1][0] = cos((e1 * pi / 180)) * sin((Vsun * pi / 180));
    sun[2][0] = sin((e1 * pi / 180)) * sin((Vsun * pi / 180));

    //Computing the distance between the sun and the Earth in AU
    dist_sun_earth = (1.00014 - 0.01671 * cos(Msun * pi / 180) - 0.00014 * cos(2 * Msun * pi / 180));
    dist_sun_earth *= 149597870.700; //converting to kms

    for (i = 0; i < 3; i++)
        rss[i][0] = dist_sun_earth * sun[i][0] - position[i][0];
    
    /// checking for ecllipse
    for (i = 0; i < 3; i++){
        res[i][0] = dist_sun_earth * sun[i][0];
    }
    solar_dist = norm(res);
    for(i=0;i<3;i++){
        res[i][0] = res[i][0]/solar_dist;
    }
    sat_dist = norm(position);
    double value = DotProduct(res,position);
    double dat_1 = 1 - (value * value);
    double dat_2 = sqrt(dat_1) * norm(position);
    double eclipse;
    if (dat_2 < 6378.137 && value < 0)
        eclipse = 0;
    else
        eclipse = 1;

    send[0][0] = eclipse;
    send[1][0] = rss[0][0];
    send[2][0] = rss[1][0];
    send[3][0] = rss[2][0];
    
    free_variable(sun,3);
    free_variable(rss,3);
    free_variable(res,3);

    return send;
}

int main(){

    double** sun_vector;
    double** position = getzeromatrix(3,1);

    double jd = 2.45673e+06;

    position[0][0] =  -5236.84633;
    position[1][0] =  4124.17773;
    position[2][0] =  -1262.94137;

    sun_vector = sun_model_purana(jd, position);

    // cout<<sun_vector[0][0]<<" "<<sun_vector[1][0]<<" "<<sun_vector[2][0]<<" "<<sun_vector[3][0]<<endl;

    return 0;
}