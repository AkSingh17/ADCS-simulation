#include<bits/stdc++.h>
#include "programs_new.h"
#define pi 3.14159265359

using namespace std;

double **eci_to_ecef(double **position_eci,double gst_angle)
{ 
    double **dcm_eci2ecef = getzeromatrix(3,3);

    gst_angle = gst_angle*(pi/180.0);

    //transformation matrix calculation
    //markley page 32
    dcm_eci2ecef[0][0] = cos(gst_angle);
    dcm_eci2ecef[0][1] = sin(gst_angle);
    dcm_eci2ecef[0][2] = 0;
    dcm_eci2ecef[1][0] = -sin(gst_angle);
    dcm_eci2ecef[1][1] = cos(gst_angle);
    dcm_eci2ecef[1][2] = 0;
    dcm_eci2ecef[2][0] = 0;
    dcm_eci2ecef[2][1] = 0;
    dcm_eci2ecef[2][2] = 1;
   
    double **position_ecef;
    position_ecef = matrixmultiply(dcm_eci2ecef,position_eci,3,3,1);
    
    for(int i=0;i<3;i++){
        position_ecef[i][0]*=1000;
    }

    free_variable(dcm_eci2ecef,3);

    return position_ecef;
}


double **ecef_to_lla(double **ecef_pos){

    // CODE PICKED UP FROM https://github.com/dhale/idh/blob/master/bench/src/gph/ecef2lla.m

    double **lla_pos = getzeromatrix(3,1);
    double alt,lat,lon;
    
    double radius_earth = 6378137;
    double ecc = 8.1819190842622e-2;

    double radius_sq = radius_earth*radius_earth;
    double e_sq = ecc*ecc;

    double b = sqrt( radius_sq * (1-e_sq));
    double b_sq = b*b;
    double ep = sqrt( (radius_sq - b_sq)/b_sq);

    double p = sqrt(pow(ecef_pos[0][0],2) + pow(ecef_pos[1][0],2));

    double th = atan2(radius_earth*ecef_pos[2][0], b*p);
 
    lon = atan2(ecef_pos[1][0],ecef_pos[0][0])*(180.0/pi);

    lat = atan2((ecef_pos[2][0] + pow(ep,2)*b*pow(sin(th),3) ), (p - e_sq*radius_earth*pow(cos(th),3)))*(180.0/pi);
   
    double N = radius_earth/(sqrt(1-e_sq*pow(sin(lat*(pi/180.0)),2)));
    alt = p / cos(lat*(pi/180.0)) - N;
   
    lla_pos[0][0] = lat;
    lla_pos[1][0] = lon;
    lla_pos[2][0] = alt/1000;
   
    return lla_pos;
}


double **ecef_to_eci(double **ecef_pos,double gst_angle,int m,int n){
    
    //markley page 32
    // calculating the trans matrix to go from eci to ecef
    double **dcm_eci2ecef = getzeromatrix(3,3);

    gst_angle = gst_angle*(pi/180.0);

    dcm_eci2ecef[0][0] = cos(gst_angle);
    dcm_eci2ecef[0][1] = sin(gst_angle);
    dcm_eci2ecef[0][2] = 0;
    dcm_eci2ecef[1][0] = -sin(gst_angle);
    dcm_eci2ecef[1][1] = cos(gst_angle);
    dcm_eci2ecef[1][2] = 0;
    dcm_eci2ecef[2][0] = 0;
    dcm_eci2ecef[2][1] = 0;
    dcm_eci2ecef[2][2] = 1;

    // calculating the inverse of this function for matrix to go from ecef to eci
    double **dcm_ecef2eci;
    dcm_ecef2eci = InverseofMatrix(dcm_eci2ecef,3,3);

    double **pos_eci;
    pos_eci = matrixmultiply(dcm_ecef2eci,ecef_pos,3,3,1);

    free_variable(dcm_ecef2eci,3);
    free_variable(dcm_eci2ecef,3);

    return pos_eci;  
}


double **ned_to_ecef_to_eci(double **magneticfeild_ned, double **pos_lla, double gst_angle){

    //https://www.mathworks.com/help/aeroblks/directioncosinematrixeceftoned.html
    // taking in values in degrees and changing to radians
    double **lla_rad = getzeromatrix(3,1);

    lla_rad[0][0] = pos_lla[0][0]*(pi/180.0);
    lla_rad[1][0] = pos_lla[1][0]*(pi/180.0);
    lla_rad[2][0] = pos_lla[2][0];
    
    // rotation matrix calculation for ecef to ned 
    double **ecef2ned = getzeromatrix(3,3);
    ecef2ned[0][0]=-sin(lla_rad[0][0])*cos(lla_rad[1][0]);
    ecef2ned[0][1]=-sin(lla_rad[0][0])*sin(lla_rad[1][0]);
    ecef2ned[0][2]=cos(lla_rad[0][0]);
    ecef2ned[1][0]=-sin(lla_rad[1][0]);
    ecef2ned[1][1]=cos(lla_rad[1][0]);
    ecef2ned[1][2]=0;
    ecef2ned[2][0]=-cos(lla_rad[0][0])*cos(lla_rad[1][0]);
    ecef2ned[2][1]=-cos(lla_rad[0][0])*sin(lla_rad[1][0]);
    ecef2ned[2][2]=-sin(lla_rad[0][0]);

    double **ned2ecef, **pos_ecef;
    double **pos_eci;

    // inverse the matrix for ned to ecef
    ned2ecef = InverseofMatrix(ecef2ned,3,3);
    
    double dcmdet = DeterminantofMatrix(ecef2ned);

    pos_ecef = matrixmultiply(ned2ecef,magneticfeild_ned,3,3,1);
    
    //finally converting it from ecef to eci
    pos_eci = ecef_to_eci(pos_ecef, gst_angle,3,1);

    free_variable(lla_rad,3);
    free_variable(ecef2ned,3);
    free_variable(ned2ecef,3);
    free_variable(pos_ecef,3);

    return pos_eci;
}


double **eci_orbit(double **keplerian_elements){
    
    double **trans = getzeromatrix(3,3);
    double i,ohm,agp,TAnrp;

    i = keplerian_elements[2][0]; //inclination
    ohm = keplerian_elements[3][0]; //raan
    agp = keplerian_elements[4][0]; //argument of perigee
    TAnrp = keplerian_elements[5][0]; //true anomaly

    trans[0][0] = -sin(agp+TAnrp)*cos(ohm)+-cos(i)*sin(ohm)*cos(agp+TAnrp);
    trans[0][1] = -sin(agp+TAnrp)*sin(ohm)+cos(i)*cos(ohm)*cos(agp+TAnrp);
    trans[0][2] = cos(agp+TAnrp)*sin(i);
    trans[1][0] = -sin(i)*sin(ohm);
    trans[1][1] = sin(i)*cos(ohm);
    trans[1][2] = -cos(i);
    trans[2][0] = -cos(agp+TAnrp)*cos(ohm)+cos(i)*sin(ohm)*sin(agp+TAnrp);
    trans[2][1] = -cos(agp+TAnrp)*sin(ohm)-cos(i)*cos(ohm)*sin(agp+TAnrp);
    trans[2][2] = -sin(agp+TAnrp)*sin(i);

    return trans;
}

double **eci_body(double **eci_vec, double** q_bi_1){

    int i;
    double** temp_quat = getzeromatrix(4,1);
    double** body_vec = getzeromatrix(3,1);
        
    for(i=0; i<3; i++){
        temp_quat[i+1][0] = eci_vec[i][0];
    }
    
    double ** q_bi_1_temp, **q_ib_1, **temp_quat_2;
    q_bi_1_temp = quatmultiply(q_bi_1,temp_quat);
    q_ib_1 = quatinv(q_bi_1);
    temp_quat_2 = quatmultiply(q_bi_1_temp,q_ib_1);
        
    for(i=0;i<3;i++){
        body_vec[i][0] = temp_quat_2[i+1][0];
    }

    free_variable(temp_quat,4);
    free_variable(q_ib_1,4);
    free_variable(q_bi_1_temp,4);
    free_variable(temp_quat_2,4);

    return body_vec;
}