#include<bits/stdc++.h>
#include "programs_new.h"

using namespace std;

ofstream pow_mag("pow_magnetorquer.txt");

double **detumble(double **mag_body, double **w_bi, double **moi_satellite){
    
    /**the bdot algorithm is a proportional controller**/

    ///vector initialisation///
    double **Tc;
    double **tau = getzeromatrix(3,1);
    double **moment;
    double power_magnetorquer;

   ///calculating detumbling torque
    for(int i=0;i<3;i++){
        tau[i][0] = moi_satellite[i][i]*(w_bi[i][0]);
    }

    ///cross product control law
    moment = CrossProduct(tau,mag_body);
    for(int i=0;i<3;i++){
        moment[i][0] = moment[i][0]/pow(norm(mag_body),2);
    }
    
    ///moment saturation
    for(int i=0;i<3;i++){

        if (moment[i][0] > 0.107472){
            moment[i][0] = 0.107472;
        }
        else if (moment[i][0] < -0.107472){
            moment[i][0] = -0.107472;
        }
    }

    ///M X B is the torque delivered by magnetorquers
    Tc = CrossProduct(moment,mag_body);

    free_variable(tau,3);
    free_variable(moment,3);

    return Tc;
}