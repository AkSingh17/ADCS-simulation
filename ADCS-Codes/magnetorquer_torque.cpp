#include <bits/stdc++.h>
#include "programs_new.h"
#define pi 3.14159265359

double** magnetorquer(double** rw_w_cur, double** mag_body, double** Tc){
    
    double **M, **torque, **temp, **dump, **mag_torque;
    double flagm;
    int i;

    torque = getzeromatrix(3,1);
    dump = getzeromatrix(3,1);
    mag_torque = getzeromatrix(3,1);

    dump[0][0] = 1000 * (2 * pi / 60);
    dump[1][0] = 1000 * (2 * pi / 60);
    dump[2][0] = 1000 * (2 * pi / 60);

    double B_sq = pow(norm(mag_body), 2);

    for (i = 0; i < 3; i++){

        if (rw_w_cur[i][0] > 0 )
            torque[i][0] = -0.1*(4.5e-6)*(rw_w_cur[i][0] - dump[i][0]);

        else if (rw_w_cur[i][0] < 0 )
            torque[i][0] = -0.1*(4.5e-6)*(rw_w_cur[i][0] + dump[i][0]);

        else if (rw_w_cur[i][0] == 0 )
            torque[i][0] = 0;
    }

    for(i=0;i<3;i++){
        torque[i][0] = torque[i][0] / B_sq;
    }

    flagm = 1;

    M = CrossProduct(mag_body, torque);
    temp = CrossProduct(M, mag_body);

    for(i=0;i<3;i++){   
        if(temp[i][0]*torque[i][0]<0){
            flagm = 0;
            break;
        }
    }

    for(i=0;i<3;i++){
        M[i][0] = flagm*M[i][0];
    }

    for(i = 0; i < 3; i++){
        if (M[i][0] > 0.107472)
            M[i][0] = 0.107472;
        
        else if (M[i][0] < -0.107472)
            M[i][0] = -0.107472;
    }

    mag_torque = CrossProduct(M, mag_body);

    free_variable(torque,3);
    free_variable(temp,3);
    free_variable(dump,3);
    free_variable(M,3);

    return mag_torque;
}