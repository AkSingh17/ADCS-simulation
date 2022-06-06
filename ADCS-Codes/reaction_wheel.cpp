#include "programs_new.h"

double**  reaction_wheels(double** tc, double rw_moi, double saturation_high, double saturation_low, double** rw_w){

    double **delta_rw_w,**true_delta_rw_w,**w_torque,**rw_w_new,**rw_state;
    double rw_moi_inv = 1/rw_moi;

    rw_state = getzeromatrix(6,1);

    delta_rw_w = MatrixScalarMultiply(tc,rw_moi_inv,3,1);
    rw_w_new = MatrixAdd(rw_w,delta_rw_w,3,1);

    for(int i=0; i<3;i++){
        if(rw_w_new[i][0] > saturation_high){
            rw_w_new[i][0] = saturation_high;
        }
        else if((rw_w_new[i][0] < saturation_low) && (rw_w_new[i][0] > 0)){
            rw_w_new[i][0] = saturation_low;
        }
        else if(rw_w_new[i][0] < -1*saturation_high){
            rw_w_new[i][0] = -1*saturation_high;
        }
        else if(rw_w_new[i][0] > -1*saturation_low && rw_w_new[i][0] < 0){
            rw_w_new[i][0] = -1*saturation_low;
        }
    }

    true_delta_rw_w = MatrixSubtract(rw_w_new,rw_w,3,1);
    w_torque = MatrixScalarMultiply(true_delta_rw_w,rw_moi,3,1);

    rw_state[0][0] = w_torque[0][0];
    rw_state[1][0] = w_torque[1][0];
    rw_state[2][0] = w_torque[2][0];
    rw_state[3][0] = rw_w_new[0][0];
    rw_state[4][0] = rw_w_new[1][0];
    rw_state[5][0] = rw_w_new[2][0];

    free_variable(delta_rw_w,3);
    free_variable(true_delta_rw_w,3);
    free_variable(w_torque,3);
    free_variable(rw_w_new,3);

    return rw_state;
}

// int main(){

//     double** rw_state;
//     double** control_torque = getzeromatrix(3,1);
//     double** wheelspeed = getzeromatrix(3,1);

//     control_torque[0][0] = 0.1;
//     control_torque[1][0] = 0.1;
//     control_torque[2][0] = 0.1;
    
//     wheelspeed[0][0] = 3500*2*pi/60;
//     wheelspeed[1][0] = 500*2*pi/60;
//     wheelspeed[2][0] = 350*2*pi/60;

//     rw_state = reaction_wheels(control_torque, 4.5e-6, 5000*2*pi/60, 150*2*pi/60, wheelspeed);

//     cout<<rw_state[0][0]<<" "<<rw_state[1][0]<<" "<<rw_state[2][0]<<endl;

//     return 0;

// }