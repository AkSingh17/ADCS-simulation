#include "programs_new.h"

double** tc_pid(double** moi_satellite, double** quaternion_error, double** w_error, double** quaternion_error_integral){
    
    double Kpx, Kpy, Kpz, Kdx, Kdy, Kdz, Kix, Kiy, Kiz;
    double** control_torque_pid = getzeromatrix(3,1);
    double** tempq = getzeromatrix(4,1);

    double step = 1.0;
    
    //pid values without failure cases
    Kpx = Kpy = Kpz = 0.030722;
    Kix = Kiy = Kiz = 0;
    Kdx = Kdy = Kdz = 0.27074;
    
    //Integral of q_error 
    for (int i = 1; i < 4; i++)
        tempq[i][0] = quaternion_error[i][0] * step;
    
    double **quaternion_error_integral_new = MatrixAdd(quaternion_error_integral, tempq, 4, 1);

    control_torque_pid[0][0] = moi_satellite[0][0]*(Kpx * quaternion_error[1][0] + Kdx  *  w_error[0][0] + Kix * quaternion_error_integral_new[1][0]);
    control_torque_pid[1][0] = moi_satellite[1][1]*(Kpy * quaternion_error[2][0] + Kdy  *  w_error[1][0] + Kiy * quaternion_error_integral_new[2][0]);
    control_torque_pid[2][0] = moi_satellite[2][2]*(Kpz * quaternion_error[3][0] + Kdz  *  w_error[2][0] + Kiz * quaternion_error_integral_new[3][0]);

    free_variable(tempq,4);
    free_variable(quaternion_error_integral_new,4);
    return control_torque_pid;
    
}