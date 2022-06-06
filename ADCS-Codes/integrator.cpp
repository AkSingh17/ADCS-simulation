#include <bits/stdc++.h>
#include "programs_new.h"
using namespace std;

// include reaction wheel parameters

double** funcEval_dwdt(double** I, double** Idash, double** Tc, double** Td, double** w, double** rw_w){

    int i, j;
    double flywheel_I = 4.5e-6;
    double **wbi, **d, **ang_momentum_satellite, **ang_momentum_rw, **total_ang_momentum;
    double** T =  getzeromatrix(3,1);

    ang_momentum_satellite = matrixmultiply(I, w , 3, 3 ,1);
 
    ang_momentum_rw = MatrixScalarMultiply(rw_w, flywheel_I, 3, 1);

    total_ang_momentum = MatrixAdd(ang_momentum_satellite, ang_momentum_rw, 3, 1);

    d = CrossProduct(w, total_ang_momentum);

    for (i = 0; i < 3; i++)
        T[i][0] = Tc[i][0] + Td[i][0] - d[i][0];

    wbi = matrixmultiply(Idash, T, 3, 3, 1);

    free_variable(ang_momentum_satellite,3);
    free_variable(ang_momentum_rw,3);
    free_variable(total_ang_momentum,3);
    free_variable(d,3);
    free_variable(T,3);

    return wbi;
}

double** funcEval_dqdt(double** y, double** w){
   
    int i;
    double** q;

    q = quatmultiply(y, w);
    double **q_final = MatrixScalarMultiply(q, 0.5, 4 ,1);
    
    free_variable(q,4);
    return q_final;
}

double **integrate_omega(double** I, double** Iinv, double** Tc, double** Td, double** w_bi, double** rw_w, double step, double h){
    float k = 0;
    int i;
    double **a, **b, **c, **d;
    double** temp = getzeromatrix(3,1);

    for (k = 0; k < step; k += h){

        //propagating w_bi
        
        a = funcEval_dwdt(I, Iinv, Tc, Td, w_bi, rw_w);
        for (i = 0; i < 3; i++)
            a[i][0] = a[i][0] * h;

        for (i = 0; i < 3; i++)
            temp[i][0] = w_bi[i][0] + a[i][0] / 2.0;

        b = funcEval_dwdt(I, Iinv, Tc, Td, temp, rw_w);
        for (i = 0; i < 3; i++)
            b[i][0] = b[i][0] * h;

        for (i = 0; i < 3; i++)
            temp[i][0] = w_bi[i][0] + b[i][0] / 2.0;

        c = funcEval_dwdt(I, Iinv, Tc, Td, temp, rw_w);
        for (i = 0; i < 3; i++)
            c[i][0] = c[i][0] * h;

        for (i = 0; i < 3; i++)
            temp[i][0] = w_bi[i][0] + c[i][0];

        d = funcEval_dwdt(I, Iinv, Tc, Td, temp, rw_w);
        for (i = 0; i < 3; i++)
            d[i][0] = d[i][0] * h;

        for (i = 0; i < 3; i++){
            w_bi[i][0] = w_bi[i][0] + ((a[i][0]/ 6.0) + (b[i][0] / 3.0) + (c[i][0] / 3.0) + (d[i][0] / 6.0));
        }

        free_variable(a,3);
        free_variable(b,3);
        free_variable(c,3);
        free_variable(d,3);        
    }  

    free_variable(temp,3);

    return w_bi;
}

double **integrate_quaternion(double** q, double** w_bi, double step, double h){
    
    double **k1, **k2, **k3, **k4, **q_normalised; 
    double** t = getzeromatrix(4,1);
    double** w = getzeromatrix(4,1);
    int j;

    for (float k = 0; k < step; k += h){
      
        for(j=0; j<3; j++){
            w[j+1][0] = w_bi[j][0];
        }

        k1 = funcEval_dqdt(q, w);

        for (j = 0; j < 4; j++)
            t[j][0] = q[j][0] + h * k1[j][0] / 2;

        k2 = funcEval_dqdt(t, w);

        for (j = 0; j < 4; j++)
            t[j][0] = q[j][0] + h * k2[j][0] / 2;

        k3 = funcEval_dqdt(t, w);

        for (j = 0; j < 4; j++)
            t[j][0] = q[j][0] + h * k3[j][0];

        k4 = funcEval_dqdt(t, w);

        for (j = 0; j < 4; j++)
            q[j][0] = q[j][0] + h / 6 * (k1[j][0] + 2 * k2[j][0] + 2 * k3[j][0] + k4[j][0]);
        
        free_variable(k1,4);
        free_variable(k2,4);
        free_variable(k3,4);
        free_variable(k4,4);
    }
    
    q_normalised = MatrixScalarMultiply(q, 1/quatnorm(q), 4, 1);
   
    free_variable(t,4);
    free_variable(w,4);

    return q_normalised;
}