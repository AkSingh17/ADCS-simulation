#include "programs_new.h"
#include <bits/stdc++.h>

double **quest(double a1, double a2, double **r1, double **r2, double **b1, double **b2){
   
    double **r1_transpose, **r2_transpose, **b1_r1transpose, **b2_r2transpose, **a1_b1_r1transpose, **a2_b2_r2transpose;
    double **B, **B_transpose, **S;
    double **b1_cross_r1, **b2_cross_r2, **a1_b1_cross_r1, **a2_b2_cross_r2, **Z;
    double **r1_cross_r2, **b1_cross_b2, **adj_S;
    double **identity_3x3, **identity_alpha, **S_beta, **S_square, **X,  **X_cal1, **X_cal2;
    double **q_optimal, **q_optimal_final;
    double sigma, eigen_max, cos_Tv_minus_Tw, kappa, delta, alpha, beta, gamma, q_opt_cal;
    q_optimal = getzeromatrix(4,1);
   
    //Calculating Attitude Profile Matrix
    r1_transpose = TransposeMatrix(r1,3,1);
    r2_transpose = TransposeMatrix(r2,3,1);
    b1_r1transpose= matrixmultiply(b1,r1_transpose, 1, 3, 3);
    b2_r2transpose= matrixmultiply(b2,r2_transpose, 1, 3, 3);
    a1_b1_r1transpose = MatrixScalarMultiply(b1_r1transpose, a1, 3, 3);
    a2_b2_r2transpose = MatrixScalarMultiply(b2_r2transpose, a2, 3, 3);
    
    B = MatrixAdd(a1_b1_r1transpose, a2_b2_r2transpose, 3, 3);

    //S Matrix
    B_transpose = TransposeMatrix(B,3,3);
    S = MatrixAdd(B, B_transpose, 3, 3);

    //Z Matrix
    b1_cross_r1 = CrossProduct(b1,r1);
    b2_cross_r2 = CrossProduct(b2,r2);
    a1_b1_cross_r1 = MatrixScalarMultiply(b1_cross_r1, a1, 3, 1);
    a2_b2_cross_r2 = MatrixScalarMultiply(b2_cross_r2, a2, 3, 1);
    Z = MatrixAdd(a1_b1_cross_r1, a2_b2_cross_r2, 3, 1);

    //Calculating maximum eigen value (reference for formula: https://sci-hub.mksa.top/10.2514/3.19717)
    r1_cross_r2 = CrossProduct(r1,r2);
    b1_cross_b2 = CrossProduct(b1,b2);
    cos_Tv_minus_Tw = ((DotProduct(r1,r2))*(DotProduct(b1,b2))) + (norm(r1_cross_r2)*norm(b1_cross_b2));   
    eigen_max = sqrt(pow(a1,2) + (2*a1*a2*cos_Tv_minus_Tw) + pow(a2,2));
    
    //Calculate terms alpha, beta, gamma to calculate the optimal quaternion 
    adj_S = AdjointofMatrix(S,3,3);
    sigma = 0.5*TraceofMatrix(S,3);
    kappa = TraceofMatrix(adj_S,3);
    delta = DeterminantofMatrix(S);

    alpha = pow(eigen_max,2) - pow(sigma,2) + kappa;
    beta = eigen_max - sigma;
    gamma = (eigen_max + sigma)*alpha - delta;

    //Calculating Vector X for optimal quaternion
    identity_3x3 = GetIdentityMatrix(3); 
    identity_alpha = MatrixScalarMultiply(identity_3x3, alpha, 3, 3);
    S_beta = MatrixScalarMultiply(S, beta, 3, 3);
    S_square = matrixmultiply(S, S, 3, 3, 3);
    X_cal1 =  MatrixAdd(identity_alpha, S_beta, 3, 3);
    X_cal2 =  MatrixAdd(X_cal1, S_square, 3, 3);
    X = matrixmultiply(X_cal2, Z, 3, 3, 1);

    //Calculating Optimal Quaternion
    q_opt_cal = (sqrt((pow(gamma,2)) + (pow(norm(X),2))));
    
    q_optimal[0][0] = gamma;

    for(int i=0; i<3; i++){
        q_optimal[i+1][0] = X[i][0];
    }
    
    q_optimal_final = MatrixScalarMultiply(q_optimal, 1.0/q_opt_cal, 4, 1);

    // free double pointer variables
    free_variable(r1_transpose,1);
    free_variable(r2_transpose,1);
    free_variable(b1_r1transpose,3);
    free_variable(b2_r2transpose,3);
    free_variable(a1_b1_r1transpose,3);
    free_variable(a2_b2_r2transpose,3);
    free_variable(B,3);
    free_variable(B_transpose,3);
    free_variable(S,3);
    free_variable(b1_cross_r1,3);
    free_variable(b2_cross_r2,3);
    free_variable(a1_b1_cross_r1,3);
    free_variable(a2_b2_cross_r2,3);
    free_variable(Z,3);
    free_variable(r1_cross_r2,3);
    free_variable(b1_cross_b2,3);
    free_variable(adj_S,3);
    free_variable(identity_3x3,3);
    free_variable(identity_alpha,3);
    free_variable(S_beta,3);
    free_variable(S_square,3);
    free_variable(X_cal1,3);
    free_variable(X_cal2,3);
    free_variable(X,3);
    free_variable(q_optimal,4);

    return q_optimal_final;
}