#include "programs_new.cpp"
#include <bits/stdc++.h>

double **quest(double a1, double a2, double **r1, double **r2, double **b1, double **b2){
   
    double **Z = getzeromatrix(3,1);
    double **q_optimal = getzeromatrix(4,1);
    double **B, **b1_cross_r1transpose, **b2_cross_r2transpose, **S, **X, **adj_S, **X_cal1, **X_cal2;
    double sigma, eigen_max, cos_Tv_minus_Tw, alpha, beta, gamma, q_opt_cal;
    
    //Calculating Attitude Profile Matrix
    b1_cross_r1transpose= matrixmultiply(b1,TransposeMatrix(r1,3,1), 1, 3, 3);
    b2_cross_r2transpose= matrixmultiply(b2,TransposeMatrix(r2,3,1), 1, 3, 3);
   
    // for(int i=0; i<3; i++){
    //     for(int j=0; j<3; j++){
    //         cout<< b1_cross_r1transpose[i][j]<<"\t";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    
    B = MatrixAdd(MatrixScalarMultiply(b1_cross_r1transpose, a1, 3, 3), MatrixScalarMultiply(b2_cross_r2transpose, a2, 3, 3), 3, 3);

    // for(int i=0; i<3; i++){
    //     for(int j=0; j<3; j++){
    //         cout<< B[i][j]<<"\t";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;

    // S Matrix
    S = MatrixAdd(B, TransposeMatrix(B,3,3), 3, 3);

    // for(int i=0; i<3; i++){
    //     for(int j=0; j<3; j++){
    //         cout<< S[i][j]<<"\t";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;

    //Z Matrix
    // Z[0][0] = B[1][2] - B[2][1];
    // Z[1][0] = B[2][0] - B[0][2];
    // Z[2][0] = B[0][1] - B[1][0];

    double** Z_temp1 = CrossProduct(b1,r1);
    double** Z_temp2 = CrossProduct(b2,r2);
    Z = MatrixAdd(MatrixScalarMultiply(Z_temp1, a1, 3, 1), MatrixScalarMultiply(Z_temp2, a2, 3, 1), 3, 1);

    //Calculating maximum eigen value (reference for formula: https://sci-hub.mksa.top/10.2514/3.19717)
    cos_Tv_minus_Tw = ((DotProduct(r1,r2))*(DotProduct(b1,b2))) + (norm(CrossProduct(r1,r2))*norm(CrossProduct(b1,b2)));   
    eigen_max = sqrt(pow(a1,2) + (2*a1*a2*cos_Tv_minus_Tw) + pow(a2,2));
    // cout<<eigen_max<<endl;
    
    //Calculate terms alpha, beta, gamma to calculate the optimal quaternion 
    adj_S = AdjointofMatrix(S,3,3);

    //   for(int i=0; i<3; i++){
    //     for(int j=0; j<3; j++){
    //         cout<< adj_S[i][j]<<"\t";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;

    // cout<<TraceofMatrix(adj_S,3)<<endl;

    sigma = 0.5*TraceofMatrix(S,3);
    double kappa = TraceofMatrix(adj_S,3);
    double delta = DeterminantofMatrix(S);

    alpha = pow(eigen_max,2) - pow(sigma,2) + kappa;
    // cout<<alpha<<endl;
    beta = eigen_max - sigma;
    // cout<<beta<<endl;
    gamma = (eigen_max + sigma)*alpha - delta;
    // cout<<gamma<<endl;

    //Calculating Vector X for optimal quaternion
    X_cal1 =  MatrixAdd((MatrixScalarMultiply(GetIdentityMatrix(3), alpha, 3, 3)), (MatrixScalarMultiply(S, beta, 3, 3)), 3, 3);
    X_cal2 =  MatrixAdd(X_cal1, (matrixmultiply(S, S, 3, 3, 3)), 3, 3);

    //    for(int i=0; i<3; i++){
    //        for(int j=0; j<3; j++){
    //            cout<< X_cal2[i][j]<<"\t";
    //        }
    //        cout<<endl;
    //    }
    //    cout<<endl;

    X = matrixmultiply(X_cal2, Z, 3, 3, 1);

    // cout<<X[0][0]<<" "<<X[1][0]<<" "<<X[2][0]<<endl;

    //Calculating Optimal Quaternion
    q_opt_cal = (sqrt((pow(gamma,2)) + (pow(norm(X),2))));

    // cout<<q_opt_cal<<endl;
    
    q_optimal[0][0] = gamma;

    for(int i=0; i<3; i++){
        q_optimal[i+1][0] = X[i][0];
    }
    
    q_optimal = MatrixScalarMultiply(q_optimal, 1/q_opt_cal, 4, 1);

    //free double pointer variables
    free_variable(B,3);
    free_variable(b1_cross_r1transpose,3);
    free_variable(b2_cross_r2transpose,3);
    free_variable(S,3);
    free_variable(Z,3);
    free_variable(X,3);
    free_variable(adj_S,3);
    free_variable(X_cal1,3);
    free_variable(X_cal2,3);
    
    return q_optimal;
}

// double** attitude_profile_matrix(double a1,double a2,double** b1,double** b2,double** r1,double** r2)
// {
//     double** B1=MatrixScalarMultiply(b1,a1,3,1);
//     double **B2=MatrixScalarMultiply(b2,a2,3,1);

//     double** B=MatrixAdd(matrixmultiply(B1,TransposeMatrix(r1,3,1),1,3,3),matrixmultiply(B2,TransposeMatrix(r2,3,1),1,3,3),3,3);

//     return B;
// }


// double max_eigen_value (double a1,double a2, double** b1, double** b2, double** r1, double** r2)
// {
//     double eig_max=sqrt(pow(a1,2)+pow(a2,2)+2*a1*a2*(DotProduct(b1,b2)*DotProduct(r1,r2)+norm(CrossProduct(b1,b2))*norm(CrossProduct(r1,r2))));
//     return eig_max;
// }


// double** quaternion_optimum(double** B, double eig_max)
// {
//     double** q;
//     double ** q_m1= getzeromatrix(6,1);

//     q = MatrixAdd(B, TransposeMatrix(B,3,3), 3, 3);

//     q_m1[0][0] = B[1][2]-B[2][1];
//     q_m1[1][0]=B[2][0]-B[0][2];
//     q_m1[2][0]=B[0][1]-B[1][0];

//     double** adj=AdjointofMatrix(q,3,3);
//     double alpha=powf(eig_max,2)-powf(TraceofMatrix(B,3),2)+TraceofMatrix(adj,3);
//     double beta=eig_max-TraceofMatrix(B,3);
//     double gamma=alpha*(eig_max+TraceofMatrix(B,3))- DeterminantofMatrix(q);

//     double** temp = MatrixAdd(MatrixScalarMultiply(GetIdentityMatrix(3),alpha,3,3),MatrixScalarMultiply(q,beta,3,3),3,3);
//     double** X = MatrixAdd(matrixmultiply(q,q,3,3,3),temp,3,3);
//     X = matrixmultiply(X,q_m1,3,3,1);

//     double** q_opt = getzeromatrix(4,1);
//     int i;
//     for(i=0; i<3; i++)
//         q_opt[i][0] = X[i][0]; 

//     double temp1 = pow(gamma,2) + pow(norm(X),2);
//     temp1 = sqrt(temp1);
//     temp1 = 1/temp1;
//     q_opt[3][0] = gamma;
//     q_opt = MatrixScalarMultiply(q_opt,temp1,4,1);
  
//     return q_opt;
// }


// double** quest(double a1, double a2, double** r1, double** r2, double** b1, double** b2)
// {
//     double** B = attitude_profile_matrix(a1, a2, b1, b2, r1, r2);
 
//     double eig_max = max_eigen_value(a1, a2, b1, b2, r1, r2);
   
//     double** q_opt = quaternion_optimum(B, eig_max);

//     double temp = q_opt[0][0];
//     q_opt[0][0] = q_opt[3][0];
//     q_opt[3][0] = q_opt[2][0];
//     q_opt[2][0] = q_opt[1][0];
//     q_opt[1][0] = temp;

//     return q_opt;
// }


int main(){

    double** b1_quest= getzeromatrix(3,1);
    double** b2_quest= getzeromatrix(3,1);
    double** r1_quest= getzeromatrix(3,1);
    double** r2_quest= getzeromatrix(3,1);
    double** temp_quat1= getzeromatrix(4,1);
    double** temp_quat2= getzeromatrix(4,1);
    double** q_br= getzeromatrix(4,1);

    // 1.95689e+07 -1.46326e+08 -1.92104e+07
    r1_quest[0][0] = 1.95689e+07 ; 
    r1_quest[1][0] = -1.46326e+08 ;
    r1_quest[2][0] = -1.92104e+07;

    // 3.34862e-06 -2.15738e-07 7.21059e-06
    r2_quest[0][0] =  3.34862e-06;
    r2_quest[1][0] = -2.15738e-07;
    r2_quest[2][0] = 7.21059e-06;

    // 1.95689e+07 -1.46326e+08 -1.92104e+07
    b1_quest[0][0] = 1.95689e+07 ; 
    b1_quest[1][0] = -1.46326e+08 ;
    b1_quest[2][0] = -1.92104e+07;

    // 3.34862e-06 -2.15738e-07 7.21059e-06
    b2_quest[0][0] =  3.34862e-06;
    b2_quest[1][0] = -2.15738e-07;
    b2_quest[2][0] = 7.21059e-06;

    //  for(int i=0; i<3; i++){
    //     temp_quat1[i+1][0] = r1_quest[i][0];
    // }
 
    // temp_quat1 = quatmultiply(quatmultiply(q_br,temp_quat1),quatinv(q_br));
        
    // for(int i=0;i<3;i++){
    //     b1_quest[i][0] = temp_quat1[i+1][0];
    // }

    // for(int i=0; i<3; i++){
    //     temp_quat2[i+1][0] = r2_quest[i][0];
    // }
 
    // temp_quat2 = quatmultiply(quatmultiply(q_br,temp_quat2),quatinv(q_br));
        
    // for(int i=0;i<3;i++){
    //     b2_quest[i][0] = temp_quat2[i+1][0];
    // }

    q_br = quatinv(quest(0.5,0.5,Unit_vector(r1_quest,3,1),Unit_vector(r2_quest,3,1),Unit_vector(b1_quest,3,1),Unit_vector(b2_quest,3,1)));

    // cout<<b1_quest[0][0]<<" "<<b1_quest[1][0]<<" "<<b1_quest[2][0]<<endl;
    // cout<<b2_quest[0][0]<<" "<<b2_quest[1][0]<<" "<<b2_quest[2][0]<<endl;
    cout<<q_br[0][0]<<" "<<q_br[1][0]<<" "<<q_br[2][0]<<" "<<q_br[3][0]<<endl;

    return 0;
}





