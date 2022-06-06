#include <bits/stdc++.h>
using namespace std;

double **getzeromatrix(int m,int n){

    double **temp;

    temp = (double**)malloc(m*sizeof(double*));

    for(int i=0;i<m;i++) {
        temp[i] = (double*)malloc(n*sizeof(double));
    }

    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            temp[i][j] = 0;
        }
    }
    return temp;
}

void free_variable(double ** variable, int m){

    for( int i=0; i<m; i++){
        
        free(variable[i]);
        variable[i] = NULL;
    }
    
    free(variable);
    variable = NULL;
}

double norm(double **v){

    double norm = 0.0;
    norm = (v[0][0] * v[0][0]) + (v[1][0] * v[1][0]) + (v[2][0] * v[2][0]);
    norm = pow(norm,0.5);
    return norm;
}

double **MatrixAdd(double **A,double **B,int m,int n){

    double **temp = getzeromatrix(m,n);

    for(int i=0;i<m;i++) {
        for(int j=0;j<n;j++) {
            temp[i][j] = A[i][j] + B[i][j];
        }
    }
    return temp;
}

double **Unit_vector( double **M, int m, int n){

    double **A = getzeromatrix(m,n);
    int i, j;
    double unit_norm = norm(M);

    for(i=0; i<3; i++){
        A[i][0]= M[i][0]/unit_norm;
    }  
    return A;
}

double **MatrixSubtract(double **A,double **B,int m,int n){

    double **temp = getzeromatrix(m,n);

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            temp[i][j] = A[i][j] - B[i][j];
        }
    }
    return temp;
}

double **MatrixScalarMultiply(double **A,double k,int m,int n){

    double **temp = getzeromatrix(m,n);

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            temp[i][j] = A[i][j] * k;
        }
    }
    return temp;
}

double **GetIdentityMatrix(int n){

    double **temp = getzeromatrix(n,n);

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i == j) 
                temp[i][j] = 1;
        }
    }
    return temp;
}

double DotProduct(double **A,double **B){

    double dotPdt = 0.0;
    dotPdt = (A[0][0]* B[0][0]) + (A[1][0]* B[1][0]) + (A[2][0] * B[2][0]);
    return dotPdt;
}

double **CrossProduct(double **u,double **v){

    double **a = getzeromatrix(3,1);

    a[0][0] = (u[1][0] * v[2][0]) - (u[2][0] * v[1][0]);
    a[1][0] = -((u[0][0] * v[2][0]) - (u[2][0] * v[0][0]));
    a[2][0] = (u[0][0] * v[1][0]) - (u[1][0] * v[0][0]);

    return a;
}

double **matrixmultiply(double **A, double **B, int A_column, int m, int n){

    double **C = getzeromatrix(m,n);

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            C[i][j] = 0;
            for(int k=0; k<A_column; k++){
                C[i][j]+=A[i][k]*B[k][j];
            }
        }
    }

    return C;
}
double **TransposeMatrix(double **A,int m, int n){

    double **A_transpose = getzeromatrix(n,m);
    int i, j;

    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            A_transpose[i][j]=A[j][i];
        }
    }

    return A_transpose;
}

double TraceofMatrix(double **A, int n){

    double trace=0;
    int i;

    for(i=0; i<n; i++){
        trace= trace+A[i][i];
    }

    return trace;
}

double **AdjointofMatrix(double** M,int m,int n){
    
    double **A = getzeromatrix(m,n);

    A[0][0] = (M[1][1] * M[2][2]) - (M[1][2] * M[2][1]);
    A[0][1] = (M[0][2] * M[2][1]) - (M[2][2] * M[0][1]);
    A[0][2] = (M[1][1] * M[1][2]) - (M[1][1] * M[0][2]);
    A[1][0] = (M[2][0] * M[1][2]) - (M[2][2] * M[1][0]);
    A[1][1] = (M[2][2] * M[0][0]) - (M[2][0] * M[0][2]);
    A[1][2] = (M[0][2] * M[1][0]) - (M[1][2] * M[0][0]);
    A[2][0] = (M[1][1] * M[2][1]) - (M[1][1] * M[2][0]);
    A[2][1] = (M[2][0] * M[0][1]) - (M[0][0] * M[2][1]);
    A[2][2] = (M[0][0] * M[1][1]) - (M[1][0] * M[0][1]);

    return A;

}

double DeterminantofMatrix( double** a){

    double determinant=0;

    determinant = a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) 
                - a[0][1] * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) 
                + a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);

    return determinant;
}

double **InverseofMatrix(double **M, int m, int n){
    
    double **A = getzeromatrix(m,n);
    double **B = AdjointofMatrix(M,m,n);
    double C = DeterminantofMatrix(M);
    int i,j;

    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            A[i][j] = B[i][j]/C;
        }
    }

    free_variable(B,m);

    return A;
}

double **calc_Tbr(double **x, double **y, double **z){

    double **temp = getzeromatrix(3,3);
    double **temp1 = getzeromatrix(3,3);

    temp = matrixmultiply(x, y, 3, 3, 3);
    temp1 = matrixmultiply(temp, z, 3, 3, 3);

    return temp1;
}
double **transformationMatrix(double **e){

    ////transform from orbital ref frame to body frame
    double **Rx, **Ry, **Rz, **Tbr;
    Rx = getzeromatrix(3,3);
    Ry = getzeromatrix(3,3);
    Rz = getzeromatrix(3,3);

    Rx[0][0] = 1;
    Rx[1][1] = cos(e[0][0]);
    Rx[1][2] = -sin(e[0][0]);
    Rx[2][1] = sin(e[0][0]);
    Rx[2][2] = cos(e[0][0]);

    Ry[0][0] = cos(e[1][0]);
    Ry[0][2] = sin(e[1][0]);
    Ry[1][1] = 1;
    Ry[2][0] = -sin(e[1][0]);
    Ry[2][2] = cos(e[1][0]);

    Rz[0][0] = cos(e[2][0]);
    Rz[0][1] = -sin(e[2][0]);
    Rz[1][0] = sin(e[2][0]);
    Rz[1][1] = cos(e[2][0]);
    Rz[2][2] = 1;

    Tbr = calc_Tbr(Rx, Ry, Rz);

    free_variable(Rx,3);
    free_variable(Ry,3);
    free_variable(Rz,3);

    return Tbr;
}

double  **trans_to_quaternion(double **m){
    
    double **q = getzeromatrix(4,1);

    double m00 = m[0][0];
    double m01 = m[0][1];
    double m02 = m[0][2];
    double m10 = m[1][0];
    double m11 = m[1][1];
    double m12 = m[1][2];
    double m20 = m[2][0];
    double m21 = m[2][1];
    double m22 = m[2][2];

    if(m11>-m22 && m00 > -m11 && m00>-m22) {
        
        q[0][0] = sqrt(1 + m00 + m11 + m22)/2;
        q[1][0] = ((m12 - m21)/sqrt(1 + m00 + m11 + m22))/2;
        q[2][0] = ((m20 - m02)/sqrt(1 + m00 + m11 + m22))/2;
        q[3][0] = ((m01 - m10)/sqrt(1 + m00 + m11 + m22))/2;
    }
    
    else if(m11<-m22 && m00>m11 && m00>m22) {
      
        q[0][0] = ((m12 - m21)/sqrt(1 + m00 - m11 - m22))/2;
        q[1][0] = sqrt(1 + m00 - m11 - m22)/2;
        q[2][0] = ((m01 + m10)/sqrt(1 + m00 - m11 - m22))/2;
        q[3][0] = ((m20 + m02)/sqrt(1 + m00 - m11 - m22))/2;
    }
    
    else if(m11>m22 && m00<m11 && m00<-m22) {

        q[0][0] = ((m20 - m02)/sqrt(1 - m00 + m11 - m22))/2;
        q[1][0] = ((m01 + m10)/sqrt(1 - m00 + m11 - m22))/2;
        q[2][0] = sqrt(1 - m00 + m11 - m22)/2;
        q[3][0] = ((m12 + m21)/sqrt(1 - m00 + m11 - m22))/2;
    }
    
    else{

        q[0][0] = ((m01 - m10)/sqrt(1 - m00 - m11 + m22))/2;
        q[1][0] = ((m20 + m02)/sqrt(1 - m00 - m11 + m22))/2;
        q[2][0] = ((m12 + m21)/sqrt(1 - m00 - m11 + m22))/2;
        q[3][0] = sqrt(1 - m00 - m11 + m22)/2;
    }

    return q;
}


double **quatmultiply(double **q, double **r){
    
    double **t = getzeromatrix(4,1);
    
    t[0][0] = r[0][0]*q[0][0] - r[1][0]*q[1][0] - r[2][0]*q[2][0] - r[3][0]*q[3][0];
    t[1][0] = r[0][0]*q[1][0] + r[1][0]*q[0][0] - r[2][0]*q[3][0] + r[3][0]*q[2][0];
    t[2][0] = r[0][0]*q[2][0] + r[1][0]*q[3][0] + r[2][0]*q[0][0] - r[3][0]*q[1][0];
    t[3][0] = r[0][0]*q[3][0] - r[1][0]*q[2][0] + r[2][0]*q[1][0] + r[3][0]*q[0][0];
    
    return t;
}

double quatsquare(double **q){

    int i;
    double sq = 0;

    for(i=0; i<4; i++) {
        sq = sq + q[i][0]*q[i][0];
    }
    return sq;
}

double **quatconj(double **q){

    int i;
    double **qconj = getzeromatrix(4,1);

    qconj[0][0] = q[0][0];

    for(i=1; i<4; i++){
        qconj[i][0] = q[i][0]*-1;
    }
    return qconj;
}

double **quatinv(double **q){

    double sq;
    int i;
    double **qconj, **qinv;
    qinv = getzeromatrix(4,1);

    sq = quatsquare(q);
    qconj = quatconj(q);

    for(i=0; i<4; i++){
        qinv[i][0] = qconj[i][0]/sq;
    }

    free_variable(qconj,4);

    return qinv;
}

double quatnorm(double** q){
    
    return sqrt(quatsquare(q));
}

double **step_differentiate(double step, double** q_current, double** q_previous){

    double **q_difference, **q_dot, **temp_q, **temp_q1, **omega;
    omega = getzeromatrix(3,1);

    step = 1/step;

    q_difference = MatrixSubtract(q_current,q_previous,4,1);
    q_dot = MatrixScalarMultiply(q_difference,step,4,1);

    double **q_current_inv;

    q_current_inv = quatinv(q_current);

    temp_q = quatmultiply(q_current_inv,q_dot);
    temp_q1 = MatrixScalarMultiply(temp_q,2,4,1);

    for(int i =0; i<3; i++){
        omega[i][0] = temp_q1[i+1][0];
    }

    free_variable(q_difference,4);
    free_variable(q_dot,4);
    free_variable(temp_q,4);
    free_variable(temp_q1,4);
    free_variable(q_current_inv,4);

    return omega;
    
}

double **compute_q_ri(double** trans, double** q_ri_prev){
       
    double **q_ri, **q_ri_final;
    double **trans_quat;

    trans_quat = trans_to_quaternion(trans);

    q_ri  = quatinv(trans_quat);

    //for q_ri smoothening  
    if(((q_ri[0][0]*q_ri_prev[0][0] < 0) && (q_ri[1][0]*q_ri_prev[1][0] < 0) & (q_ri[2][0]*q_ri_prev[2][0] < 0)) || 
       ((q_ri[0][0]*q_ri_prev[0][0] < 0) && (q_ri[2][0]*q_ri_prev[2][0] < 0) & (q_ri[3][0]*q_ri_prev[3][0] < 0)) || 
       ((q_ri[0][0]*q_ri_prev[0][0] < 0) && (q_ri[1][0]*q_ri_prev[1][0] < 0) & (q_ri[3][0]*q_ri_prev[3][0] < 0)) || 
       ((q_ri[1][0]*q_ri_prev[1][0] < 0) && (q_ri[2][0]*q_ri_prev[2][0] < 0) & (q_ri[3][0]*q_ri_prev[3][0] < 0))){
       
        q_ri_final = MatrixScalarMultiply(q_ri,-1,4,1);
    }

    else{
        q_ri_final  = quatinv(trans_quat);
    }

    free_variable(q_ri,4);
    free_variable(trans_quat,4);

    return q_ri_final;
}