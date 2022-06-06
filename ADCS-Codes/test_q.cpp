#include "programs_new.cpp"

double** funcEval_dqdt(double** y, double** w){
   
    int i;
    double** q;

    q = quatmultiply(y, w);
    q = MatrixScalarMultiply(q, 0.5, 4 ,1);

    // cout<<"dqdt "<<q[0][0]<<" "<<q[1][0]<<" "<<q[2][0]<<" "<<q[3][0]<<endl;


    return q;
}

double **integrate_quaternion(double** q, double** w_bi, double step, double h){
    
    double **k1, **k2, **k3, **k4; 
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
    }
    
    q = MatrixScalarMultiply(q, 1/quatnorm(q), 4, 1);
   
    free_variable(k1,4);
    free_variable(k2,4);
    free_variable(k3,4);
    free_variable(k4,4);
    free_variable(t,4);
    free_variable(w,4);

    return q;
}

int main(){

    double ** q, **w;

    q = getzeromatrix(4,1);
    w = getzeromatrix(3,1);

    q[0][0] = 1.0;
    q[1][0] = 0;
    q[2][0] = 0;
    q[3][0] = 0;

    w[0][0] = 1.0;
    w[1][0] = 0;
    w[2][0] = 0;

    for(int i=0; i<10; i++){
        q = integrate_quaternion(q,w,1.0,0.2);
        cout<<q[0][0]<<" "<<q[1][0]<<" "<<q[2][0]<<" "<<q[3][0]<<endl<<endl;
    }


}