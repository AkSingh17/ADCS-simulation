// #ifndef PROGRAMS_NEW
// #define PROGRAMS_NEW
double **getzeromatrix(int,int);
double **MatrixAdd(double **,double **,int,int);
double **MatrixSubtract(double **,double **,int,int);
double **MatrixScalarMultiply(double **,double,int,int);
double **GetIdentityMatrix(int);
double DotProduct(double **,double **);
double **CrossProduct(double **,double **);
double norm(double **);
double **matrixmultiply(double **, double **, int, int, int);
double **TransposeMatrix(double **,int , int );
double TraceofMatrix(double **, int );
double **AdjointofMatrix(double** ,int ,int );
double DeterminantofMatrix( double** );
double **InverseofMatrix(double **, int , int );
double **transformationMatrix(double **);
double** trans_to_quaternion(double**);
double** quatmultiply(double**, double**);
double quatsquare(double **);
double **quatconj(double **);
double **quatinv(double **);
double quatnorm(double**);
double **Unit_vector( double **, int , int );
void free_variable(double ** , int );
double **angle2quat(double**);
double** step_differentiate(double, double**, double**);
double **compute_q_ri(double**, double**);
// #endif