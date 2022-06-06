#include "programs_new.h"
#include "frame_conversion_new.h"
#include <bits/stdc++.h>

/********************Disturbance Torques***********************/
/*
 * Function:  drag_torque
 * --------------------
 *   Calculates the drag torque acting on the satellite
 *   returns: Ndrag (VECTOR)
 */
double  **drag_torque(double **keplerian_elements, double **rnew, double **vnew, double **q_br){

    double argumentOfPerigee = keplerian_elements[4][0];
    double TAnrp = keplerian_elements[5][0];
    double omega =  keplerian_elements[3][0];
    double angleOfIncidence = keplerian_elements[2][0];
    int i;
    /* converting vectors to meters*/
    /*
    Cd=coefficient of drag
    row=atmospheric density
    oe=rate of rotation of earth
    */
    double cd = 2;
    double row = 9.59e-13;
    double oe = 7.292115856e-5;

    /* Satellites dimensions in meters */
    double x = 227e-3;
    double y = 100e-3;
    double z = 100e-3;

    /* Calculating velocity of atmosphere wrt satellite*/
    //Source= fundamentals of astrodynamics chapter 9
    double **vat_eci = getzeromatrix(3,1); 
    
    vat_eci[0][0] = vnew[0][0] + oe * rnew[1][0];
    vat_eci[1][0] = vnew[1][0] - oe * rnew[0][0];
    vat_eci[2][0] = vnew[2][0];

    /* vat in orbit frame */
    double **trans, **vat_orb ; 
    trans = eci_orbit(keplerian_elements);
    vat_orb = matrixmultiply(trans, vat_eci,3,3,1);
    
    /* vat in body frame */
    double **vat_body = getzeromatrix(3,1); 
    double **temp_00 =  getzeromatrix(4,1); 
    temp_00[0][0] = 0;
    for (i = 0; i < 3; i++)
        temp_00[i + 1][0] = vat_orb[i][0];
    double ** q_rb, **q_br_temp;
    q_rb = quatinv(q_br);
    q_br_temp = quatmultiply(q_br,temp_00);
    double **temp_01 = quatmultiply(q_br_temp,q_rb);

    for (i = 0; i < 3; i++)
        vat_body[i][0] = temp_01[i + 1][0];

    double **vat_body_unit;
    vat_body_unit = Unit_vector(vat_body,3,1);

    double dotv = vat_body_unit[0][0];
    double areav = dotv * x * z;
    if (areav < 0)
        areav = -areav;

    double dotc = vat_body_unit[2][0];
    double areac = dotc * y * z;
    if (areac < 0)
        areac = -areac;

    double dotg = vat_body_unit[1][0];
    double areag = dotg * y * x;
    if (areag < 0)
        areag = -areag;

    double area;
    area = areav + areac + areag;

    double **FDrag;
    double tx = -0.5 * cd * row * area * norm(vat_body) * norm(vat_body);

    FDrag = MatrixScalarMultiply(vat_body_unit,tx,3,1);

    double **temp = getzeromatrix(3,1);
    temp[0][0] = 2.8817e-3;
    temp[1][0] = -3.3e-3;
    temp[2][0] = 1.3601e-3;

    double **NDrag;
    NDrag = CrossProduct(temp, FDrag);
    double nomm = norm(NDrag);

    free_variable(trans,3);
    free_variable(vat_eci,3);
    free_variable(vat_orb,3);
    free_variable(vat_body,3);
    free_variable(temp_00,4);
    free_variable(temp,3);
    free_variable(vat_body_unit,3);
    free_variable(FDrag,3);
    free_variable(q_rb,4);
    free_variable(q_br_temp,4);
    free_variable(temp_01,4);

    return NDrag;
}

/*
 * Function:  gravity_torque
 * --------------------
 *   Calculates the gravity torque acting on the satellite
 *   returns: ngg (VECTOR)
 */

double **gravity_torque(double **keplerian_elements, double **rn_mat,double **q_br){

    double argumentOfPerigee = keplerian_elements[4][0];
    double TAnrp = keplerian_elements[5][0];
    double omega =  keplerian_elements[3][0];
    double angleOfIncidence = keplerian_elements[2][0];

    int i;
    double ixx = 2.33869e-12;
    double iyy = 6.6726503e-12;
    double izz = 5.938327e-12;

    double **mk = getzeromatrix(3,3);
    mk[0][0] = ixx;
    mk[1][1] = iyy;
    mk[2][2] = izz;

    /* r in orbit frame */
    double **trans, **rorbit;
    trans = eci_orbit(keplerian_elements);
    rorbit = matrixmultiply(trans, rn_mat,3,3,1);
    
    /* r in body frame */
    double **r_body = getzeromatrix(3,1);
    double **temp_00 = getzeromatrix(4,1);
    temp_00[0][0] = 0;
    for (i = 0; i < 3; i++)
        temp_00[i + 1][0] = rorbit[i][0];

    double ** q_rb, **q_br_temp;
    q_rb = quatinv(q_br);
    q_br_temp = quatmultiply(q_br,temp_00);
    double **temp_01 = quatmultiply(q_br_temp,q_rb);

    for (i = 0; i < 3; i++)
        r_body[i][0] = temp_00[i + 1][0];

    double **t1, **t_2;
    t1 =  matrixmultiply(mk, r_body,3,3,1);
    t_2 = CrossProduct(r_body, t1);

    /* Gravitational Constant */
    double muu = 398603e+9;

    double d = norm(r_body);
    double dd = d * d * d * d * d;
    double d1 = 3 * muu / dd;
    double **ngg;
    ngg =  MatrixScalarMultiply(t_2, d1,3,1);

    free_variable(mk,3);
    free_variable(trans,3);
    free_variable(rorbit,3);
    free_variable(r_body,3);
    free_variable(temp_00,4);
    free_variable(temp_01,4);
    free_variable(t1,3);
    free_variable(t_2,3);
    free_variable(q_rb,4);
    free_variable(q_br_temp,4);
    return ngg;
}

double **solar_torque(double **keplerian_elements, double **rss_body){

    double argumentOfPerigee = keplerian_elements[4][0];
    double TAnrp = keplerian_elements[5][0];
    double omega =  keplerian_elements[3][0];
    double angleOfIncidence = keplerian_elements[2][0];

    int i;
    double eclipse = 1;
    double x = 227e-3;
    double y = 100e-3;
    double z = 100e-3;
    double x_area = x * z;
    double y_area = y * x;
    double z_area = y * z;

    /*for worst case senario crk is 1(black body),unit given is wrong in bate and muller..see the correct unit from AAUSAT 3  */
    /*%reference for Prmf page 645 bate and muller*/

    double crk = 1;
    double prmf = 4.5565e-6;

    double **unit_rss_body;
    unit_rss_body =  Unit_vector(rss_body,3,1);
    double x1 = unit_rss_body[0][0] * eclipse;
    if (x1 < 0)
        x1 = 0;

    double x2 = -unit_rss_body[0][0] * eclipse;
    if (x2 < 0)
        x2 = 0;

    double y1 = unit_rss_body[1][0] * eclipse;
    if (y1 < 0)
        y1 = 0;

    double y2 = -unit_rss_body[1][0] * eclipse;
    if (y2 < 0)
        y2 = 0;

    double z1 = unit_rss_body[2][0] * eclipse;
    if (z1 < 0)
        z1 = 0;

    double z2 = -unit_rss_body[2][0] * eclipse;
    if (z2 < 0)
        z2 = 0;

    double x1_area = x1 * (x_area);
    double x2_area = x2 * (x_area);
    double z1_area = z1 * (y_area);
    double z2_area = z2 * (y_area);
    double y1_area = y1 * (z_area);
    double y2_area = y2 * (z_area);

    double area = x1_area + x2_area + y1_area + y2_area + z1_area + z2_area;
    double dum = crk * area * prmf;
    double **fsolar;
    fsolar =  MatrixScalarMultiply(unit_rss_body,dum,3,1);

    double **mmm = getzeromatrix(3,1);
    mmm[0][0] = 2.8817e-3;
    mmm[1][0] = -3.3e-3;
    mmm[2][0] = 1.3601e-3;

    double **nsolar;
    nsolar =  CrossProduct(mmm, fsolar);

    free_variable(unit_rss_body,3);
    free_variable(fsolar,3);
    free_variable(mmm,3);

    return nsolar;
}

double **calc_disturbance_torque(double **keplerian_elements, double **state_sim, double **q_br, double **rss_body){

    /// Disturbance Torques
    double argumentOfPerigee = keplerian_elements[4][0];
    double TAnrp = keplerian_elements[5][0];
    double omega =  keplerian_elements[3][0];
    double angleOfIncidence = keplerian_elements[2][0];
    double **rn,**vn;
    
    rn=getzeromatrix(3,1);
    vn=getzeromatrix(3,1);

    rn[0][0]=state_sim[0][0];
    rn[1][0]=state_sim[1][0];
    rn[2][0]=state_sim[2][0];

    vn[0][0]=state_sim[3][0];
    vn[1][0]=state_sim[4][0];
    vn[2][0]=state_sim[5][0];

    double **aerodrag, **gravity, **solar, **Td;

    
    gravity = gravity_torque(keplerian_elements, rn, q_br);
    solar = solar_torque(keplerian_elements,rss_body);

    for(int i=0;i<3;i++){
        rn[i][0] = 1000*rn[i][0];
        vn[i][0] = 1000*vn[i][0];
    }
    aerodrag =  drag_torque(keplerian_elements, rn, vn, q_br);
    double **aero_gravity;
    aero_gravity = MatrixAdd(aerodrag, gravity,3,1);
    Td =  MatrixAdd(aero_gravity, solar,3,1);

    free_variable(aerodrag,3);
    free_variable(gravity,3);
    free_variable(solar,3);
    free_variable(aero_gravity,3);
    free_variable(rn,3);
    free_variable(vn,3);
    return Td;
}