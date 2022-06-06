#include <bits/stdc++.h>
#include "programs_new.h"
using namespace std;

#define mu 398600.4418
#define pi 3.14159265359
#define R 6378.137 //radius of the earth
#define J2 0.00108263

double **cartesian_to_keplerian(double **state){

    double **position, **velocity, **angular_momentum, **mean_motion;
    double **unit_position, **v_cross_H, **v_cross_H_by_mu, **eccentricity_vector, **keplerian_elements;  
    double semimajor_axis,eccentricity,inclination,right_ascension_of_ascending_node,argument_of_perigee,true_anomaly;

    position = getzeromatrix(3,1);
    velocity = getzeromatrix(3,1);
    mean_motion = getzeromatrix(3,1);
    keplerian_elements = getzeromatrix(6,1);

    for(int i=0;i<3;i++){
        position[i][0] = state[i][0]; 
        velocity[i][0] = state[i+3][0];
    }
    
    angular_momentum = CrossProduct(position,velocity);

    mean_motion[0][0] = -angular_momentum[1][0];
    mean_motion[1][0] = angular_momentum[0][0];
    mean_motion[2][0] = 0;

    unit_position =  MatrixScalarMultiply(position, 1/norm(position), 3, 1);
    v_cross_H = CrossProduct(velocity,angular_momentum);
    v_cross_H_by_mu = MatrixScalarMultiply(v_cross_H, 1/mu, 3, 1);
    eccentricity_vector = MatrixSubtract(v_cross_H_by_mu, unit_position, 3, 1);

    eccentricity = norm(eccentricity_vector);
    semimajor_axis = -mu*norm(position)/(norm(position)*pow(norm(velocity),2) - 2*mu);
    inclination = acos(angular_momentum[2][0]/norm(angular_momentum));

    /*RAAN considerations*/
    right_ascension_of_ascending_node = acos(mean_motion[0][0]/norm(mean_motion));
    if(mean_motion[1][0] >= 0){
        right_ascension_of_ascending_node = right_ascension_of_ascending_node;
    }
    else{
        right_ascension_of_ascending_node = 2*pi - right_ascension_of_ascending_node;
    }

    /*AGP considerations*/
    argument_of_perigee = acos(DotProduct(mean_motion,eccentricity_vector)/((norm(mean_motion))*(norm(eccentricity_vector))));

    if(eccentricity_vector[2][0] >= 0){
        argument_of_perigee = argument_of_perigee;
    }
    else{
        argument_of_perigee = 2*pi - argument_of_perigee;
    }
    
    if(DotProduct(position,velocity)<0){
        true_anomaly = 2*pi - acos(DotProduct(position,eccentricity_vector)/(norm(eccentricity_vector)*norm(position)));
    }
    else{
        true_anomaly = acos(DotProduct(position,eccentricity_vector)/(norm(eccentricity_vector)*norm(position)));
    }

    keplerian_elements[0][0] = semimajor_axis;
    keplerian_elements[1][0] = eccentricity;
    keplerian_elements[2][0] = inclination;
    keplerian_elements[3][0] = right_ascension_of_ascending_node;
    keplerian_elements[4][0] = argument_of_perigee;
    keplerian_elements[5][0] = true_anomaly;

    free_variable(position,3);
    free_variable(velocity,3);
    free_variable(angular_momentum,3);
    free_variable(mean_motion,3);
    free_variable(unit_position,3);
    free_variable(v_cross_H,3);
    free_variable(v_cross_H_by_mu,3);
    free_variable(eccentricity_vector,3);

    return keplerian_elements;
}

double **J2_accelerations(double **position){

    double **J2_acc = getzeromatrix(3,1);
    double norm_r,J_constant;

    norm_r = norm(position);
    J_constant= (-1.5 * mu * J2 * R * R) / pow(norm_r,5);

    J2_acc[0][0] = J_constant * (1 - (5 * pow((position[2][0] / norm_r), 2))) * position[0][0];
    J2_acc[1][0] = J_constant * (1 - (5 * pow((position[2][0] / norm_r), 2))) * position[1][0];
    J2_acc[2][0] = J_constant * (3 - (5 * pow((position[2][0] / norm_r), 2))) * position[2][0];

    return J2_acc;
}

double **state_derivative(double **state){

    double **position, **j2_acc, **acc_twobody, **acc_total, **cartesian_derivative;    
    double mod_r;

    position = getzeromatrix(3,1);
    acc_twobody = getzeromatrix(3,1); 
    cartesian_derivative = getzeromatrix(6,1); 

    for(int i=0;i<3;i++){
        position[i][0] = state[i][0];
    }
    
    mod_r = norm(state);

    for (int i = 0; i < 3; i++){
        acc_twobody[i][0] = (-mu / pow(mod_r, 3)) * position[i][0];
    }

    j2_acc = J2_accelerations(position);

    acc_total = MatrixAdd(acc_twobody, j2_acc, 3, 1);

    for (int i = 0; i < 3; i++){
        cartesian_derivative[i][0] = state[3 + i][0];
        cartesian_derivative[3 + i][0] = acc_total[i][0];
    }

    free_variable(acc_total,3);
    free_variable(acc_twobody,3);
    free_variable(position,3);
    free_variable(j2_acc,3);

    return cartesian_derivative;
}

double **RK4(double **state, double initial_time, double final_time, double step){

    double time = 0;
    double **state_der, **k1, **k2, **k3, **k4, **k1_div_2, **k2_div_2, **change_state;

    k1 = getzeromatrix(6,1);
    k2 = getzeromatrix(6,1);
    k3 = getzeromatrix(6,1);
    k4 = getzeromatrix(6,1); 

    time = initial_time;

    if(final_time < initial_time)
        step = -step;

    while (fabs(final_time - time) >= 0.000001) {

        if (fabs(final_time - time) < fabs(step))
            step = final_time - time;

        state_der = state_derivative(state);

        for (int i = 0; i < 6; i++) 
            k1[i][0] = step * state_der[i][0];

        k1_div_2 = MatrixScalarMultiply(k1, 0.5, 6 ,1);
        change_state = MatrixAdd(state, k1_div_2, 6, 1);

        free_variable(state_der,6);
        free_variable(k1_div_2,6);
            
        state_der = state_derivative(change_state);
         
        free_variable(change_state,6);

        for (int i = 0; i < 6; i++) 
            k2[i][0] = step * state_der[i][0];
           
        free_variable(state_der,6);

        k2_div_2 = MatrixScalarMultiply(k2, 0.5, 6, 1);        
        change_state = MatrixAdd(state, k2_div_2, 6, 1);

        state_der = state_derivative(change_state);

        free_variable(k2_div_2,6);
        free_variable(change_state,6);

        for (int i = 0; i < 6; i++) 
            k3[i][0] = step * state_der[i][0];

        free_variable(state_der,6);

        change_state = MatrixAdd(state, k3, 6, 1);
        state_der = state_derivative(change_state);

        free_variable(change_state,6);

        for (int i = 0; i < 6; i++) 
            k4[i][0] = step * state_der[i][0];

        free_variable(state_der,6);
        
        for (int i = 0; i < 6; i++) 
            state[i][0] = state[i][0] + (k1[i][0] + 2 * k2[i][0] + 2 * k3[i][0] + k4[i][0]) / 6;
        
        time = time + step;
    }
    
    free_variable(k1,6);
    free_variable(k2,6);
    free_variable(k3,6);
    free_variable(k4,6);

    return state;
}

double **orbit_propagate(double **state, double initial_time, double final_time, double step){
    
    double** state_new;

    state_new = RK4(state,initial_time,final_time,step);

    return state_new;
}