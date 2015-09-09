#ifndef COLLISION_H
#define COLLISION_H

#include "data_types.h"

// SUPPORT FUNCTIONS
double u_square(Node *current_node);
double e_mul_u(Node *current_node, DomainConstant, int *i);
void turbulent_viscosity(Node *current_node, DomainConstant *domain_constants, double *f_eq, double *tau);

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
void bgk_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void bgk_guo_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void bgk_ntpor_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void bgk_ntpor_guo_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void mrt_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void mrt_guo_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void mrt_ntpor_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void mrt_ntpor_guo_collision(Node *current_node, DomainConstant *domain_constants, double *tau);
void meq_d2q9(Node *current_node, double *meq);
void meq_d3q15(Node *current_node, double *meq);
void bounceback(Node *current_node, DomainConstant *domain_constants, double *tau);

#endif