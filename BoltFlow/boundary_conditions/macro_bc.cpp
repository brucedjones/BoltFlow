#ifndef MACRO_BC
#define MACRO_BC

// Necessary includes
#include "../macros.cpp"
#include "macro_bc.h"

// These files are only included to remove squiggly red lines in VS2010
#include "../data_types.h"

const macro_condition macro_conditions[15] = { presc_ux,		 presc_uy,			presc_uz,
																presc_uy_ux,	 presc_uz_ux,
																presc_uz_uy,	 presc_ux_uy_uz,
																presc_rho,		 presc_rho_ux,		presc_rho_uy,		presc_rho_uz, 
																presc_rho_ux_uy, presc_rho_ux_uz,	presc_rho_uy_uz,
																presc_rho_ux_uy_uz};

void presc_ux(Node *current_node, Domain *domain)
{
	current_node->u[0] = domain->u[0][current_node->ixd];
}

void presc_uy(Node *current_node, Domain *domain)
{
	current_node->u[1] = domain->u[1][current_node->ixd];
}

void presc_uz(Node *current_node, Domain *domain)
{
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_uy_ux(Node *current_node, Domain *domain)
{
	current_node->u[0] = domain->u[0][current_node->ixd];
	current_node->u[1] = domain->u[1][current_node->ixd];
}

void presc_uz_ux(Node *current_node, Domain *domain)
{
	current_node->u[0] = domain->u[0][current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_uz_uy(Node *current_node, Domain *domain)
{
	current_node->u[1] = domain->u[1][current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_ux_uy_uz(Node *current_node, Domain *domain)
{
	current_node->u[0] = domain->u[0][current_node->ixd];
	current_node->u[1] = domain->u[1][current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_rho(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
}

void presc_rho_ux(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	current_node->u[0] = domain->u[0][current_node->ixd];
}

void presc_rho_uy(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	current_node->u[1] = domain->u[1][current_node->ixd];
}

void presc_rho_uz(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_rho_ux_uy(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	current_node->u[0] = domain->u[0][current_node->ixd];
	current_node->u[1] = domain->u[1][current_node->ixd];
}

void presc_rho_ux_uz(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	current_node->u[0] = domain->u[0][current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_rho_uy_uz(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	current_node->u[1] = domain->u[1][current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}

void presc_rho_ux_uy_uz(Node *current_node, Domain *domain)
{
	current_node->rho = domain->rho[current_node->ixd];
	current_node->u[0] = domain->u[0][current_node->ixd];
	current_node->u[1] = domain->u[1][current_node->ixd];
	#if DIM > 2
		current_node->u[2] = domain->u[2][current_node->ixd];
	#endif
}
#endif
