#ifndef D3Q15_SF_DEFS_H
#define D3Q15_SF_DEFS_H

//#include "d3q15_sf_defs.cpp"

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
void sf_x(Node *current_node, Lattice *lattice);
void sf_X(Node *current_node, Lattice *lattice);
void sf_y(Node *current_node, Lattice *lattice);
void sf_Y(Node *current_node, Lattice *lattice);
void sf_z(Node *current_node, Lattice *lattice);
void sf_Z(Node *current_node, Lattice *lattice);

#endif