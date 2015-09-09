#ifndef D3Q15_ZH_DEFS_H
#define D3Q15_ZH_DEFS_H

#include "d3q15_zh_defs.cpp"

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
void zh_pressure_x(Node *current_node, Lattice *lattice);
void zh_pressure_X(Node *current_node, Lattice *lattice);
void zh_pressure_y(Node *current_node, Lattice *lattice);
void zh_pressure_Y(Node *current_node, Lattice *lattice);
void zh_pressure_z(Node *current_node, Lattice *lattice);
void zh_pressure_Z(Node *current_node, Lattice *lattice);

#endif