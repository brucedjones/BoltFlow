#ifndef MICRO_BC_H
#define MICRO_BC_H

#ifdef D2Q9
	#include "d2q9/d2q9_zh_defs.h"
	//#include "d2q9/d2q9_sf_defs.h"
	const micro_condition micro_conditions[4] = { zh_pressure_x, zh_pressure_X, NULL, NULL};/*,
																	sf_x, sf_X, sf_y, sf_Y};*/
#endif

#ifdef D3Q15
	#include "d3q15/d3q15_zh_defs.h"
	//#include "d3q15/d3q15_sf_defs.h"
	const micro_condition micro_conditions[6] = { zh_pressure_x, zh_pressure_X, zh_pressure_y, zh_pressure_Y, zh_pressure_z, zh_pressure_Z };/* ,
																		sf_x, sf_X, sf_y, sf_Y, sf_z, sf_Z};*/
#endif

#endif