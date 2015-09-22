#ifndef KERNEL
#define KERNEL
////////////////////////////////////////////////////////////////////////////////
//
// LBM-C
// A lattice Boltzmann fluid flow solver written using cppDA
//
// Copyright (C) 2011  Bruce Jones
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTIcppLAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// D2Q9 Lattice configuration:
//
//       6   2   5
//        \  |  /
//         \ | /
//          \|/
//       3---0---1
//          /|\
//         / | \
//        /  |  \
//       7   4   8
//
///////////////////////////////////////////////////////////////////////////////

#ifdef OS_WINDOWS
   #ifdef _WIN64
		//#pragma comment(lib, "cgns/x64/lib/cgns.lib")
		#include "cgns\x64\include\cgnslib.h"
		//#pragma comment(lib, "HDF5/x64/lib/hdf5.lib")
		#include "HDF5/x64/include/hdf5.h"
		//#pragma comment(lib, "HDF5/x64/lib/libszip.lib")
		#include "HDF5/x64/include/szlib.h"
		//#pragma comment(lib, "HDF5/x64/lib/libzlib.lib")
		#include "HDF5/x64/include/zlib.h"
	#elif _WIN32
		#pragma comment(lib, "cgns/x86/lib/cgns.lib")
		#include "cgns\x86\include\cgnslib.h"
		#pragma comment(lib, "HDF5/x86/lib/hdf5.lib")
		#include "HDF5/x86/include/hdf5.h"
		#pragma comment(lib, "HDF5/x86/lib/libszip.lib")
		#include "HDF5/x86/include/szlib.h"
		#pragma comment(lib, "HDF5/x86/lib/libzlib.lib")
		#include "HDF5/x86/include/zlib.h"
	#endif
#else //UNIX
	#include <cgnslib.h>
	#include <hdf5.h>
	#include <szlib.h>
	#include <zlib.h>
#endif


#include <stdio.h>
#include <time.h>
#include "data_types.h"
#include "macros.cpp"
#include "solver.cpp"
#include "index.h"
#include "model_builder.cpp"
#include "cgns/cgns_output_handler.cpp"

// HOST VARIABLE DECLARATION
Lattice *lattice;
Domain *domain;
DomainConstant *domain_constants;
OutputController *output_controller;
Timing *times;
ProjectStrings *project;
ModelBuilder model_builder;


// SCALAR DECLARATION (PLATFORM AGNOSTIC)
bool store_macros = false;

// DECLARE OUTPUT HANDLER
CGNSOutputHandler output_handler;

int main(int argc, char **argv)
{	
	// Initialise memory for LBM model
	setup(argv[1]);

	// Report domain configuration
	printf("X-Length:		%d\n", domain_constants->length[0]);
	printf("Y-Length:		%d\n", domain_constants->length[1]);
	#if DIM > 2
		printf("Z-Length:		%d\n", domain_constants->length[2]);
	#endif
	printf("Relaxation Time (Tau):	%f\n", domain_constants->tau);
	printf("\nPress return to continue...");
	if (output_controller->interactive == true) getchar();

	// Get current clock cycle number
	clock_t t1=clock();
	clock_t t1Mess = t1;
	clock_t t2Mess = t1;

	int domain_size=1;
	int stop=0;
	for(int d = 0; d<DIM ;d++)
	{
		domain_size = domain_size*domain_constants->length[d];
	}

	for(int i = 1; i<times->max+1; i++)
	{
		if((times->plot>0 && i%times->plot == 0) ||
		   (times->steady_check>0 && i%times->steady_check) || 
		   (times->screen>0 && i%times->screen)) store_macros = true;

		iterate(i-1);

		if(times->plot>0 && i%times->plot == 0)
		{
			output_macros(i);
			store_macros = false;
		}

		if(times->screen>0 && i%times->screen == 0)
		{
			t2Mess = clock();
			double timePerIter = (((double)t2Mess - (double)t1Mess) / (double)CLOCKS_PER_SEC) / (double)times->screen;
			double lups = domain_size / timePerIter;
			t1Mess = t2Mess;
			screen_mess(i,output_controller->screen_node,lups);
			store_macros = false;
		}

		if(times->steady_check>0 && i%times->steady_check == 0)
		{
			compute_residual(i);
			
			for(int resid=0;resid<NUM_RESIDS;resid++)
			{
				if(domain_constants->residual[resid]<domain_constants->tolerance) stop += 1;
			}
			if(isIndeterminate(domain_constants->residual[i%NUM_RESIDS]))
			{
				output_macros(i);
				exit(1);
			} else if(stop==NUM_RESIDS)
			{
				output_macros(i);
				break;
			}
			stop = 0;
			store_macros = false;
		}
	}

	// Get current clock cycle number
	clock_t t2=clock();
	// Compare and report global execpption time
	double cputime = ((double)t2-(double)t1)/(double)CLOCKS_PER_SEC;
	printf("\n\nTotal Run Time: %fs",cputime);
	printf("\nPress return to finish");
	if (output_controller->interactive == true) getchar();
}


// EXEcppTES ALL ROUTINES REQUIRED FOR THE MODEL SET UP
void setup(char *data_file)
{
	// Allocate container structures
	lattice = (Lattice*)malloc(sizeof(Lattice));
	domain = (Domain*)malloc(sizeof(Domain));
	domain_constants = (DomainConstant*)malloc(sizeof(DomainConstant));
	output_controller = (OutputController*)malloc(sizeof(OutputController));
	times = (Timing *)malloc(sizeof(Timing));
	project = (ProjectStrings *)malloc(sizeof(ProjectStrings));

	ModelBuilder tmpmb(data_file, lattice, domain_constants, domain, output_controller, times, project);
	model_builder = tmpmb;

	int z_len = 1;
	#if DIM > 2
		z_len = domain_constants->length[2];
	#endif
	CGNSOutputHandler tmp(project->output_fname,domain_constants->length[0],domain_constants->length[1],z_len);
	output_handler = tmp;
}



// COPIES f_i DATA FROM DEVICE TO HOST AND COMPUTERS MACROSCOPIC VALUES ON HOST, THIS DATA
// IS THEN WRITTEN TO THE OUTPUT FILE
//
// Note:	A computationally more efficient implementation would compute macroscopic
//			value's on the gpu and then just copy that data, this would however consume
//			more memory
void output_macros(int time)
{
	int domain_size = domain_constants->length[0]*domain_constants->length[1];
	#if DIM > 2
		domain_size = domain_size*domain_constants->length[2];
	#endif

	int num_fields = 0;
	if (output_controller->u[0] == true) num_fields++;
	if (output_controller->u[1] == true) num_fields++;
#if DIM > 2
	if (output_controller->u[2] == true) num_fields++;
#endif
	if (output_controller->rho == true) num_fields++;

	char **labels;
	double **data;

	labels = (char **)malloc(num_fields * sizeof (char *));
	data = (double **)malloc(num_fields * sizeof(double));

	for(int i = 0; i<num_fields;i++)
	{
		labels[i] = (char *)malloc(STR_LENGTH*sizeof(char));
	}

	int counter = 0;

	if (output_controller->u[0] == true)
	{
		data[counter] = domain->u[0];
		strcpy(labels[counter],"VelocityX");
		counter++;
	}

	if (output_controller->u[1] == true)
	{
		data[counter] = domain->u[1];
		strcpy(labels[counter],"VelocityY");
		counter++;
	}
#if DIM > 2
	if (output_controller->u[2] == true)
	{
		data[counter] = domain->u[2];
		strcpy(labels[counter],"VelocityZ");
		counter++;
	}
#endif	
	if (output_controller->rho == true)
	{
		data[counter] = domain->rho;
		strcpy(labels[counter],"Density");
		counter++;
	}

	output_handler.append_solution_output(time,num_fields,data,labels);
}

// CONFIGURES THE KERNEL CONFIGURATION AND LAUNCHES KERNEL
void iterate(int t)
{
	// ITERATE ONCE
#pragma omp parallel for
	for (int k = 0; k < domain_constants->length[2]; k++)
	{
		for (int j = 0; j < domain_constants->length[1]; j++)
		{
			for (int i = 0; i < domain_constants->length[0]; i++)
			{
				iterate_kernel(lattice, domain, domain_constants, store_macros, t, i, j, k);
			}
		}
	}
}

double current_RMS(Domain *domain)
{
	double curr_RMS = 0;
	for (int k = 0; k < domain_constants->length[2]; k++)
	{
		for (int j = 0; j < domain_constants->length[1]; j++)
		{
			for (int i = 0; i < domain_constants->length[0]; i++)
			{
				int ixd = i + j*domain_constants->length[0] + k*domain_constants->length[0] * domain_constants->length[1];
				double rho = domain->rho[ixd];
				double u = 0;
				for (int d = 0; d < DIM; d++)
				{
					u += domain->u[d][ixd] * domain->u[d][ixd];
				}
				double energy = 0.5*rho*u; // u is velocity squared here
				curr_RMS += energy;	
			}
		}
	}
	return sqrt(curr_RMS / (domain_constants->length[0] * domain_constants->length[1] * domain_constants->length[2]));
}

double prev_RMS = 0;

double error_RMS(Domain *domain)
{
	double curr_RMS = current_RMS(domain);
	double tmp = ((abs(curr_RMS-prev_RMS)/times->steady_check))/curr_RMS;

	prev_RMS = curr_RMS;

	return tmp;
}

void compute_residual(int time)
{
	domain_constants->residual[time%NUM_RESIDS] = error_RMS(domain);
}

void screen_mess(int iter, int coord[DIM], double lups)
{
	int idx = coord[0]+coord[1]*domain_constants->length[0];
	#if DIM > 2
		idx += coord[2]*domain_constants->length[0]*domain_constants->length[1];
	#endif

		cout << "time = " << iter << "; Lup/s = " << lups << "; rho = " << domain->rho[idx] << "; uX = " << domain->u[0][idx] << "; uY = " << domain->u[1][idx] << "; ";
	#if DIM>2
		cout << "uZ = " << domain->u[2][idx] << "; ";
	#endif
	cout << "resid = " << domain_constants->residual[iter%NUM_RESIDS] << endl;
}

bool isIndeterminate(const double pV)
{
    return (pV != pV);
} 

#endif