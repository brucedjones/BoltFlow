#ifndef MODEL_BUILDER
#define MODEL_BUILDER

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
#include "infile_reader.cpp"
#include "cgns/cgns_input_handler.cpp"
using namespace std;

class ModelBuilder
{
	int length[DIM];
	// HOST VARIABLE DECLARATION
	Timing *time_t;
	ProjectStrings *project_t;
	OutputController *output_controller;
	Lattice *lattice;
	Domain *domain;
	DomainConstant *domain_constants;
	RuntimeDomain *runtime_domain;
	double **f, *rho, **u, *geometry, **force;
	int *micro_bc;
	int *macro_bc;

	// SCALAR DECLARATION (PLATFORM AGNOSTIC)
	double tau, residual;
	double tolerance;
	int domain_size, maxT, saveT, steadyT, collision_type;

	// CONFIG FLAGS AND STRINGS
	char *fname_config;
	bool zhoue;
	bool forcing;
	bool is2D;

// Allocates memory for variables which are constant in size
	void constant_size_allocator()
	{
		// Allocate container structures
		//combi_malloc<Lattice>(&lattice, &lattice_d, sizeof(Lattice));
		//combi_malloc<Domain>(&domain, &domain_d, sizeof(Domain));
		//combi_malloc<DomainConstant>(&domain_constants, &domain_constants_d, sizeof(DomainConstant));
		//combi_malloc<OutputController>(&output_controller, &output_controller_d, sizeof(OutputController));
		//domain_constants = (DomainConstant *)malloc(sizeof(DomainConstant));
		//time_t = (Timing *)malloc(sizeof(Timing));
		//project_t = (ProjectStrings *)malloc(sizeof(ProjectStrings));
	}

	void constant_loader()
	{
		// LOAD CONSTANTS FROM FILE
		InfileReader infile_reader(fname_config, project_t, domain_constants, time_t, output_controller, runtime_domain);
		
		// LOAD LATTICE CONSTANTS
		LOAD_E(domain_constants->e);
		LOAD_OMEGA(domain_constants->omega);
		LOAD_OPP(domain_constants->opp);
		LOAD_M(domain_constants->M);
		LOAD_M_INV(domain_constants->M_inv);
		for(int i =0;i<NUM_RESIDS;i++)
		{
			domain_constants->residual[i] = 1;
		}

		domain_size = 1;
		for(int d = 0; d<DIM; d++)
		{
			domain_size = domain_size*domain_constants->length[d];
		}

		// Check for runtime domain and flip switches accordingly
		if(domain_constants->runtime_domain)
		{
			bool forcing = false;
			for(int i = 0; i<DIM; i++)
			{
				if(runtime_domain->gravity[i]!=0) forcing = true;
			}
			domain_constants->forcing = forcing;

			bool micro_bc = false;
			for(int i = 0; i<2*DIM; i++)
			{
				if(runtime_domain->micro_bc[i]>0) micro_bc = true;
			}
			domain_constants->micro_bc = micro_bc;

			bool macro_bc = false;
			for(int i = 0; i<2*DIM; i++)
			{
				if(runtime_domain->macro_bc_spec[i]>0) macro_bc = true;
			}
			domain_constants->macro_bc = macro_bc;

		}
	}

// Allocates memory for variables which have variable size due to problem geometry
	void variable_size_allocator()
	{
		int domain_data_size = domain_size*sizeof(double);

		// Allocate required arrays
		// PDFS
		double *f_tmp[Q];
		f = (double**)malloc(sizeof(double*)*Q);
		for(int i=0;i<Q;i++)
		{
			f[i] = (double*)malloc(domain_data_size);
		}
		
		// RHO
		rho = (double*)malloc(domain_data_size);
		
		// VELOCITY
		double *u_tmp[DIM];
		u = (double**)malloc(sizeof(double*)*DIM);
		for(int i=0;i<DIM;i++)
		{
			u[i] = (double*)malloc(domain_data_size);
		}

		// GEOMETRY
		geometry = (double*)malloc(domain_data_size);
		
		// ALLOCATE OPTION ARRAYS
		// FORCING
		if(domain_constants->forcing == true)
		{
			force = (double**)malloc(sizeof(double*)*DIM);
			for(int i=0;i<DIM;i++)
			{
				force[i] = (double*)malloc(domain_data_size);
			}
		}

		// MICRO BC
		if(domain_constants->micro_bc == true)
		{
			micro_bc = (int*)malloc(domain_size*sizeof(int));
		}

		// MACRO BC
		if(domain_constants->macro_bc == true)
		{
			macro_bc = (int*)malloc(domain_size*sizeof(int));
		}
	}

	void variable_assembler()
	{
		lattice->f = f;

		domain->micro_bc = micro_bc;
		domain->macro_bc = macro_bc;
		domain->geometry = geometry;
		domain->force = force;
		domain->u = u;
		domain->rho = rho;
	}

	void runtime_domain_loader()
	{
		// Gravity
		if(domain_constants->forcing)
		{
			for(int idx=0; idx<domain_size; idx++)
			{
				for(int d=0; d<DIM; d++)
				{
					domain->force[d][idx] = runtime_domain->gravity[d];
				}
			}
		}

		// MICRO BOUNDARY CONDITION SPECIFICATION
		if(domain_constants->micro_bc)
		{
#if DIM < 3
			for(int j=0; j<domain_constants->length[1]; j++)
			{
				for(int i=0; i<domain_constants->length[0]; i++)
				{	
					int caseApplied = -1;

					if(j==0 && runtime_domain->micro_bc[2]>0) caseApplied = 2;
					if(j==domain_constants->length[1]-1 && runtime_domain->micro_bc[3]>0) caseApplied = 3;
					if(i==0 && runtime_domain->micro_bc[0]>0) caseApplied = 0;
					if(i==domain_constants->length[0]-1 && runtime_domain->micro_bc[1]>0) caseApplied = 1;

					int idx = i+j*domain_constants->length[0];

					if(caseApplied>=0)
					{
						domain->micro_bc[idx] = runtime_domain->micro_bc[caseApplied];
					} else {
						domain->micro_bc[idx] = 0;
					}
				}
			}
#else
			for(int k=0; k<domain_constants->length[2]; k++)
			{
				for(int j=0; j<domain_constants->length[1]; j++)
				{
					for(int i=0; i<domain_constants->length[0]; i++)
					{
						int caseApplied = -1;

						if(k==0 && runtime_domain->micro_bc[4]>0) caseApplied = 4;
						if(k==domain_constants->length[2]-1 && runtime_domain->micro_bc[5]>0) caseApplied = 5;
						if(j==0 && runtime_domain->micro_bc[2]>0) caseApplied = 2;
						if(j==domain_constants->length[1]-1 && runtime_domain->micro_bc[3]>0) caseApplied = 3;
						if(i==0 && runtime_domain->micro_bc[0]>0) caseApplied = 0;
						if(i==domain_constants->length[0]-1 && runtime_domain->micro_bc[1]>0) caseApplied = 1;
						
						int idx = i+j*domain_constants->length[0]+ k*domain_constants->length[0]*domain_constants->length[1];

						if(caseApplied>=0)
						{
							domain->micro_bc[idx] = runtime_domain->micro_bc[caseApplied];
						} else {
							domain->micro_bc[idx] = 0;
						}				
					}
				}
			}
#endif
		}

		// MACRO BOUNDARY CONDITION SPECIFICATION
		if(domain_constants->macro_bc)
		{
#if DIM < 3
			for(int j=0; j<domain_constants->length[1]; j++)
			{
				for(int i=0; i<domain_constants->length[0]; i++)
				{	
					int caseApplied = -1;

					if(j==0 && runtime_domain->macro_bc_spec[2]>0) caseApplied = 2;
					if(j==domain_constants->length[1]-1 && runtime_domain->macro_bc_spec[3]>0) caseApplied = 3;
					if(i==0 && runtime_domain->macro_bc_spec[0]>0) caseApplied = 0;
					if(i==domain_constants->length[0]-1 && runtime_domain->macro_bc_spec[1]>0) caseApplied = 1;

					int idx = i+j*domain_constants->length[0];

					if(caseApplied>=0)
					{
						domain->macro_bc[idx] = runtime_domain->macro_bc_spec[caseApplied];
						domain->rho[idx] = runtime_domain->macro_bc_val[0+caseApplied*(DIM+1)];
						domain->u[0][idx] = runtime_domain->macro_bc_val[1+caseApplied*(DIM+1)];
						domain->u[1][idx] = runtime_domain->macro_bc_val[2+caseApplied*(DIM+1)];
					} else {
						domain->macro_bc[idx] = 0;
					}
				}
			}
#else
			for(int k=0; k<domain_constants->length[2]; k++)
			{
				for(int j=0; j<domain_constants->length[1]; j++)
				{
					for(int i=0; i<domain_constants->length[0]; i++)
					{
						int caseApplied = -1;

						if(k==0 && runtime_domain->macro_bc_spec[4]>0) caseApplied = 4;
						if(k==domain_constants->length[2]-1 && runtime_domain->macro_bc_spec[5]>0) caseApplied = 5;
						if(j==0 && runtime_domain->macro_bc_spec[2]>0) caseApplied = 2;
						if(j==domain_constants->length[1]-1 && runtime_domain->macro_bc_spec[3]>0) caseApplied = 3;
						if(i==0 && runtime_domain->macro_bc_spec[0]>0) caseApplied = 0;
						if(i==domain_constants->length[0]-1 && runtime_domain->macro_bc_spec[1]>0) caseApplied = 1;
						
						int idx = i+j*domain_constants->length[0]+ k*domain_constants->length[0]*domain_constants->length[1];

						if(caseApplied>=0)
						{
							domain->macro_bc[idx] = runtime_domain->macro_bc_spec[caseApplied];
							domain->rho[idx] = runtime_domain->macro_bc_val[0+caseApplied*(DIM+1)];
							domain->u[0][idx] = runtime_domain->macro_bc_val[1+caseApplied*(DIM+1)];
							domain->u[1][idx] = runtime_domain->macro_bc_val[2+caseApplied*(DIM+1)];
							domain->u[2][idx] = runtime_domain->macro_bc_val[3+caseApplied*(DIM+1)];
						} else {
							domain->macro_bc[idx] = 0;
						}				
					}
				}
			}
#endif
		}

		// Domain Walls
#if DIM < 3
		for(int j=0; j<domain_constants->length[1]; j++)
		{
			for(int i=0; i<domain_constants->length[0]; i++)
			{	
				double geometry = 0.0;
				
				if(j==0 && runtime_domain->domain_walls[2]>0) geometry=1.0;
				if(j==domain_constants->length[1]-1 && runtime_domain->domain_walls[3]>0) geometry=1.0;
				if(i==0 && runtime_domain->domain_walls[0]>0) geometry=1.0;
				if(i==domain_constants->length[0]-1 && runtime_domain->domain_walls[1]>0) geometry=1.0;
				
				int idx = i+j*domain_constants->length[0];
				domain->geometry[idx] = geometry;
			}
		}
#else
		for(int k=0; k<domain_constants->length[2]; k++)
		{
			for(int j=0; j<domain_constants->length[1]; j++)
			{
				for(int i=0; i<domain_constants->length[0]; i++)
				{
					double geometry = 0.0;

					if(k==0 && runtime_domain->domain_walls[4]>0) geometry=1.0;
					if(k==domain_constants->length[2]-1 && runtime_domain->domain_walls[5]>0) geometry=1.0;
					if(j==0 && runtime_domain->domain_walls[2]>0) geometry=1.0;
					if(j==domain_constants->length[1]-1 && runtime_domain->domain_walls[3]>0) geometry=1.0;
					if(i==0 && runtime_domain->domain_walls[0]>0) geometry=1.0;
					if(i==domain_constants->length[0]-1 && runtime_domain->domain_walls[1]>0) geometry=1.0;

					int idx = i+j*domain_constants->length[0];
					domain->geometry[idx] = geometry;
				}
			}
		}
#endif
	}

	void binary_domain_loader()
	{
		// LOAD GEOMETRY
		CGNSInputHandler input_handler(project_t->domain_fname, domain_constants->length);
		input_handler.read_field(domain->geometry, "Porosity");

		// LOAD FORCES IF REQUIRED
		if(domain_constants->forcing == true)
		{
			char force_labels[3][33];
			strcpy(force_labels[0], "ForceX");
			strcpy(force_labels[1], "ForceY");
			strcpy(force_labels[2], "ForceZ");

			for(int d=0;d<DIM;d++)
			{
				input_handler.read_field(domain->force[d], force_labels[d]);
			}
		}

		// LOAD MICRO BOUNDARY CONDITIONS IF REQUIRED
		if(domain_constants->micro_bc == true)
		{
			input_handler.read_field(domain->micro_bc, "MicroBC");
		}

		// LOAD MACRO BOUNDARY CONDITIONS IF REQUIRED
		if(domain_constants->macro_bc == true)
		{
			input_handler.read_field(domain->macro_bc, "MacroBC");
			
			char vel_labels[3][33];
			strcpy(vel_labels[0], "VelocityX");
			strcpy(vel_labels[1], "VelocityY");
			strcpy(vel_labels[2], "VelocityZ");

			for(int d=0;d<DIM;d++)
			{
				input_handler.read_field(domain->u[d], vel_labels[d]);
			}

			input_handler.read_field(domain->rho, "Rho");
		}
	}

	void load_static_IC()
	{
		double omega[Q];
		LOAD_OMEGA(omega);
		for(int i=0;i<Q;i++)
		{
			for(int index=0;index<(domain_size);index++)
			{
				lattice->f[i][index] = 1.0*omega[i];
			}
		}
	}

public:
	ModelBuilder (char *, Lattice*, DomainConstant*, Domain*, OutputController*, Timing*, ProjectStrings*, RuntimeDomain*);

	ModelBuilder ();

};

ModelBuilder::ModelBuilder (char *input_filename, Lattice *lattice, DomainConstant *domain_constants, Domain *domain, OutputController *output_controller, Timing *time, ProjectStrings *project, RuntimeDomain *runtime_domain) 
{
	this->lattice = lattice;
	this->domain_constants = domain_constants;
	this->domain = domain;
	this->runtime_domain = runtime_domain;
	this->output_controller= output_controller;
	time_t = time;
	project_t = project;

	fname_config = input_filename;
	constant_size_allocator();
	constant_loader();
	variable_size_allocator();
	variable_assembler();
	cout << "variable assembler complete" << endl;

	if(domain_constants->init_type == 0)
	{
		load_static_IC();
		cout << "initialised to static domain" << endl;
	}

	if(this->domain_constants->runtime_domain)
	{
		runtime_domain_loader();
		cout << "runtime domain initialisation complete" << endl;
	}

	if(this->domain_constants->binary_domain)
	{
		binary_domain_loader();
		cout << "binary domain initialisation complete" << endl;
	}
	cout << "variable loader complete" << endl;
}

ModelBuilder::ModelBuilder (){}

#endif