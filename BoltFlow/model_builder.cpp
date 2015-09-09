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
		InfileReader infile_reader(fname_config, project_t, domain_constants, time_t, output_controller);
		
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
	}

// Allocates memory for variables which have variable size due to problem geometry
	void variable_size_allocator()
	{
		domain_size = 1;
		for(int d = 0; d<DIM; d++)
		{
			domain_size = domain_size*domain_constants->length[d];
		}
		int domain_data_size;
		domain_data_size = domain_size*sizeof(double);

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

	void variable_loader()
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

		if(domain_constants->init_type == 0)
		{
			load_static_IC();
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
	ModelBuilder (char *, Lattice*, DomainConstant*, Domain*, OutputController*, Timing*, ProjectStrings*);

	ModelBuilder ();

};

ModelBuilder::ModelBuilder (char *input_filename, Lattice *lattice, DomainConstant *domain_constants, Domain *domain, OutputController *output_controller, Timing *time, ProjectStrings *project) 
{
	this->lattice= lattice;
	this->domain_constants= domain_constants;
	this->domain= domain;
	this->output_controller= output_controller;
	time_t = time;
	project_t = project;

	fname_config = input_filename;
	constant_size_allocator();
	constant_loader();
	variable_size_allocator();
	variable_assembler();
	cout << "variable assembler complete" << endl;
	variable_loader();
	cout << "variable loader complete" << endl;
}

ModelBuilder::ModelBuilder (){}

#endif