#include "data_types.h"

//#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class InfileReader {
	
	char *fname;
	ifstream data_file;

	DomainConstant *domain_constants;
	Timing *timer;
	OutputController *output_controller;
	ProjectStrings *project;
	RuntimeDomain *runtime_domain;

	void initialise()
	{
		domain_constants->c_smag = 0;
	}

	void open_file() {
		data_file.open(fname);
	}

	void close_file() {
		data_file.close();
	}

	void parse_file() {

		string line_tot;
		string line_part;
		while(getline(data_file, line_tot))
		{
			//std::getline(stringstream(line_tot),line_part,'#');
			process_data_field(line_tot);

		}
	}

	void process_data_field(string line)
	{
		// checks current line in file for data field and processes input accordingly
		string field_name;
		istringstream lineSS(line);
		if (line.find("ProjName")!=string::npos)
		{
			lineSS >> field_name >> project->name;
			cout << field_name << " = " << project->name << endl;
		}
		else if (line.find("DomainFile")!=string::npos)
		{
			lineSS >> field_name >> project->domain_fname;
			cout << field_name << " = " << project->domain_fname << endl;
		}
		else if (line.find("OutputFile")!=string::npos)
		{
			lineSS >> field_name >> project->output_fname;
			cout << field_name << " = " << project->output_fname << endl;
		}
		else if (line.find("TauMRT")!=string::npos)
		{
			string screen_label = "S_";
#if DIM < 3
			lineSS >> field_name >> domain_constants->tau_mrt[0] >> domain_constants->tau_mrt[1] >> domain_constants->tau_mrt[2] >> domain_constants->tau_mrt[3] >> domain_constants->tau_mrt[4] >> domain_constants->tau_mrt[5] >> domain_constants->tau_mrt[6] >> domain_constants->tau_mrt[7] >> domain_constants->tau_mrt[8];
#else
			lineSS >> field_name >> domain_constants->tau_mrt[0] >> domain_constants->tau_mrt[1] >> domain_constants->tau_mrt[2] >> domain_constants->tau_mrt[3] >> domain_constants->tau_mrt[4] >> domain_constants->tau_mrt[5] >> domain_constants->tau_mrt[6] >> domain_constants->tau_mrt[7] >> domain_constants->tau_mrt[8] >> domain_constants->tau_mrt[9] >> domain_constants->tau_mrt[10] >> domain_constants->tau_mrt[11] >> domain_constants->tau_mrt[12] >> domain_constants->tau_mrt[13] >> domain_constants->tau_mrt[14];
#endif		
			for(int i =0; i<Q; i++)
			{
				//lineSS >> domain_constants->tau_mrt[i];
				cout << screen_label << i << " = " << domain_constants->tau_mrt[i] << endl;
			}
		}
		else if (line.find("Tau")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->tau;
			cout << field_name << " = " << domain_constants->tau << endl;
		}
		else if (line.find("Lx")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->length[0];
			cout << field_name << " = " << domain_constants->length[0] << endl;
		}
		else if (line.find("Ly")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->length[1];
			cout << field_name << " = " << domain_constants->length[1] << endl;
		}
		else if (line.find("Lz")!=string::npos)
		{
			int tmp;
			lineSS >> field_name >> tmp;
			//check for 2d
			#if DIM > 2
				domain_constants->length[2] = tmp;
				cout << field_name << " = " << tmp << endl;
			#endif
		}
		else if (line.find("DeltaX")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->h;
			cout << field_name << " = " << domain_constants->h << endl;
		}
		else if (line.find("DeltaT")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->dt;
			cout << field_name << " = " << domain_constants->dt << endl;
		}
		else if (line.find("C_smag")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->c_smag;
			cout << field_name << " = " << domain_constants->c_smag << endl;
		}
		else if (line.find("ColType")!=string::npos)
		{
			string col_type;

			if(line.find("BGK")!=string::npos)
			{
				col_type = "BGK";
				domain_constants->collision_type = 0;
			}
			if(line.find("NTPOR")!=string::npos)
			{
				col_type = "NTPOR";
				domain_constants->collision_type = 1;
			}
			if(line.find("MRT")!=string::npos)
			{
				col_type = "MRT";
				domain_constants->collision_type = 2;
			}
			if(line.find("MRTPOR")!=string::npos)
			{
				col_type = "MRTPOR";
				domain_constants->collision_type = 3;
			}
			lineSS >> field_name;
			cout << field_name << " = " << col_type << endl;
		}
		else if (line.find("Force")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->forcing;
			cout << field_name << " = " << domain_constants->forcing << endl;
		}
		else if (line.find("MicroBC")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->micro_bc;
			cout << field_name << " = " << domain_constants->micro_bc << endl;
		}
		else if (line.find("MacroBC")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->macro_bc;
			cout << field_name << " = " << domain_constants->macro_bc << endl;
		}
		else if (line.find("Tolerance")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->tolerance;
			cout << field_name << " = " << domain_constants->tolerance << endl;
		}
		else if (line.find("Init")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->init_type;
			cout << field_name << " = " << domain_constants->init_type << endl;
		}
		else if (line.find("MaxT")!=string::npos)
		{
			lineSS >> field_name >> timer->max;
			cout << field_name << " = " << timer->max << endl;
		}
		else if (line.find("FileOut")!=string::npos)
		{
			lineSS >> field_name >> timer->plot;
			cout << field_name << " = " << timer->plot << endl;
		}
		else if (line.find("ScreenMes")!=string::npos)
		{
			lineSS >> field_name >> timer->screen;
			cout << field_name << " = " << timer->screen << endl;
		}
		else if (line.find("SteadyCheck")!=string::npos)
		{
			lineSS >> field_name >> timer->steady_check;
			cout << field_name << " = " << timer->steady_check << endl;
		}
		else if (line.find("OutputVars")!=string::npos)
		{
			bool tmp;
			lineSS >> field_name >> output_controller->u[0] >> output_controller->u[1] >> tmp >> output_controller->rho >> output_controller->pressure;
			#if DIM > 2
				if(domain_constants->length[2]>0) output_controller->u[2] = tmp;
			#endif
			if(output_controller->u[0]==true) cout << "Output VelocityX" << " = " << output_controller->u[0] << endl;
			if(output_controller->u[1]==true) cout << "Output VelocityY" << " = " << output_controller->u[1] << endl;
			#if DIM > 2
				if(output_controller->u[2]==true) cout << "Output VelocityZ" << " = " << output_controller->u[2] << endl;
			#endif
			if(output_controller->rho==true) cout << "Output Density" << " = " << output_controller->rho << endl;
			if(output_controller->pressure==true) cout << "Output Pressure" << " = " << output_controller->pressure << endl;
		}
		else if (line.find("ScreenNode")!=string::npos)
		{
			int tmp;
			lineSS >> field_name >> output_controller->screen_node[0] >> output_controller->screen_node[1] >> tmp;
			#if DIM > 2
				if(domain_constants->length[2]>0) output_controller->screen_node[2] = tmp;
			#endif
			cout << "Screen Node X" << " = " << output_controller->screen_node[0] << endl;
			cout << "Screen Node Y" << " = " << output_controller->screen_node[1] << endl;
			#if DIM > 2
				cout << "Screen Node Z" << " = " << output_controller->screen_node[2] << endl;
			#endif
		}

		else if (line.find("Interactive")!=string::npos)
		{
			lineSS >> field_name >> output_controller->interactive;
			cout << field_name << " = " << output_controller->interactive << endl;
		}

		else if (line.find("RuntimeDomain")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->runtime_domain;
			cout << field_name << " = " << domain_constants->runtime_domain << endl;
		}

		else if (line.find("Gravity")!=string::npos)
		{
			lineSS >> field_name;
			cout << field_name << " = ";

			for(int i =0; i<DIM; i++)
			{
				lineSS >> runtime_domain->gravity[i];
				cout << runtime_domain->gravity[i] << " ";
			}

			cout << endl;
		}

		else if (line.find("MicroBound")!=string::npos)
		{
			lineSS >> field_name;
			cout << field_name << " = ";

			for(int i =0; i<2*DIM; i++)
			{
				lineSS >> runtime_domain->micro_bc[i];
				cout << runtime_domain->micro_bc[i] << " ";
			}

			cout << endl;
		}

		else if (line.find("MacroBoundSpec")!=string::npos)
		{
			lineSS >> field_name;
			cout << field_name << " = ";

			for(int i =0; i<2*DIM; i++)
			{
				lineSS >> runtime_domain->macro_bc_spec[i];
				cout << runtime_domain->macro_bc_spec[i] << " ";
			}

			cout << endl;
		}

		else if (line.find("MacroBoundVal")!=string::npos)
		{
			lineSS >> field_name;
			cout << field_name << " = ";

			for(int i =0; i<2*DIM*(DIM+1); i++)
			{
				lineSS >> runtime_domain->macro_bc_val[i];
				cout << runtime_domain->macro_bc_val[i] << " ";
			}

			cout << endl;
		}

		else if (line.find("DomainWalls")!=string::npos)
		{
			lineSS >> field_name;
			cout << field_name << " = ";

			for(int i =0; i<2*DIM; i++)
			{
				lineSS >> runtime_domain->domain_walls[i];
				cout << runtime_domain->domain_walls[i] << " ";
			}

			cout << endl;
		}

		else if (line.find("GeomType")!=string::npos)
		{
			lineSS >> field_name >> runtime_domain->geom_type;
			cout << field_name << " = " << runtime_domain->geom_type << endl;
		}

		else if (line.find("GeomFile")!=string::npos)
		{
			lineSS >> field_name >> runtime_domain->geom_fname;
			cout << field_name << " = " << runtime_domain->geom_fname << endl;
		}

		else if (line.find("BinaryDomain")!=string::npos)
		{
			lineSS >> field_name >> domain_constants->binary_domain;
			cout << field_name << " = " << domain_constants->binary_domain << endl;
		}
	}

	public:
		InfileReader(char*, ProjectStrings*, DomainConstant *,Timing *,OutputController *, RuntimeDomain *);
};

InfileReader::InfileReader(char *input_filename, ProjectStrings *project_in, DomainConstant *domain_constants_in, Timing *timer_in, OutputController *output_controller_in, RuntimeDomain *runtime_domain_in) {

			fname = input_filename;
			project = project_in;
			domain_constants = domain_constants_in;
			timer = timer_in;
			output_controller = output_controller_in;
			runtime_domain = runtime_domain_in;

			cout << endl << "Reading configuration data: " << endl << endl;

			initialise();
			open_file();
			parse_file();
			close_file();

			cout << endl << "Finished reading configuration data." << endl;
}