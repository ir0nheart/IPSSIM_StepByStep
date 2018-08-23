#include "stdafx.h"
#include "Storage.h"

Storage * Storage::m_pInstance = nullptr;

Storage * Storage::instance()
{
	if (m_pInstance == nullptr)
	{
		m_pInstance = new Storage;
	}
	return m_pInstance;
}

std::vector<char> Storage::get_mesh_dim_string() const
{
	return mesh_dim_string;
}

std::vector<char> Storage::get_mesh_type_string() const
{
	return mesh_type_string;
}

std::vector<char> Storage::get_nn1_string() const
{
	return nn1_string;
}
std::vector<char> Storage::get_nn2_string() const
{
	return nn2_string;
}

std::vector<char> Storage::get_nn3_string() const
{
	return nn3_string;
}

std::vector<char> Storage::get_flow_type_string() const
{
	return flow_type_string;
}

std::vector<char> Storage::get_transport_type_string() const
{
	return transport_type_string;
}

std::vector<char> Storage::get_simulation_condition_string() const
{
	return simulation_condition_string;
}

std::vector<char> Storage::get_simulation_start_string() const
{
	return simulation_start_string;
}

int Storage::get_nn() const
{
	return NN;
}

void Storage::set_nn(int nn)
{
	NN = nn;
}

int Storage::get_ne() const
{
	return NE;
}

void Storage::set_ne(int ne)
{
	NE = ne;
}

int Storage::get_nsop() const
{
	return NSOP;
}

void Storage::set_nsop(int nsop)
{
	NSOP = nsop;
}

int Storage::get_nsou() const
{
	return NSOU;
}

void Storage::set_nsou(int nsou)
{
	NSOU = nsou;
}

int Storage::get_npbc() const
{
	return NPBC;
}

void Storage::set_npbc(int npbc)
{
	NPBC = npbc;
}

int Storage::get_nubc() const
{
	return NUBC;
}

void Storage::set_nubc(int nubc)
{
	NUBC = nubc;
}

int Storage::get_nobs() const
{
	return NOBS;
}

void Storage::set_nobs(int nobs)
{
	NOBS = nobs;
}

int Storage::get_solution_storage() const
{
	return solution_storage;
}

void Storage::set_solution_storage(int val)
{
	solution_storage = val;
}

double Storage::get_up() const
{
	return UP;
}

double Storage::get_gnup() const
{
	return GNUP;
}

double Storage::get_gnuu() const
{
	return GNUU;
}
int Storage::get_npcyc() const
{
	return NPCYC;
}
int Storage::get_nsch() const
{
	return NSCH;
}

int Storage::get_nucyc() const
{
	return NUCYC;
}

int Storage::get_rpmax() const
{
	return RPMAX;
}

int Storage::get_rumax() const
{
	return RUMAX;
}

int Storage::get_itrmax() const
{
	return ITRMAX;
}

int Storage::get_max_p_iterations() const
{
	return max_p_iterations;
}

double Storage::get_p_tolerance() const
{
	return p_tolerance;
}

double Storage::get_u_tolerance() const
{
	return u_tolerance;
}

int Storage::get_max_u_iterations() const
{
	return max_u_iterations;
}

std::vector<char> Storage::get_p_solver_string() const
{
	return p_solver_string;
}

std::vector<char> Storage::get_u_solver_string() const
{
	return u_solver_string;
}

std::vector<char> Storage::get_temporal_control(int indx)
{
	return temporalControl[indx];
}

std::vector<char> Storage::get_simulation_output_controls(int indx) const
{
	return simulation_output_controls[indx];
}

std::vector<char> Storage::get_node_output_controls(int indx) const
{
	return node_output_headers[indx];
}

std::vector<char> Storage::get_element_output_controls(int indx) const
{
	return element_output_headers[indx];
}

int Storage::get_element_output_every() const
{
	return element_output_every;
}

int Storage::get_node_output_every() const
{
	return node_output_every;
}

int Storage::get_simulation_output_every() const
{
	return simulation_output_every;
}

int Storage::get_observation_first_line() const
{
	return observation_first_line;
}

int Storage::get_nbcfpr() const
{
	return NBCFPR;
}

int Storage::get_nbcspr() const
{
	return NBCSPR;
}

int Storage::get_nbcppr() const
{
	return NBCPPR;
}

int Storage::get_nbcupr() const
{
	return NBCUPR;
}

char Storage::get_cinact() const
{
	return CINACT;
}

double Storage::get_compfl() const
{
	return COMPFL;
}

double Storage::get_cw() const
{
	return CW;
}

double Storage::get_sigmaw() const
{
	return SIGMAW;
}

double Storage::get_rhow0() const
{
	return RHOW0;
}

double Storage::get_urhow0() const
{
	return URHOW0;
}

double Storage::get_visc0() const
{
	return VISC0;
}

double Storage::get_drwdu() const
{
	return DRWDU;
}

double Storage::get_compma() const
{
	return COMPMA;
}

double Storage::get_cs() const
{
	return CS;
}

double Storage::get_rhos() const
{
	return RHOS;
}

std::vector<char> Storage::get_adsorption_string() const
{
	return adsorption_string;
}

double Storage::get_sigmas() const
{
	return SIGMAS;
}

double Storage::get_prodf0() const
{
	return PRODF0;
}

double Storage::get_prodf1() const
{
	return PRODF1;
}

double Storage::get_prods1() const
{
	return PRODS1;
}

double Storage::get_gravx() const
{
	return GRAVX;
}

double Storage::get_gravy() const
{
	return GRAVY;
}

double Storage::get_gravz() const
{
	return GRAVZ;
}

double Storage::get_scaly() const
{
	return SCALY;
}

double Storage::get_scalz() const
{
	return SCALZ;
}

double Storage::get_scalx() const
{
	return SCALX;
}

double Storage::get_porfac() const
{
	return PORFAC;
}

std::vector<std::vector<char>> Storage::get_element_props() const
{
	return element_props;
}

double Storage::get_prods0() const
{
	return PRODS0;
}




std::vector<char> Storage::get_sutra_string() const
{
	return sutra_string;
}

std::vector<char> Storage::get_version_string() const
{
	return version_string;
}

std::vector<char> Storage::get_version_num_string() const
{
	return version_num_string;
}

std::vector<char> Storage::get_simulation_type_string() const
{
	return simulation_type_string;
}

std::vector<char> Storage::get_transport_string() const
{
	return transport_string;
}

Storage::Storage()
{
	std::vector<std::string> dataSets{ "DataSet_1", "DataSet_2A", "DataSet_2B", "DataSet_3", "DataSet_4", "DataSet_5", "DataSet_6",
		"DataSet_7A", "DataSet_7B", "DataSet_7C", "DataSet_8ABC", "DataSet_8D", "DataSet_8E_9_10_11", "DataSet_12_13_14A",
		"DataSet_14B", "DataSet_15A", "DataSet_15B", "DataSet_17", "DataSet_18", "DataSet_19", "DataSet_20", "DataSet_22" };
	for (int i = 0; i < 22; i++){
		DataSet* dataSet = new DataSet(dataSets[i]);
		dataSetMap[dataSets[i]] = dataSet;
	}
}


Storage::~Storage()
{
}

void Storage::addTittle(std::string str)
{
	titles.push_back(str);
}

std::string Storage::getTittle(int index)
{
	if (index >= 0)
	{
		if (index < titles.size())
		{
			return titles[index];
		} else
		{
			return "There is no such title.";
		}
	} else
	{
		return "Title index cannot be negative.";
	}
}