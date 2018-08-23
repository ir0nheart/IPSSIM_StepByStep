#ifndef STORAGE_H
#define STORAGE_H
#pragma once
#include "DataSet.h"

class Storage
{
public:
	static Storage * instance();
	void addTittle(std::string str);
	std::string getTittle(int index);
	void set_sutra_string(std::vector<char> str){ sutra_string = str; }
	void set_version_string(std::vector<char> str){ version_string = str; }
	void set_version_num_string(std::vector<char> str){ version_num_string = str; }
	void set_simulation_type_string(std::vector<char> str){ simulation_type_string = str; }
	void set_transport_string(std::vector<char> str){ transport_string = str; }
	void set_mesh_dim_string(std::vector<char> str){ mesh_dim_string = str; }
	void set_mesh_type_string(std::vector<char> str){ mesh_type_string = str; }
	void set_nn1_string(std::vector<char> str){ nn1_string = str; }
	void set_nn2_string(std::vector<char> str){ nn2_string = str; }
	void set_nn3_string(std::vector<char> str){ nn3_string = str; }
	void set_flow_type_string(std::vector<char> str){flow_type_string = str;}
	void set_transport_type_string(std::vector<char> str){ transport_type_string = str; }
	void set_simulation_condition_string(std::vector<char> str){ simulation_condition_string = str; }
	void set_simulation_start_string(std::vector<char> str){ simulation_start_string = str; }
	std::vector<char> get_sutra_string() const;
	std::vector<char> get_version_string() const;
	std::vector<char> get_version_num_string() const;
	std::vector<char> get_simulation_type_string() const;
	std::vector<char> get_transport_string() const;
	std::vector<char> get_mesh_dim_string() const;
	std::vector<char> get_mesh_type_string() const;
	std::vector<char> get_nn1_string() const;
	std::vector<char> get_nn2_string() const;
	std::vector<char> get_nn3_string() const;
	std::vector<char> get_flow_type_string() const;;
	std::vector<char> get_transport_type_string() const;;
	std::vector<char> get_simulation_condition_string() const;;
	std::vector<char> get_simulation_start_string() const;;
	int get_nn() const;
	void set_nn(int nn);
	int get_ne() const;
	void set_ne(int ne);
	int get_nsop() const;
	void set_nsop(int nsop);
	int get_nsou() const;
	void set_nsou(int nsou);
	int get_npbc() const;
	void set_npbc(int npbc);
	int get_nubc() const;
	void set_nubc(int nubc);
	int get_nobs() const;
	void set_nobs(int nobs);
	int get_solution_storage() const;
	void set_solution_storage(int val);
	void set_up(double up){UP = up;}
	void set_gnup(double gnup){GNUP = gnup;}
	void set_gnuu(double gnuu){GNUU = gnuu;}
	void reserveNodeData(){ nodeData.reserve(NN); }
	double get_up() const;
	double get_gnup() const;
	double get_gnuu() const;
	void set_nsch(int nsch){ NSCH = nsch; }
	void set_npcyc(int npcyc){ NPCYC = npcyc; }
	void set_nucyc(int nucyc){ NUCYC = nucyc; }
	int get_nsch() const;
	int get_npcyc() const;
	int get_nucyc() const;
	int get_rpmax() const;
	int get_rumax() const;
	int get_itrmax() const;
	int get_max_p_iterations() const;
	double get_p_tolerance() const;
	int get_max_u_iterations() const;
	double get_u_tolerance() const;
	void set_max_p_iterations(int val){ max_p_iterations = val; }
	void set_p_tolerance(double val){ p_tolerance = val; }
	void set_max_u_iterations(int val){ max_u_iterations = val; }
	void set_u_tolerance(double val){ u_tolerance = val; }
	void set_rpmax(int rpmax){ RPMAX = rpmax; }
	void set_rumax(int rumax){ RUMAX = rumax; }
	void set_itrmax(int itrmax){ ITRMAX = itrmax; }
	void add_temporal_control(std::vector<char> str){ temporalControl.push_back(str); }
	void set_p_solver_string(std::vector<char>str){ p_solver_string = str; }
	void set_u_solver_string(std::vector<char>str){ u_solver_string = str; }
	std::vector<char> get_p_solver_string() const;
	std::vector<char> get_u_solver_string() const;
	std::vector<char> get_temporal_control(int indx);
	void add_simulation_output_controls(std::vector<char> str){ simulation_output_controls.push_back(str); }
	void add_node_output_headers(std::vector<char> str){ node_output_headers.push_back(str); }
	void add_element_output_headers(std::vector<char> str){ element_output_headers.push_back(str); }
	std::vector<char> get_simulation_output_controls(int indx) const;
	std::vector<char> get_node_output_controls(int indx) const;
	std::vector<char> get_element_output_controls(int indx) const;
	void set_element_output_every(int val){ element_output_every = val; }
	void set_node_output_every(int val){ node_output_every = val; }
	void set_simulation_output_every(int val){ simulation_output_every = val; }
	void set_adsorption_string(std::vector<char> str){ adsorption_string = str; }
	std::vector<char> get_adsorption_string() const;
	int get_element_output_every() const;
	int get_node_output_every() const;
	int get_simulation_output_every() const;
	int get_observation_first_line() const;
	void set_observation_first_line(int val){ observation_first_line = val; }
	void add_obsData(std::vector<std::vector<char>> &obs){ obsData.push_back(obs); }
	void add_nodeData(std::vector<char> &node){ nodeData.push_back(node); }
	void add_elementData(std::vector<char> &element){ elementData.push_back(element); }
	void set_nbcfpr(int val){ NBCFPR = val; }
	void set_nbcspr(int val){ NBCSPR = val; }
	void set_nbcppr(int val){ NBCPPR = val; }
	void set_nbcupr(int val){ NBCUPR = val; }
	void set_cinact(char c){ CINACT = c; }
	int get_nbcfpr() const;
	int get_nbcspr() const;
	int get_nbcppr() const;
	int get_nbcupr() const;
	char get_cinact() const;
	void set_compfl(double compfl){COMPFL = compfl;}
	void set_cw(double cw){CW = cw;}
	void set_sigmaw(double sigmaw){SIGMAW = sigmaw;}
	void set_rhow0(double rhow0){RHOW0 = rhow0;}
	void set_urhow0(double urhow0){URHOW0 = urhow0;}
	void set_visc0(double visc0){VISC0 = visc0;}
	void set_drwdu(double drwdu){DRWDU = drwdu;}
	void set_compma(double compma){ COMPMA = compma; }
	void set_cs(double cs){ CS = cs; }
	void set_rhos(double rhos){ RHOS = rhos; }
	void set_sigmas(double sigmas){ SIGMAS = sigmas; }
	void set_prodf0(double prodf0){ PRODF0 = prodf0; }
	void set_prods0(double prods0){ PRODS0 = prods0; }
	void set_prodf1(double prodf1){ PRODF1 = prodf1; }
	void set_prods1(double prods1){ PRODS1 = prods1; }
	void set_gravx(double gravx){ GRAVX = gravx; }
	void set_gravy(double gravy){ GRAVY = gravy; }
	void set_gravz(double gravz){ GRAVZ = gravz; }
	void set_scalx(double scalx){ SCALX = scalx; }
	void set_scaly(double scaly){ SCALX = scaly; }
	void set_scalz(double scalz){ SCALX = scalz; }
	void set_porfac(double porfac){ PORFAC=porfac; }
	double get_compfl() const;
	double get_cw() const;
	double get_sigmaw() const;
	double get_rhow0() const;
	double get_urhow0() const;
	double get_visc0() const;
	double get_drwdu() const;
	double get_compma() const;
	double get_cs() const;
	double get_rhos() const;
	double get_sigmas() const;
	double get_prodf0() const;
	double get_prods0() const;
	double get_prodf1() const;
	double get_prods1() const;
	double get_gravx() const;
	double get_gravy() const;
	double get_gravz() const;
	double get_scalx() const;
	double get_scaly() const;
	double get_scalz() const;
	double get_porfac() const;
	void set_element_props(std::vector<std::vector<char>> str){ element_props = str; }
	std::vector<std::vector<char>> get_element_props() const;
	void reserveElementData(){ elementData.reserve(NE); }
	void reserve_npbc_data(){ npbcData.reserve(NPBC); }
	void reserve_nubc_data(){ nubcData.reserve(NUBC); }
	void reserve_nsop_data(){ nsopData.reserve(NSOP); }
	void reserve_nsou_data(){ nsouData.reserve(NSOU); }
	void reserve_incidence_data(){ incidenceData.reserve(NE+1); }
	void add_incidence_data(std::vector<char> &str){ incidenceData.push_back(str); }
	void add_nsop_data(std::vector<char> &str){ nsopData.push_back(str); }
	void add_nsou_data(std::vector<char> &str){ nsouData.push_back(str); }
	void add_npbc_data(std::vector<char> &str){ npbcData.push_back(str); }
	void add_nubc_data(std::vector<char> &str){ nubcData.push_back(str); }






private:
	static Storage * m_pInstance;
	std::unordered_map<std::string, DataSet *> dataSetMap;
	std::vector<std::string> titles;
	std::vector<char> sutra_string;
	std::vector<char> version_string;
	std::vector<char> version_num_string;
	std::vector<char> simulation_type_string;
	std::vector<char> transport_string;
	std::vector<char> mesh_dim_string;
	std::vector<char> mesh_type_string;
	std::vector<char> nn1_string;
	std::vector<char> nn2_string;
	std::vector<char> nn3_string;
	std::vector<char> flow_type_string;
	std::vector<char> transport_type_string;
	std::vector<char> simulation_condition_string;
	std::vector<char> simulation_start_string;
	std::vector<std::vector<char>> simulation_output_controls;
	std::vector<std::vector<char>> node_output_headers;
	std::vector<std::vector<char>> element_output_headers;
	std::vector<std::vector<char>> element_props;
	std::vector<char> p_solver_string;
	std::vector<char> u_solver_string;
	std::vector<char> adsorption_string;
	std::vector<std::vector<char>> temporalControl;
	std::vector<std::vector<std::vector<char>>> obsData;
	std::vector<std::vector<char>> nodeData;
	std::vector<std::vector<char>> elementData;
	std::vector<std::vector<char>> incidenceData;
	std::vector<std::vector<char>> nsopData;
	std::vector<std::vector<char>> nsouData;
	std::vector<std::vector<char>> npbcData;
	std::vector<std::vector<char>> nubcData;
	int solution_storage;
	int NN, NE, NSOP, NSOU, NPBC, NUBC,NOBS;
	int NSCH, NPCYC, NUCYC;
	int ITRMAX, RPMAX, RUMAX;
	int max_p_iterations;
	int max_u_iterations;
	int simulation_output_every;
	int node_output_every;
	int element_output_every;
	int observation_first_line;
	int NBCFPR, NBCSPR, NBCPPR, NBCUPR;
	char CINACT;
	double u_tolerance;
	double p_tolerance;
	double UP, GNUP, GNUU;
	double COMPFL;
	double CW;
	double SIGMAW;

	double RHOW0;
	double URHOW0;
	double VISC0;
	double DRWDU;

	double COMPMA;
	double CS;
	double RHOS;
	double SIGMAS;
	double PRODF0;
	double PRODS0;
	double PRODF1;
	double PRODS1;

	double GRAVX;
	double GRAVY;
	double GRAVZ;

	double SCALX;
	double SCALY;
	double SCALZ;
	double PORFAC;
	Storage();
	~Storage();
};
#endif


