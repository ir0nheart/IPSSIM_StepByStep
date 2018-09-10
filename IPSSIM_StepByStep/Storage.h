#ifndef STORAGE_H
#define STORAGE_H
#pragma once
#include "DataSet.h"
#include "Schedule.h"
#include "Node.h"
#include "Element.h"
#include "obsPoint.h"
#include "Bcs.h"

class Storage
{
public:
	static Storage * instance();
	void BASIS3_Simple(int L, double XLOC, double YLOC, double ZLOC, double& DET, double CJ[]);

	void BUBSAT(double& SWB, double& RELKB, double PRES, double CNUB, double & RELKT, double& SWT, double SW, double RELK);
	void DISPR3(double vx, double vy, double vz, double vmag, double ang1, double ang2, double ang3,
		double ALMAX, double ALMID, double ALMIN, double ATMAX, double ATMID, double ATMIN,
		double& dxx, double & dxy, double& dxz, double & dyx, double & dyy, double & dyz, double& dzx, double & dzy, double & dzz);
	void GLOCOL(int L,double vole[],double bflowe[8][8],double dflowe[],double btrane[8][8],double dtrane[8][8]);
	void BCTIME();
	void BCSTEP();
	void ADSORB();
	void ELEMN3();
	void NODAL();
	void BC();
	void ROTATE(std::vector<double>& vec, double& v1, double& v2, double& v3, std::vector<double>& out_vec);
	void BASIS3(int ICALL, int L, double XLOC, double YLOC, double ZLOC, double F[],double W[], double& DET,double CJ[],
		double DFDXG[],double DFDYG[],double DFDZG[],double DWDXG[],double DWDYG[],double DWDZG[],double & swbg,double& relkbg,
		double & vxg,double&vyg,double & vzg,double & vgmag,double& swtg,double&relktg,double &viscg,double & rhog,double &rgxg,double&rgyg,double&rgzg,double &porg);
	void ELEMN2();
	void set_bcs_defined(bool v){ bcs_defined = v; }
	bool get_bcs_defined()const { return bcs_defined; }
	void simulation();
	void output_initial_starting_if_transient();
	void determine_tmax();
	void check_restart();
	void addTittle(std::string str);
	void check_data_sets();
	void set_steady_state_switches();
	void init_a_val(double * vec, int ssize, double val){ for (int i = 0; i < ssize; i++) vec[i] = val; }
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
	void add_layerData(std::vector<char> &layer){ layerData.push_back(layer); }
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
	void set_scaly(double scaly){ SCALY = scaly; }
	void set_scalz(double scalz){ SCALZ = scalz; }
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
	void set_flags();
	void set_nn1(int val){ NN1 = val; }
	void set_nn2(int val){ NN2 = val; }
	void set_nn3(int val){ NN3 = val; }
	int get_nn1() const;
	int get_nn2() const;
	int get_nn3() const;
	void set_pstar(double val){ PSTAR = val; }
	void set_gconst(double val){GCONST = val; }
	void set_temp(double val){ TEMP= val; }
	void set_smwh(double val){ SMWH = val; }
	void set_water_table(double val){ water_table = val; }
	void set_time_step_divide(double val){ time_step_divide = val; }
	void set_number_of_layers(int val){ number_of_layers = val; }
	void add_p_ics(double val){ p_ics.push_back(val); }
	void add_u_ics(double val){ u_ics.push_back(val); }
	void set_t_ics(double val){ t_ics = val; }
	void set_p_ics_string(std::vector<char> str){ p_ics_string = str; }
	void set_u_ics_string(std::vector<char> str){ u_ics_string = str; }
	void reserve_p_ics(){ p_ics.reserve(NN); };
	void reserve_u_ics(){ u_ics.reserve(NN); };
	int FRCSTP(double time);
	int FINDL3(int el_no,obsPoint& obs,double& xsi_,double& eta_,double& zet_);
	void set_starting_time();
	void PTRSET();
	void BANWID();
	void allocate_element_arrays();
	void de_allocate_element_arrays();
	void allocate_node_arrays();
	void de_allocate_node_arrays();
	void ROTMAT(double& a1, double& a2, double& a3, std::vector<double>& vec);
	void TENSYM(double& pmax, double& pmid, double& pmin, std::vector<double>& rotMat, double& permxx, double& permxy, double &permxz, double& permyx, double&permyy, double& permyz, double& permzx, double& permzy, double& permzz);
	std::vector<Schedule *> get_schedule_list(){ return schedule_list; }
	void add_bcs(Bcs * bcs){ bcsContainer.push_back(bcs); }
	std::vector<Bcs *> get_bcs_container(){ return bcsContainer; }
	double DNRM2(int N, double * X, int INCX);

	void re_orient_matrix(int jmper_size, int vals_size, double vals[], std::vector<int>&jmper, std::vector<int>& indices, double * new_vals, int * new_jmper, int * new_indices);
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
	std::vector<std::vector<char>> layerData;
	std::vector<std::vector<char>> nsopData;
	std::vector<std::vector<char>> nsouData;
	std::vector<std::vector<char>> npbcData;
	std::vector<std::vector<char>> nubcData;

	std::vector<char> p_ics_string;
	std::vector<char> u_ics_string;

	std::vector<int> IQSOP;
	std::vector<int> IQSOU;
	std::vector<int>IPBC;
	std::vector<int>IUBC;
	std::vector<int> incidence_vector;
	std::vector<int> IA;
	std::vector<int> JA;

	std::vector<double> p_ics;
	std::vector<double> u_ics;

	std::vector<Schedule *> schedule_list;
	bool ONCEP;
	bool SETBCS;
	bool INTIM;
	char CINACT;

	int IBCT;
	int solution_storage;
	int NN;
	int NE, NSOP, NSOU, NPBC, NUBC, NOBS;
	int NSCH, NPCYC, NUCYC;
	int NN1, NN2, NN3;
	int ITRMAX, RPMAX, RUMAX;
	int max_p_iterations;
	int max_u_iterations;
	int simulation_output_every;
	int node_output_every;
	int element_output_every;
	int observation_first_line;
	int NBCFPR, NBCSPR, NBCPPR, NBCUPR;
	int number_of_layers;
	int NOBCYC;
	
	
	int KSOLVP;
	int ISTOP;
	int ME; // -1 For Solute , + 1 for ENERGY
	int IUNSAT;
	int ISSFLO;
	int ISSTRA;
	int IREAD;
	int NSCHAU;
	int time_steps_index;
	int NSAVEP;
	int NSAVEU;
	int NBCN;
	int NOBSN;
	int N48;
	int NEX;
	int NIN;
	int KNODAL;
	int KELMNT;
	int KINCID;
	int KPANDS;
	int KVEL;
	int KCORT;
	int KBUDG;
	int KSCRN;
	int KPAUSE;
	int NSOPI;
	int NSOUI;
	int IQSOPT;
	int IQSOUT;
	int IPBCT;
	int IUBCT;
	int ITMAX;
	int NDIMJA;
	int NELT;
	int NB;
	int NBHALF;
	int NBI;
	int max_bandwidth_element;
	int ITRST;
	int IT;
	int ITBCS;
	int ITREL;
	int ML;
	int NOUMAT;
	int ITER;

	double TSECM1;
	double RELCHG;
	double DELTM1;
	double DLTPM1;
	double DLTUM1;
	double BDELP1;
	double BDELU1;
	double BDELP;
	double BDELU;
	double TELAPS;
	double TSECP0;
	double TSECU0;
	double TMIN;
	double THOUR;
	double TDAY;
	double TWEEK;
	double TMONTH;
	double TYEAR;
	double DIT;
	double DELTLC;
	double TMAX;
	double TEMAX;
	double DELTP;
	double DELTU;
	double TSEC;
	double TSTART;
	double DELT;
	double CHI1;
	double CHI2;
	double PMAXFA;
	double PMIDFA;
	double PMINFA;
	double ANG1FA;
	double ANG2FA;
	double ANG3FA;
	double ALMAXF;
	double ALMIDF;
	double ALMINF;
	double ATMXF;
	double ATMDF;
	double ATMNF;
	
	//double GNUP1;
	//double GNUU1;
	double SCALX;
	double SCALY;
	double SCALZ;
	double PORFAC;
	double PSTAR;
	double GCONST;
	double TEMP;
	double SMWH;
	double time_step_divide;
	double water_table;
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
	double t_ics;

	std::vector<bool> BCSFL;
	std::vector<bool> BCSTR;
	std::vector<int> KTYPE; // KTYPE[0] for MESH TYPE 2D = 2, 3D=3  --- KTYPE[1] for MESH TYPE IRREGULAR=0,LAYERED=1,REGULAR=2,BLOCKWISE =3
	std::vector<Node> nodeContainer;
	std::vector<Element> elementContainer;

	std::vector<obsPoint> obsContainer;
	std::vector < std::pair<int, std::vector<int>>> incidenceContainer;
	// Element Relevant Data
	int * el_num;
	int * el_lreg;
	double * el_pmax;
	double * el_pmid;
	double * el_pmin;
	double * el_ang1;
	double * el_ang2;
	double * el_ang3;
	double * el_almax;
	double * el_almid;
	double * el_almin;
	double * el_atmax;
	double * el_atmid;
	double * el_atmin;
	double * el_pangl1;
	double * el_pangl2;
	double * el_pangl3;
	double * el_permxx;
	double * el_permxy;
	double * el_permxz;
	double * el_permyx;
	double * el_permyy;
	double * el_permyz;
	double * el_permzx;
	double * el_permzy;
	double * el_permzz;
	double * el_vmag;
	double * el_vang1;
	double * el_vang2;

	std::vector<std::vector<double>> el_gxsi;
	std::vector<std::vector<double>> el_geta;
	std::vector<std::vector<double>> el_gzet;
	std::vector<std::vector<double>> el_det;

	// Node Relevant Data

	int * node_num;
	int * node_nreg;
	double * node_x;
	double * node_y;
	double * node_z;
	double * node_por;
	double * node_sop;
	double * node_qin;
	double * node_uin;
	double * node_quin;
	double * node_pbc;
	double * node_ubc;
	double * node_pvec;
	double * node_uvec;
	double * node_pm1;
	double * node_um1;
	double * node_um2;
	double * node_rcit;
	double * node_sw;
	double * node_swt;
	double * node_swb;
	double * node_cnub;
	double * node_dswdp;
	double * node_cs1;
	double * node_cs2;
	double * node_cs3;
	double * node_sl;
	double * node_sr;
	double * node_dpdtitr;
	double * node_piter;
	double * node_pvel;
	double * node_uiter;
	double * node_rcitm1;
	double * node_qinitr;
	double * node_cnubm1;
	double * node_vol;
	double * node_p_rhs;
	double * node_u_rhs;
	double * node_p_solution;
	double * node_u_solution;
	double * node_rho;
	double * node_relk;
	double * node_relkb;
	double * node_relkt;
	double * PMAT;
	double * UMAT;

	int * row_jumper;
	int * col_indices;
	double * new_MAT;

	int * IBCPBC;
	int * IBCUBC;
	int * IBCSOP;
	int * IBCSOU;
	std::vector<std::vector<int>> node_neighbors;
	std::vector<double> QPLITR;
	std::vector<double> GNUP1;
	std::vector<double> GNUU1;
	std::vector<Bcs *> bcsContainer;
	std::vector<double> GXLOC;
	std::vector<double> GYLOC;
	std::vector<double> GZLOC;
	bool switch_set;
	bool bcs_defined;
	Storage();
	~Storage();
};
#endif


