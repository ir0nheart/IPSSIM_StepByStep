#include "stdafx.h"
#include "Storage.h"
#include "SimulationControl.h"
#include "obsPoint.h"
#include "Timer.h"
#include <algorithm>
#include "Parser.h"
#include "InputFiles.h"
#include "Writer.h"
Storage * Storage::m_pInstance = nullptr;
std::string Storage::K5SYM[] = { "N", "X", "Y", "Z", "P", "U", "S", "ES", "RU" };
std::string Storage::K6SYM[] = { "E","X","Y","Z","VX","VY","VZ" };
std::string Storage::VARNK5[] = { "NODE NUMBER", "X-COORDINATE", "Y-COORDINATE", "Z-COORDINATE", "PRESSURE", "CONCENTRATION", "SATURATION", "EFFECTIVE STRESS", "STRESS RATIO" };
std::string Storage::VARNK6[] = { "ELEMENT NUMBER", "X-COORDINATE OF CENTROID", "Y-COORDINATE OF CENTROID", "Z-COORDINATE OF CENTROID", "X-VELOCITY", "Y-VELOCITY", "Z-VELOCITY" };
int Storage::J6COL[] = {0,0,0,0,0,0,0};
std::string Storage::LCOL[] = { "X origin","Y origin","Z origin","X velocity", "Y velocity","Z velocity" };
int Storage::J5COL[] = {0,0,0,0,0,0,0,0,0};
//std::string Storage::NCOL[] = { "Node","X","Y","Z","Pressure","Concentration","Saturation","Eff. Stress","Stress Rat." };
std::string Storage::NCOL[] = { "X", "Y", "Z", "Pressure", "Concentration", "Saturation" };
template<typename MatrixT, typename VectorT>
struct monitor_user_data
{
	monitor_user_data(MatrixT const & A, VectorT const & b, VectorT const & guess) : A_ptr(&A), b_ptr(&b), guess_ptr(&guess) {}

	MatrixT const *A_ptr;
	VectorT const *b_ptr;
	VectorT const *guess_ptr;
};

/**
*  The actual callback-routine takes the current approximation to the result as the first parameter, and the current estimate for the relative residual norm as second argument.
*  The third argument is a pointer to our user-data, which in a first step cast to the correct type.
*  If the monitor returns true, the iterative solver stops. This is handy for defining custom termination criteria, e.g. one-norms for the result change.
*  Since we do not want to terminate the iterative solver with a custom criterion here, we always return 'false' at the end of the function.
*
*  Note to type-safety evangelists: This void*-interface is designed to be later exposed through a shared library ('libviennacl').
*  Thus, user types may not be known at the point of compilation, requiring a void*-approach.
**/
template<typename VectorT, typename NumericT, typename MatrixT>
bool my_custom_monitor(VectorT const & current_approx, NumericT residual_estimate, void *user_data)
{
	// Extract residual:
	monitor_user_data<MatrixT, VectorT> const *data = reinterpret_cast<monitor_user_data<MatrixT, VectorT> const*>(user_data);

	// Form residual r = b - A*x, taking an initial guess into account: r = b - A * (current_approx + x_initial)
	VectorT x = current_approx + *data->guess_ptr;
	VectorT residual = *data->b_ptr - viennacl::linalg::prod(*data->A_ptr, x);
	VectorT initial_residual = *data->b_ptr - viennacl::linalg::prod(*data->A_ptr, *data->guess_ptr);

	std::cout << "Residual estimate vs. true residual: " << residual_estimate << " vs. " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(initial_residual) << std::endl;

	return false; // no termination of iteration
}

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

double Storage::get_rpmax() const
{
	return RPMAX;
}

double Storage::get_rumax() const
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
	INTIM = true;
	ISTOP = 0;
	KSOLVP = 1;
	KSOLVU = 1;


	GXLOC.assign({ -1, 1, 1, -1, -1, 1, 1, -1 });
	GYLOC.assign({ -1, -1, 1, 1, -1, -1, 1, 1 });
	GZLOC.assign({ -1, -1, -1, -1, 1, 1, 1, 1 });

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

int Storage::get_nn1() const
{
	return NN1;
}


int Storage::get_nn2() const
{
	return NN2;
}

int Storage::get_nn3() const
{
	return NN3;
}

int Storage::FRCSTP(double time)
{
	int val = -1;
	for (int i = 0; i < schedule_list[time_steps_index]->get_step_list_size(); i++)
	{
		if (time == schedule_list[time_steps_index]->get_step_time()[i].second){
			val = schedule_list[time_steps_index]->get_step_time()[i].first;
			break;
		}
	}
	return val;
}

int Storage::FINDL3(int el_no,obsPoint& obs,double& xsi_,double& eta_,double& zet_)
{
	int INOUT = 0; // -1 out , 0 == failed to converge, 1 = in
	int itrmax_ = 25;
	double tol_ = 0.001;
	double epsilon_ = 0.001;
	double ope = 1.0 + epsilon_;

	// element corner coordinates
	double x1, x2, x3, x4, x5, x6, x7, x8;
	double y1, y2, y3, y4, y5, y6, y7, y8;
	double z1, z2, z3, z4, z5, z6, z7, z8;
	// coefficients
	double ax, bx, cx, dx, ex, fx, gx, hx;
	double ay, by, cy, dy, ey, fy, gy, hy;
	double az, bz, cz, dz, ez, fz, gz, hz;
//	double xsi_, zet_,eta_;

	double f10, f20, f30, fp11, fp12, fp13, fp21, fp22, fp23, fp31, fp32, fp33;
	double s11, s12, s13, cf12, cf34, cf43, cf56;
	double detxsi, deteta, detzet,determ;
	double delxsi, deleta, delzet;
	int m0 = (el_no - 1) * 8;
	x1 = node_x[incidence_vector[m0]-1];
	x2 = node_x[incidence_vector[m0+1]-1];
	x3 = node_x[incidence_vector[m0+2]-1];
	x4 = node_x[incidence_vector[m0+3]-1];
	x5 = node_x[incidence_vector[m0+4]-1];
	x6 = node_x[incidence_vector[m0+5]-1];
	x7 = node_x[incidence_vector[m0+6]-1];
	x8 = node_x[incidence_vector[m0+7]-1];
	y1 = node_y[incidence_vector[m0]-1];
	y2 = node_y[incidence_vector[m0 + 1]-1];
	y3 = node_y[incidence_vector[m0 + 2]-1];
	y4 = node_y[incidence_vector[m0 + 3]-1];
	y5 = node_y[incidence_vector[m0 + 4]-1];
	y6 = node_y[incidence_vector[m0 + 5]-1];
	y7 = node_y[incidence_vector[m0 + 6]-1];
	y8 = node_y[incidence_vector[m0 + 7]-1];
	z1 = node_z[incidence_vector[m0]-1];
	z2 = node_z[incidence_vector[m0 + 1]-1];
	z3 = node_z[incidence_vector[m0 + 2]-1];
	z4 = node_z[incidence_vector[m0 + 3]-1];
	z5 = node_z[incidence_vector[m0 + 4]-1];
	z6 = node_z[incidence_vector[m0 + 5]-1];
	z7 = node_z[incidence_vector[m0 + 6]-1];
	z8 = node_z[incidence_vector[m0 + 7]-1];


	ax = +x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;
	bx = -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8;
	cx = -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8;
	dx = -x1 - x2 - x3 - x4 + x5 + x6 + x7 + x8;
	ex = +x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8;
	fx = +x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8;
	gx = +x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8;
	hx = -x1 + x2 - x3 + x4 + x5 - x6 + x7 - x8;

	ay = +y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8;
	by = -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8;
	cy = -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8;
	dy = -y1 - y2 - y3 - y4 + y5 + y6 + y7 + y8;
	ey = +y1 - y2 + y3 - y4 + y5 - y6 + y7 - y8;
	fy = +y1 - y2 - y3 + y4 - y5 + y6 + y7 - y8;
	gy = +y1 + y2 - y3 - y4 - y5 - y6 + y7 + y8;
	hy = -y1 + y2 - y3 + y4 + y5 - y6 + y7 - y8;

	az = +z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8;
	bz = -z1 + z2 + z3 - z4 - z5 + z6 + z7 - z8;
	cz = -z1 - z2 + z3 + z4 - z5 - z6 + z7 + z8;
	dz = -z1 - z2 - z3 - z4 + z5 + z6 + z7 + z8;
	ez = +z1 - z2 + z3 - z4 + z5 - z6 + z7 - z8;
	fz = +z1 - z2 - z3 + z4 - z5 + z6 + z7 - z8;
	gz = +z1 + z2 - z3 - z4 - z5 - z6 + z7 + z8;
	hz = -z1 + z2 - z3 + z4 + z5 - z6 + z7 - z8;

	xsi_ = eta_ = zet_ = 0;

	for (int i = 0; i < itrmax_; i++)
	{
		f10 = ax - 8 * obs.get_x() + bx*xsi_ + cx*eta_ + dx*zet_ + ex*xsi_*eta_
			+ fx*xsi_*zet_ + gx*eta_*zet_ + hx * xsi_ * eta_*zet_;
		f20 = ay - 8 * obs.get_y() + by*xsi_ + cy*eta_ + dy*zet_ + ey*xsi_*eta_
			+ fy*xsi_*zet_ + gy*eta_*zet_ + hy* xsi_ * eta_*zet_;
		f30 = az - 8 * obs.get_z() + bz*xsi_ + cz*eta_ + dz*zet_ + ez*xsi_*eta_
			+ fz*xsi_*zet_ + gz*eta_*zet_ + hz* xsi_ * eta_*zet_;

		fp11 = bx + ex*eta_ + fx*zet_ + hx*eta_*zet_;
		fp12 = cx + ex*xsi_ + gx*zet_ + hx*xsi_*zet_;
		fp13 = dx + fx*xsi_ + gx*eta_ + hx*xsi_*eta_;

		fp21 = by + ey*eta_ + fy*zet_ + hy*eta_*zet_;
		fp22 = cy + ey*xsi_ + gy*zet_ + hy*xsi_*zet_;
		fp23 = dy + fy*xsi_ + gy*eta_ + hy*xsi_*eta_;

		fp31 = bz + ez*eta_ + fz*zet_ + hz*eta_*zet_;
		fp32 = cz + ez*xsi_ + gz*zet_ + hz*xsi_*zet_;
		fp33 = dz + fz*xsi_ + gz*eta_ + hz*xsi_*eta_;

		s11 = fp22*fp33 - fp32*fp23;
		s12 = fp21*fp33 - fp31*fp23;
		s13 = fp21*fp32 - fp31*fp22;
		cf12 = -f20*fp33 + f30*fp23;
		cf34 = -f20*fp32 + f30*fp22;
		cf43 = -cf34;
		cf56 = -f30*fp21 + f20*fp31;

		detxsi = -f10*s11 - fp12*cf12 + fp13*cf34;
		deteta = fp11*cf12 + f10*s12 + fp13*cf56;
		detzet = fp11*cf43 - fp12*cf56 - f10*s13;
		determ = fp11*s11 - fp12*s12 + fp13*s13;
		delxsi = detxsi / determ;
		deleta = deteta / determ;
		delzet = detzet / determ;

		xsi_ += delxsi;
		eta_ += deleta;
		zet_ += delzet;


		if (abs(delxsi) < tol_ && abs(deleta) < tol_ && abs(delzet) <tol_)
		{
			INOUT = 1;
			if (abs(xsi_) > ope || abs(eta_) > ope || abs(zet_) > ope)
				INOUT = -1;
			break;
		}
	}
	return INOUT;
}

void Storage::PTRSET()
{
	for (int j = 0; j < incidenceContainer.size(); j++)
	{
		int nod_, neighbor_;
		for (int k = 0; k < N48; k++)
		{
			nod_ = incidenceContainer[j].second[k];
			for (int t = 0; t < N48; t++)
			{
				neighbor_ = incidenceContainer[j].second[t];
				if (neighbor_ != nod_)
				{
				//	nodeContainer[nod_ - 1].add_neighbor(neighbor_);
					node_neighbors[nod_ - 1].push_back(neighbor_);
				}
			}
		}
		
	}
	NELT = 0;
	for (int i = 0; i < NN; i++)
	{
		//nodeContainer[i].make_neighbors_unique();
		sort(node_neighbors[i].begin(), node_neighbors[i].end());
		node_neighbors[i].erase(std::unique(node_neighbors[i].begin(), node_neighbors[i].end()), node_neighbors[i].end());
		//NELT += nodeContainer[i].get_neighbor_size();
		NELT += node_neighbors[i].size();
	}
	NELT = NELT + NN;

	IA.reserve(NELT);
	JA.reserve(NDIMJA);

	for (int i = 0; i < NN; i++)
	{
		JA.push_back(IA.size());
		IA.push_back(i); // it was i+1
		for (int k = 0; k < node_neighbors[i].size(); k++)
		{
			IA.push_back(node_neighbors[i][k] -1);
		}
	}
	JA.push_back(NELT);

	// ALLOCATE PMAT UMAT
	PMAT = new double[NELT];
	UMAT = new double[NELT];
	new_MAT = new double[NELT]{};
	row_jumper = new int[NN + 1]{};
	col_indices = new int[NELT]{};
}

void Storage::BANWID()
{
	// find element with maximum difference in node numbers
	int ielo;
	int iehi;
	int ndif = 0;
	int ndiff;
	int element;
	for (int i = 0; i < incidenceContainer.size(); i++)
	{
		ielo = incidenceContainer[i].second[0];
		iehi = incidenceContainer[i].second[0];
		for (int j = 1; j < N48; j++)
		{
			if (incidenceContainer[i].second[j] < ielo) ielo = incidenceContainer[i].second[j];
			if (incidenceContainer[i].second[j] > iehi) iehi = incidenceContainer[i].second[j];
		}
		ndiff = iehi - ielo;
		if (ndiff > ndif)
		{
			ndif = ndiff;
			max_bandwidth_element = incidenceContainer[i].first;
		}
	}
	NB = 2 * ndif + 1;
	NBHALF = ndif + 1;
	NBI = NB;

}


void Storage::check_data_sets()
{
	Writer * lstWriter = Writer::instance("LST");
	std::string logLine;
	// Check Titles
	if (titles.size() != 2)
	{
		std::cout << "Titles are not properly defined." << std::endl;
		SimulationControl::exitOnError();
	}


	if (!strcmp(sutra_string.data(), "SUTRA"))
	{
		SimulationControl::exitOnError("INP-2A-1");
	}

	if (!strncmp(version_string.data(), "VERSION", 7))
	{
		if (version_num_string.data() == "2D3D.1")
		{
			version_num_string = std::vector<char>{'2', '.', '0'};
		}
		else if (!strcmp(version_num_string.data(), "2.0") &&
			!strcmp(version_num_string.data(), "2.1") &&
			!strcmp(version_num_string.data(), "2.2"))
		{
			SimulationControl::exitOnError("INP-2A-4");
		}
	}
	else
	{
		version_num_string = std::vector<char>{'2', '.', '0'};
	}

	if (strcmp(simulation_type_string.data(), "SOLUTE"))
	{
		ME = -1;
		logLine.append("\n\n");
		logLine.append(std::string(132, '*'));
		logLine.append("\n\n\n");
		logLine.append(std::string(20, ' '));
		logLine.append("* * * * *   I P S S I M   S O L U T E   T R A N S P O R T   S I M U L A T I O N   * * * * *\n\n\n");
		logLine.append(std::string(132, '*'));
		lstWriter->add_line(logLine);
		logLine.clear();
	}
	else if (strcmp(simulation_type_string.data(), "ENERGY"))
	{
		ME = +1;
		//std::cout << "ENERGY TRANSPORT " << std::endl;
		logLine.append("\n\n");
		logLine.append(std::string(132, '*'));
		logLine.append("\n\n\n");
		logLine.append(std::string(20, ' '));
		logLine.append("* * * * *   I P S S I M   E N E R G Y   T R A N S P O R T   S I M U L A T I O N   * * * * *\n\n\n");
		logLine.append(std::string(132, '*'));
		lstWriter->add_line(logLine);
		logLine.clear();
	}
	else
	{
		//std::cout << "Simulation TYPE is not Solute or Energy Transport " << std::endl;
		SimulationControl::exitOnError("INP-2A-2");
	}

	// Output Titles

	logLine.append("\n\n\n\n ");
	logLine.append(std::string(131, '-'));
	logLine.append("\n\n");
	logLine.append(std::string(26, ' '));
	logLine.append(titles[0]);
	logLine.append("\n\n");
	logLine.append(std::string(26, ' '));
	logLine.append(titles[1]);
	logLine.append("\n\n ");
	logLine.append(std::string(131, '-'));
	lstWriter->add_line(logLine);
	logLine.clear();

	//Write File Assignments in LST file
	InputFiles::instance()->printInputFilesToLST();

	if (!strncmp(mesh_dim_string.data(), "2D", 2))
	{
		KTYPE.push_back(2);
	}
	else if (!strncmp(mesh_dim_string.data(), "3D", 2))
	{
		KTYPE.push_back(3);
	}
	else
	{
		std::cout << "Error in Mesh type definition." << std::endl;
		std::cout << "\t Mesh type should be either '2D' or '3D'" << std::endl;
		SimulationControl::exitOnError("INP-2B-1");
	}

	// This program runs only REGULAR and IRREGULAR MESH
	if (!strncmp(mesh_type_string.data(), "REGULAR", mesh_type_string.size()))
	{
		KTYPE.push_back(2);
		if (KTYPE[0] == 2)
		{
			NN1 = std::stoi(std::string(nn1_string.begin(), nn1_string.end()));
			NN2 = std::stoi(std::string(nn2_string.begin(), nn2_string.end()));
			NN3 = 1;
		}
		else if (KTYPE[0] == 3)
		{
			NN1 = std::stoi(std::string(nn1_string.begin(), nn1_string.end()));
			NN2 = std::stoi(std::string(nn2_string.begin(), nn2_string.end()));
			NN3 = std::stoi(std::string(nn3_string.begin(), nn3_string.end()));
		}
		if (NN1 < 2 || NN2 < 2 || (KTYPE[0] == 3 && NN3 < 2))
		{
			SimulationControl::exitOnError("INP-2B-3");
		}

	}
	else if (!strncmp(mesh_type_string.data(), "IRREGULAR", mesh_type_string.size()))
	{
		KTYPE.push_back(0);
	}
	else
	{
		std::cout << "Mesh is not REGULAR or IRREGULAR" << std::endl;
		SimulationControl::exitOnError("INP-2B-4");
	}


	if (KTYPE[1] > 1)
	{
		int NN123 = NN1*NN2*NN3;
		int NE123 = 0;
		if (NN123 != NN)
		{
			SimulationControl::exitOnError("INP-2B,3-1");
		}
		if (KTYPE[0] == 3)
		{
			NE123 = (NN1 - 1)*(NN2 - 1)*(NN3 - 1);
		}
		else
		{
			NE123 = (NN1 - 1)*(NN2 - 1);
		}
		if (NE123 != NE)
		{
			SimulationControl::exitOnError("INP-2B,3-2");
		}
	}


	logLine.append("\n\n\n\n           S I M U L A T I O N   M O D E   O P T I O N S\n\n");

	if (!strncmp(simulation_condition_string.data(), "UNSATURATED", simulation_condition_string.size())){
		IUNSAT = +1;
		logLine.append("           - ALLOW UNSATURATED AND SATURATED FLOW: UNSATURATED PROPERTIES ARE USER-PROGRAMMED IN UNSAT \n");
	}
	else if (!strncmp(simulation_condition_string.data(), "SATURATED", simulation_condition_string.size())){
		IUNSAT = 0;
		logLine.append("           - ASSUME SATURATED FLOW ONLY \n");
	}
	else{
		SimulationControl::exitOnError("INP-4-1");
	}


	if (!strncmp(flow_type_string.data(), "TRANSIENT", flow_type_string.size())){
		ISSFLO = 0;
		logLine.append("           - ALLOW TIME-DEPENDENT FLOW FIELD \n");
	}
	else if (!strncmp(flow_type_string.data(), "STEADY", flow_type_string.size())){
		ISSFLO = +1;
		if (ME == -1)
			logLine.append("           - ASSUME STEADY-STATE FLOW FIELD CONSISTENT WITH INITIAL CONCENTRATION CONDITIONS \n");
		if (ME == +1)
			logLine.append("           - ASSUME STEADY-STATE FLOW FIELD CONSISTENT WITH INITIAL TEMPERATURE CONDITIONS \n");
	}
	else{
		SimulationControl::exitOnError("INP-4-2");
	}



	if (!strncmp(transport_type_string.data(), "TRANSIENT", transport_type_string.size())){
		ISSTRA = 0;
		logLine.append("           - ALLOW TIME-DEPENDENT TRANSPORT \n");
	}
	else if (!strncmp(transport_type_string.data(), "STEADY", transport_type_string.size())){
		ISSTRA = +1;
		logLine.append("           - ASSUME STEADY-STATE TRANSPORT \n");
	}
	else{
		SimulationControl::exitOnError("INP-4-3");
	}


	if (!strncmp(simulation_start_string.data(), "COLD", simulation_start_string.size())){
		IREAD = +1;
		logLine.append("           - COLD START - BEGIN NEW SIMULATION \n");
	}
	else if (!strncmp(simulation_start_string.data(), "WARM", simulation_start_string.size())){
		IREAD = -1;
		logLine.append("           - WARM START - SIMULATION IS TO BE CONTINUED FROM PREVIOUSLY-STORED DATA \n");
	}
	else{
		SimulationControl::exitOnError("INP-4-4");
	}

	if (solution_storage > 0)
	{
		logLine.append("           - STORE RESULTS AFTER EVERY ");
		logLine.append(std::to_string(solution_storage));
		logLine.append(" TIME STEPS IN RESTART FILE AS BACKUP AND FOR USE IN A SIMULATION RESTART\n");
	}
	else if (solution_storage == 0)
	{
		logLine.append("           - DO NOT STORE RESULTS FOR USE IN A RESTART OF SIMULATION \n");
	}
	
	lstWriter->add_line(logLine);
	logLine.clear();

	if (ME == -1)
	{
		char buff[1024];
		logLine.append("\n\n\n\n\n           S I M U L A T I O N   C O N T R O L   N U M B E R S\n\n");
		_snprintf(buff, 1024, "        %9d     NUMBER OF NODES IN FINITE-ELEMENT MESH\n", NN);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %9d     NUMBER OF ELEMENTS IN MESH\n\n", NN);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %9d     EXACT NUMBER OF NODES IN MESH AT WHICH PRESSURE IS A SPECIFIED CONSTANT OR FUNCTION OF TIME\n", NPBC);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %9d     EXACT NUMBER OF NODES IN MESH AT WHICH SOLUTE CONCENTRATION IS A SPECIFIED CONSTANT OR FUNCTION OF TIME\n", NUBC);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %9d     EXACT NUMBER OF NODES IN MESH AT WHICH FLUID INFLOW OR OUTFLOW IS A SPECIFIED CONSTANT OR FUNCTION OF TIME\n", NSOP);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %9d     EXACT NUMBER OF NODES IN MESH AT WHICH A SOURCE OR SINK OF SOLUTE MASS IS A SPECIFIED CONSTANT OR FUNCTION OF TIME\n\n", NSOU);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %9d     EXACT NUMBER OF NODES IN MESH AT WHICH PRESSURE AND CONCENTRATION WILL BE OBSERVED\n", NOBS);
		logLine.append(buff);
		
		logLine.append("\n\n\n\n\n           N U M E R I C A L   C O N T R O L   D A T A\n\n");
		_snprintf(buff, 1024, "        %+15.5f     'UPSTREAM WEIGHTING' FACTOR\n", UP);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %+15.4e     SPECIFIED PRESSURE BOUNDARY CONDITION FACTOR\n", GNUP);
		logLine.append(buff);
		_snprintf(buff, 1024, "        %+15.4e     SPECIFIED CONCENTRATION BOUNDARY CONDITION FACTOR\n", GNUU);
		logLine.append(buff);

		lstWriter->add_line(logLine);
		logLine.clear();
	}
	else
	{
		SimulationControl::exitOnError("ENERGY TRANSPORT IS NOT IMPLEMENTED IN IPSSIM. ");
	}

	if (ISSTRA == 1 && ISSFLO != 1)
	{
		SimulationControl::exitOnError("INP-4-5");
	}

	if (NPCYC < 1 || NUCYC <1)
	{
		SimulationControl::exitOnError("INP-6-1");
	} else if (NPCYC != 1 && NUCYC != 1)
	{
		SimulationControl::exitOnError("INP-6-2");
	}

	NSCHAU = 3;

	// Evde Bakacam
	logLine.append("\n\n\n\n\n           T E M P O R A L   C O N T R O L   A N D   S O L U T I O N   C Y C L I N G   D A T A\n\n");



	std::vector<std::pair<int, double>> step_time;
	if (ISSTRA == 1)
	{
		TSTART = t_ics;
		DELT = max(0.1*abs(TSTART), 1.0);
		int step=0;
		double time=TSTART;
		schedule_list.reserve(4);
		schedule_list.push_back(new Schedule("TIME_STEPS"));
		schedule_list.push_back(new Schedule("STEPS_1&UP"));
		schedule_list.push_back(new Schedule("STEP_0"));
		schedule_list.push_back(new Schedule("STEP_1"));

		schedule_list[0]->add_step_time(step, time);
		schedule_list[2]->add_step_time(step, time);
		time = time + DELT;
		step = 1;
		schedule_list[0]->add_step_time(step, time);
		schedule_list[1]->add_step_time(step, time);
		schedule_list[3]->add_step_time(step, time);

		ITMAX = 1;
		char buff[1024];
		logLine.append("\n             NOTE: BECAUSE FLOW AND TRANSPORT ARE STEADY-STATE, USER-DEFINED SCHEDULES ARE NOT IN EFFECT.");
		logLine.append("\n             STEADY-STATE RESULTS WILL BE WRITTEN TO THE APPROPRIATE OUTPUT FILES.");
		_snprintf(buff, 1024, "\n\n          THE FOLLOWING %d SCHEDULES CAN BE USED TO CONTROL SPECIFICATIONS OF STEADY-STATE BOUNDARY", NSCH);
		logLine.append(buff);
		logLine.append("\n             CONDITIONS IN (OPTIONAL) .BCS FILES:");
		logLine.append("\n                SCHEDULE TIME_STEPS   CONSISTS OF TIME STEPS 0 (STEADY FLOW) AND 1 (STEADY TRANSPORT);\n");
		logLine.append(std::string(41, ' '));
		logLine.append("THIS SCHEDULE IS DEFINED AUTOMATICALLY BY IPSSIM\n");
		lstWriter->add_line(logLine);
		goto _l846;
	}

	NSCH = NSCH + NSCHAU;
	char buff[1024];
	_snprintf(buff, 1024, "\n             THE %5d SCHEDULES ARE LISTED BELOW. SCHEDULE 'TIME_STEPS' CONTROLS TIME STEPPING.\n",NSCH);
	logLine.append(buff);
	schedule_list.reserve(NSCH);
	// READ SCHEDULE NAME AND TYPE
	// BASED ON TYPE CONSTRUCT SCHEDULE
	for (std::vector<char> sch : temporalControl)
	{
		std::string sch_nam;
		std::string sch_type;
		sch.push_back('\0');
		sch_nam = strtok(sch.data(), " ");
		sch_type = strtok(NULL, " ");
		sch_type.append(" ");
		sch_type.append(strtok(NULL, " "));
		schedule_list.push_back(new Schedule(sch_nam, sch_type));
	}

	for (int i = 0; i < schedule_list.size();i++)
	{
		if (schedule_list[i]->get_type() == "STEP CYCLE")
		{
			schedule_list[i]->set_sbased(false);
			std::vector<char> str = temporalControl[i];
			strtok(str.data(), " ");
			strtok(NULL, " ");
			strtok(NULL, " ");
			int NSMAX, ISTEPI, ISTEPL, ISTEPC;
			int NSTEP, NDSTEP;
			double TIME;
			int STEP;
			NSMAX = ISTEPI = ISTEPL = ISTEPC = 0;
			NSMAX = std::stoi(strtok(NULL, " "));
			ISTEPI = std::stoi(strtok(NULL, " "));
			ISTEPL = std::stoi(strtok(NULL, " "));
			ISTEPC = std::stoi(strtok(NULL, " "));

			NSTEP = ISTEPI;
			NDSTEP = ISTEPC;
			STEP = NSTEP;
			TIME = STEP;
			schedule_list[i]->add_step_time(STEP, TIME);
			for (int ii = 0; ii < NSMAX; ii++)
			{
				NSTEP = NSTEP + NDSTEP;
				STEP = NSTEP;
				TIME = STEP;
				schedule_list[i]->add_step_time(STEP, TIME);
			}

			logLine.append("\n");
			logLine.append(std::string(16, ' '));
			logLine.append("SCHEDULE ");
			logLine.append(schedule_list[i]->get_name());
			logLine.append("   STEP CYCLE WITH THE FOLLOWING SPECIFICATIONS:\n");
			logLine.append(std::string(40, ' '));
			_snprintf(buff, 1024, "%8d     MAXIMUM NUMBER OF TIME AFTER INITIAL TIME STEP NUMBER\n",NSMAX);
			logLine.append(buff);
			logLine.append(std::string(40, ' '));
			_snprintf(buff, 1024, "%8d     INITIAL TIME STEP NUMBER\n",ISTEPI);
			logLine.append(buff);
			logLine.append(std::string(40, ' '));
			_snprintf(buff, 1024, "%8d     LIMITING TIME STEP NUMBER\n", ISTEPL);
			logLine.append(buff);
			logLine.append(std::string(40, ' '));
			_snprintf(buff, 1024, "%8d     TIME STEP INCREMENT\n", ISTEPC);
			logLine.append(buff);
			lstWriter->add_line(logLine);
			logLine.clear();

		} else if (schedule_list[i]->get_type() == "TIME CYCLE")
		{
			schedule_list[i]->set_sbased(false);
			std::vector<char> str = temporalControl[i];
			str.push_back('\0');
			strtok(str.data(), " ");
			strtok(NULL, " ");
			strtok(NULL, " ");
			std::string CREFT;
			double SCALT;
			double TIMEI, TIMEL, TIMEC, TCMULT, TCMIN, TCMAX;
			int NTMAX, NTCYC;
			CREFT = strtok(NULL, " ");
			SCALT = std::stod(strtok(NULL, " "));
			NTMAX = std::stoi(strtok(NULL, " "));
			TIMEI = std::stod(strtok(NULL, " "));
			TIMEL = std::stod(strtok(NULL, " "));
			TIMEC = std::stod(strtok(NULL, " "));
			NTCYC = std::stoi(strtok(NULL, " "));
			TCMULT = std::stod(strtok(NULL, " "));
			TCMIN = std::stod(strtok(NULL, " "));
			TCMAX = std::stod(strtok(NULL, " "));

			if (CREFT == "ELAPSED")
			{
				if (schedule_list[i]->get_name() == "TIME_STEPS" && TIMEI != 0.0)
				{
					SimulationControl::exitOnError("INP-6-7");
				}
				schedule_list[i]->set_elapsed(true);

			} else if (CREFT == "ABSOLUTE")
			{
				schedule_list[i]->set_elapsed(false);
			} else
			{
				SimulationControl::exitOnError("INP-6-6");
			}

			TIMEI *= SCALT;
			TIMEL *= SCALT;
			TIMEC *= SCALT;
			TCMIN *= SCALT;
			TCMAX *= SCALT;

			double time, dtime;
			int step;
			time = TIMEI;
			step = time;
			dtime = TIMEC;
			schedule_list[i]->add_step_time(step, time);
			int ctr = 0;
			for (int j = 0; j < NTMAX; j++)
			{
				if (ctr == NTCYC && j > 1)
				{
					dtime = dtime*TCMULT;
					ctr = 0;
				}

				if (time_step_divide > 0 && j== time_step_divide-1)
				{
					dtime = TIMEC;
					ctr = 0;
				}

				if (dtime > TCMAX)
					dtime = TCMAX;

				if (dtime < TCMIN)
					dtime = TCMIN;

				time += dtime;
				ctr++;
				step = time;
				schedule_list[i]->add_step_time(step, time);

				if (time >= TIMEL)
					break;
			}
			logLine.append("\n");
			logLine.append(std::string(16, ' '));
			logLine.append("SCHEDULE ");
			logLine.append(schedule_list[i]->get_name());
			logLine.append("   TIME CYCLE WITH THE FOLLOWING SPECIFICATIONS IN TERMS OF ");
			logLine.append(CREFT);
			logLine.append(" TIMES : \n");
			logLine.append(std::string(46, ' '));
			_snprintf(buff, 1024, "%8d     MAXIMUM NUMBER OF TIMES AFTER INITIAL TIME\n", NTMAX);
			logLine.append(buff);
			logLine.append(std::string(39, ' '));
			_snprintf(buff, 1024, "%15.7e     INITIAL TIME\n", TIMEI);
			logLine.append(buff);
			logLine.append(std::string(39, ' '));
			_snprintf(buff, 1024, "%15.7e     LIMITING TIME\n", TIMEL);
			logLine.append(buff);
			logLine.append(std::string(39, ' '));
			_snprintf(buff, 1024, "%15.7e     INITIAL TIME INCREMENT\n", TIMEC);
			logLine.append(buff);

			logLine.append(std::string(46, ' '));
			_snprintf(buff, 1024, "%8d     TIME INCREMENT CHANGE CYCLE\n", NTCYC);
			logLine.append(buff);
			logLine.append(std::string(39, ' '));
			_snprintf(buff, 1024, "%15.7e     TIME INCREMENT MULTIPLIER\n", TIMEI);
			logLine.append(buff);
			logLine.append(std::string(39, ' '));
			_snprintf(buff, 1024, "%15.7e     MINIMUM TIME INCREMENT\n", TIMEL);
			logLine.append(buff);
			logLine.append(std::string(39, ' '));
			_snprintf(buff, 1024, "%15.7e     MAXIMUM TIME INCREMENT\n", TIMEC);
			logLine.append(buff);
			lstWriter->add_line(logLine);
			logLine.clear();
		}
		else if (schedule_list[i]->get_type() == "TIME LIST")
		{
			schedule_list[i]->set_sbased(false);
			std::vector<char> str = temporalControl[i];
			str.push_back('\0');
			strtok(str.data(), " ");
			strtok(NULL, " ");
			strtok(NULL, " ");
			std::string CREFT;
			double SCALT;
			int NTLIST;
			CREFT = strtok(NULL, " ");
			SCALT = std::stod(strtok(NULL, " "));
			NTLIST = std::stoi(strtok(NULL, " "));
			std::vector<double> timeList;
			timeList.reserve(NTLIST);

			for (int j = 0; j < NTLIST; j++)
			{
				timeList.push_back(std::stoi(strtok(NULL, " ")));
			}

			if (CREFT == "ELAPSED")
			{
				if (schedule_list[i]->get_name() == "TIME_STEPS" && timeList[0] != 0.0)
				{
					SimulationControl::exitOnError("INP-6-7");
				}
				schedule_list[i]->set_elapsed(true);
			} else if (CREFT == "ABSOLUTE")
			{
				schedule_list[i]->set_elapsed(false);
			} else
			{
				SimulationControl::exitOnError("INP-6-6");
			}


			for (int j = 0; j < NTLIST; j++)
			{
				timeList[j] = timeList[j] * SCALT;
			}

			for (int j = 0; j < NTLIST; j++)
			{
				schedule_list[i]->add_step_time(timeList[j], timeList[j]);
			}
			logLine.append("\n\n");
			logLine.append(std::string(16, ' '));
			logLine.append("SCHEDULE ");
			logLine.append(schedule_list[i]->get_name());
			logLine.append("   TIME LIST THAT INCLUDES THE FOLLOWING ");
			logLine.append(CREFT);
			logLine.append(" TIMES (SEC):\n");

			int ctr = 1;
			bool newline = true;
			for (int j = 0; j < NTLIST; j++)
			{
				if (newline){
					logLine.append(std::string(40, ' '));
					newline = false;
				}
				_snprintf(buff, 1024, "%15.7f  ", schedule_list[i]->get_step_time()[j].second);
				logLine.append(buff);
				ctr++;
				if (ctr == 8 || i==NTLIST-1)
				{
					logLine.append("\n");
					newline = true;
					ctr = 1;
				}

			}
			lstWriter->add_line(logLine);
			logLine.clear();

		}
		else if (schedule_list[i]->get_type() == "STEP LIST")
		{
			schedule_list[i]->set_sbased(true);
			schedule_list[i]->set_elapsed(false);
			std::vector<char> str = temporalControl[i];
			str.push_back('\0');
			strtok(str.data(), " ");
			strtok(NULL, " ");
			strtok(NULL, " ");
			int NSLIST = std::stoi(strtok(NULL, " "));
			int step;
			for (int j = 0; j < NSLIST; j++)
			{
				step = std::stoi(strtok(NULL, " "));
				schedule_list[i]->add_step_time(step, step);
			}
			logLine.append("\n\n");
			logLine.append(std::string(16, ' '));
			logLine.append("SCHEDULE ");
			logLine.append(schedule_list[i]->get_name());
			logLine.append("   STEP LIST THAT INCLUDES THE FOLLOWING TIME STEPS\n");
			int ctr = 1;
			bool newline = true;
			for (int j = 0; j < NSLIST; j++)
			{
				if (newline){
					logLine.append(std::string(40, ' '));
					newline = false;
				}
				_snprintf(buff, 1024, "%8d  ", schedule_list[i]->get_step_time()[j].first);
				logLine.append(buff);
				ctr++;
				if (ctr == 8 || i==NSLIST-1 )
				{
					logLine.append("\n");
					newline = true;
					ctr = 1;
				}
				
			}
			lstWriter->add_line(logLine);
			logLine.clear();
				


		}
		else
		{
			SimulationControl::exitOnError("INP-6-9");
		}
	}


	bool TSYES = false;
	for (int j = 0; j < schedule_list.size(); j++)
	{
		if (schedule_list[j]->get_name() == "TIME_STEPS")
		{
			TSYES = true;
			time_steps_index = j;
			break;
		}
	}

	if ((ISSTRA == 0) && !TSYES)
	{
		SimulationControl::exitOnError("INP-6-14");
	}


	if (schedule_list[time_steps_index]->get_step_list_size() <= 1)
	{
		SimulationControl::exitOnError("INP-6-10");
	}

	int NSMAX = schedule_list[time_steps_index]->get_step_list_size();

	schedule_list.push_back(new Schedule("STEPS_1&UP"));
	double TREF;
	if (schedule_list[time_steps_index]->get_elapsed())
	{
		if (IREAD == 1){
			TREF = t_ics;
		}
		else
		{
			std::cout << "IPSSIM FOR COLD START ON THIS STAGE" << std::endl;
			SimulationControl::exitOnError();
		}

	}
	else
	{
		TREF = 0;
	}

	ITMAX = NSMAX - 1;
	step_time = schedule_list[time_steps_index]->get_step_time();
	double tstart_, tfinish_, delt_;
	double time;
	tstart_ = TREF + step_time[0].second;
	tfinish_ = TREF + step_time[ITMAX].second;
	delt_ = step_time[1].second - step_time[0].second;
	DELT = delt_;
	for (int j = 0; j < NSMAX; j++)
	{
		if (j >0)
		{
			if (step_time[j].second == step_time[j - 1].second)
			{
				SimulationControl::exitOnError("INP-6-12");
			}
		}

		time = TREF + step_time[0].second;
		schedule_list[time_steps_index]->set_step_at(j, j);
		if (j > 0)
			schedule_list[schedule_list.size() - 1]->add_step_time(j, time);

	}

	schedule_list.push_back(new Schedule("STEP_0"));
	schedule_list[schedule_list.size() - 1]->add_step_time(0, schedule_list[time_steps_index]->get_step_time()[0].second);
	schedule_list.push_back(new Schedule("STEP_1"));
	schedule_list[schedule_list.size() - 1]->add_step_time(0, schedule_list[time_steps_index]->get_step_time()[1].second);

	for (int j = 0; j < NSCH - NSCHAU; j++)
	{
		if (time_steps_index == j)
			continue;

		NSMAX = schedule_list[j]->get_step_list_size();
		std::vector<std::pair<int, double>> sstep_time = schedule_list[j]->get_step_time();

		if (schedule_list[j]->get_elapsed())
			TREF = tstart_;
		else
			TREF = 0.0;

		if (schedule_list[j]->get_sbased())
		{
			for (int k = 0; k < NSMAX; k++)
			{
				if (k > 0)
				{
					if (sstep_time[k].first == sstep_time[k - 1].first)
						SimulationControl::exitOnError("INP-6-12");
				}
				int step = sstep_time[k].first;
				if (step <0 || step > ITMAX)
					continue;

				time = step_time[step].second;
				schedule_list[j]->set_time_at(k, time);
			}

		}
		else
		{
			for (int k = 0; k < NSMAX; k++)
			{
				if (k>0)
				{
					if (sstep_time[k].second == sstep_time[k - 1].second)
						SimulationControl::exitOnError("INP-6-12");
				}
				time = TREF + sstep_time[k].second;
				if ((time<tstart_) || (time > tfinish_))
					continue;
				int step = FRCSTP(time);
				schedule_list[j]->set_step_at(k, step);
			}
		}
	}

	_l846:
	// Writes

	logLine.append("\n\n");
	logLine.append(std::string(16, ' '));
	logLine.append("SCHEDULE STEPS_1&UP   IDENTICAL TO SCHEDULE 'TIME STEPS', EXCEPT THAT IT OMITS TIME STEP 0;\n");
	logLine.append(std::string(41, ' '));
	logLine.append("THIS SCHEDULE IS DEFINED AUTOMATICALLY BY IPSSIM\n");
	logLine.append("\n\n");
	logLine.append(std::string(16, ' '));
	logLine.append("SCHEDULE STEP_0   CONSISTS ONLY TIME STEP 0;\n");
	logLine.append(std::string(41, ' '));
	logLine.append("THIS SCHEDULE IS DEFINED AUTOMATICALLY BY IPSSIM\n");
	logLine.append("\n\n");
	logLine.append(std::string(16, ' '));
	logLine.append("SCHEDULE STEP_1   CONSISTS ONLY TIME STEP 1;\n");
	logLine.append(std::string(41, ' '));
	logLine.append("THIS SCHEDULE IS DEFINED AUTOMATICALLY BY IPSSIM\n");

	logLine.append("\n\n");
	logLine.append(std::string(13, ' '));
	logLine.append("SOLUTION CYCLINC DATA:\n\n");
	logLine.append(std::string(11, ' '));
	_snprintf(buff, 1024, "%15d     FLOW SOLUTION CYCLE (IN TIME STEPS)\n", NPCYC);
	logLine.append(buff);
	logLine.append(std::string(11, ' '));
	_snprintf(buff, 1024, "%15d     TRANSPORT SOLUTION CYCLE (IN TIME STEPS)\n\n", NUCYC);
	logLine.append(buff);
	logLine.append(std::string(16, ' '));
	logLine.append("FLOW AND TRANSPORT SOLUTIONS ARE ALSO COMPUTED AUTOMATICALLY ON TIME STEPS ON WHICH FLOW-RELATED\n");
	logLine.append(std::string(16, ' '));
	logLine.append("AND TRANSPORT-RELATED BOUNDARY CONDITIONS, RESPECTIVELY, ARE SET IN (OPTIONAL) BCS FILES.\n");
	lstWriter->add_line(logLine);
	logLine.clear();

	//print iteration and temporal
	logLine.append("\n\n\n\n           I T E R A T I O N   C O N T R O L   D A T A\n");
	if (ITRMAX - 1 <= 0){
		
		logLine.append("\n\n");
		logLine.append(std::string(11, ' '));
		logLine.append("  NON-ITERATIVE SOLUTION\n");
	} else
	{
		logLine.append("\n\n");
		_snprintf(buff, 1024, "%15d     MAXIMUM NUMBER OF ITERATIONS PER TIME STEP\n", ITRMAX);
		logLine.append(buff);
		_snprintf(buff, 1024, "%15.4e     ABSOLUTE CONVERGENCE CRITERION FOR FLOW SOLUTION\n", RPMAX);
		logLine.append(buff);
		_snprintf(buff, 1024, "%15.4e     ABSOLUTE CONVERGENCE CRITERION FOR TRANSPORT SOLUTION\n", RUMAX);
		logLine.append(buff);

	}
	lstWriter->add_line(logLine);
	logLine.clear();

	logLine.append("\n\n\n\n           S O L V E R - R E L A T E D   P A R A M E T E R S\n\n");
	if (std::string(p_solver_string.begin(),p_solver_string.end()) != "DIRECT"){
		logLine.append(std::string(13, ' '));
		logLine.append("SOLVER FOR P: ");
		logLine.append(std::string(p_solver_string.begin(),p_solver_string.end()));
		logLine.append("\n\n");
		logLine.append(std::string(20, ' '));
		_snprintf(buff, 1024, "%6d     MAXIMUM NUMBER OF MATRIX SOLVER ITERATIONS DURING P SOLUTION\n", max_p_iterations);
		logLine.append(buff);
		logLine.append(std::string(11, ' '));
		_snprintf(buff, 1024, "%15.4e     CONVERGENCE TOLERANCE FOR MATRIX SOLVER ITERATIONS DURING P SOLUTION\n\n", p_tolerance);
		logLine.append(buff);
		logLine.append(std::string(13, ' '));
		logLine.append("SOLVER FOR U: ");
		logLine.append(std::string(u_solver_string.begin(), u_solver_string.end()));
		logLine.append("\n\n");
		logLine.append(std::string(20, ' '));
		_snprintf(buff, 1024, "%6d     MAXIMUM NUMBER OF MATRIX SOLVER ITERATIONS DURING U SOLUTION\n", max_u_iterations);
		logLine.append(buff);
		logLine.append(std::string(11, ' '));
		_snprintf(buff, 1024, "%15.4e     CONVERGENCE TOLERANCE FOR MATRIX SOLVER ITERATIONS DURING U SOLUTION\n", u_tolerance);
		logLine.append(buff);
	}
	else
	{
		logLine.append("\n             SOLVER FOR P AND U : DIRECT");
	}
	lstWriter->add_line(logLine);
	logLine.clear();

	if (ISSFLO == 1)
	{
		NPCYC = ITMAX + 1;
		NUCYC = 1;
	}

	NSAVEP = 10;
	NSAVEU = 10;

	NBCN = NPBC + NUBC + 1;
	NSOP = NSOP + 1;
	NSOU = NSOU + 1;
	NOBSN = NOBS + 1;

	

	if (KTYPE[0] == 3)
	{
		N48 = 8;
		NEX = NE;
	} else
	{
		N48 = 4;
		NEX = 1;
	}

	NIN = NE * N48;

	if (bcs_defined)
	{
		Parser::instance()->mapFile(InputFiles::instance()->getFilesForReading()["BCS"]);
		Parser::instance()->extractBCS();
		SimulationControl::wConsole("\t BCS are extracted.", 10);
	}

	// Create Nodes
	nodeContainer.reserve(NN);
	allocate_node_arrays();
	if (KTYPE[0] == 3){
		NRTEST = 1;
		for (int j = 0; j < nodeData.size();j++)
		{
			
			node_num[j] = std::stoi(strtok(nodeData[j].data(), " "));
			node_nreg[j] = std::stoi(strtok(NULL, " "));
			node_x[j] = std::stod(strtok(NULL, " "));
			node_y[j] = std::stod(strtok(NULL, " "));
			node_z[j] = std::stod(strtok(NULL, " "));
			node_por[j] = std::stod(strtok(NULL, " "));
			if (node_nreg[j] != NROLD)
				NRTEST++;
			NROLD = node_nreg[j];
		}
	} else
	{
		//2D Nodes
	}

	if (KTYPE[0] == 3){

		for (int j = 0; j < nodeData.size(); j++)
		{
			node_x[j] = node_x[j] * SCALX;
			node_y[j] = node_y[j] * SCALY;
			node_z[j] = node_z[j] * SCALZ;
			node_por[j] = node_por[j] * PORFAC;
			node_sop[j] = (1.0 - node_por[j])*COMPMA + node_por[j] * COMPFL;
		}
	}
	else
	{
		//2D Nodes
	}



	// define incidence
	incidenceContainer.reserve(NE);
		for (int j = 1; j < incidenceData.size();j++)
		{
			std::vector<int> nodes;
			int el_numm = std::stoi(strtok(incidenceData[j].data(), " "));
			for (int k = 0; k < N48; k++)
			{
				nodes.push_back(std::stoi(strtok(NULL, " ")));
			}
			incidenceContainer.push_back(std::pair<int, std::vector<int>>(el_numm, nodes));
		}

	
		IBCUBC = new int[NBCN];
		IBCPBC = new int[NBCN];
		IBCSOP = new int[NSOP];
		IBCSOU = new int[NSOU];
	// define iqsop,iqsou,ipbc,iubc,iidpbc,iidubc,iidsop,iidsou,ibcpbc,ibcubc,ibcsop,ibcsou
	 // !!TODO
	BCSTR = std::vector<bool>(ITMAX + 1,true);
	BCSFL = std::vector<bool>(ITMAX + 1, true);

	logLine.append("\n\n\n\n           O U T P U T   C O N T R O L S   A N D   O P T I O N S\n\n             .LST FILE\n             ---------\n\n");
	_snprintf(buff, 1024, "             %8d   PRINTED OUTPUT CYCLE (IN TIME STEPS)\n",get_node_output_every());
	logLine.append(buff);

	if (simulation_output_controls[0][0] == 'Y'){
		KNODAL = +1;
		logLine.append("\n            - PRINT NODE COORDINATES, THICKNESSES AND POROSITIES");
	}
	else{
		KNODAL = 0;
		logLine.append("\n            - CANCEL PRINT OF NODE COORDINATES, THICKNESSES AND POROSITIES");
	}

	if (simulation_output_controls[1][0] == 'Y'){
		KELMNT = +1;
		logLine.append("\n            - PRINT ELEMENT PERMEABILITIES AND DISPERSIVITIES");
	}
	else{
		KELMNT = 0;
		logLine.append("\n            - CANCEL PRINT OF ELEMENT PERMEABILITIES AND DISPERSIVITIES");
	}

	if (simulation_output_controls[2][0] == 'Y'){
		KINCID = +1;
		logLine.append("\n            - PRINT NODE INCIDENCES IN EACH ELEMENT\n");
	}
	else{
		KINCID = 0;
		logLine.append("\n            - CANCEL PRINT OF NODE INCIDENCES IN EACH ELEMENT\n");
	}

	if (simulation_output_controls[3][0] == 'Y'){
		KPANDS = +1;
		logLine.append("\n            - PRINT PRESSURES AND SATURATIONS AT NODES ON EACH TIME STEP WITH OUTPUT");
	}
	else{
		KPANDS = 0;
		logLine.append("\n            - CANCEL PRINT OF PRESSURES AND SATURATIONS");
	}

	if (simulation_output_controls[4][0] == 'Y'){
		KVEL = +1;
		logLine.append("\n            - CALCULATE AND PRINT VELOCITIES AT ELEMENT CENTROIDS ON EACH TIME STEP WITH OUTPUT");

	}
	else{
		KVEL = 0;
		logLine.append("\n            - CANCEL PRINT OF VELOCITIES");
	}

	if (simulation_output_controls[5][0] == 'Y'){
		KCORT = +1;
		logLine.append("\n            - PRINT CONCENTRATIONS AT NODES AT EACH TIME STEP WITH OUTPUT\n");
	}
	else{
		KCORT = 0;
		logLine.append("\n            - CANCEL PRINT OF CONCENTRATIONS\n");
	}

	if (simulation_output_controls[6][0] == 'Y'){
		KBUDG = +1;
		logLine.append("\n            - CALCULATE AND PRINT FLUID AND SOLUTE BUDGETS ON EACH TIME STEP WITH OUTPUT\n");
	}
	else{
		KBUDG = 0;
		logLine.append("\n            - CANCEL PRINT OF BUDGETS\n");
	}

	if (simulation_output_controls[7][0] == 'Y'){
		KSCRN = +1;
	}
	else{
		KSCRN = 0;
	}

	if (simulation_output_controls[8][0] == 'Y'){
		KPAUSE = +1;
	}
	else{
		KPAUSE = 0;
	}
	lstWriter->add_line(logLine);
	logLine.clear();

	// Set NODAL OUTPUT HEADERS
	// !!TODO

	//	node_output_headers  
	//
	NCOLS5 = node_output_headers.size() - 1;
	LCOLS6 = element_output_headers.size() - 1;

	logLine.append("\n\n            .NOD FILE\n             ---------\n\n");
	_snprintf(buff, 1024, "%8d   PRINTED OUTPUT CYCLE (IN TIME STEPS)\n",node_output_every);

	for (int i = 0; i < NCOLS5; i++)
	{
		bool found = false;
		for (int j = 0; j < 9; j++)
		{
			if (std::string(node_output_headers[i].begin(),node_output_headers[i].end()) == K5SYM[j])
			{
				if (j == 0 && i != 0)
				{
					SimulationControl::exitOnError("INP-8B-1");
				}
				if (j == 4 && KTYPE[0] == 2)
				{
					SimulationControl::exitOnError("INP-8B-2");
				}
				J5COL[i] = j;
				found = true;
				break;
			}
		}
		if (!found)
		SimulationControl::exitOnError("INP-8B-3");
	}

	for (int j = 0; j < NCOLS5; j++)
	{
		_snprintf(buff, 1024, "             COLUMN %1d :  ", j+1);
		logLine.append(buff);
		logLine.append(VARNK5[J5COL[j]]);
		logLine.append("\n");
	}

	//Set ELE OUTPUT HEADERS
	// !!TODO
	logLine.append("\n\n            .ELE FILE\n             ---------\n\n");
	_snprintf(buff, 1024, "%8d   PRINTED OUTPUT CYCLE (IN TIME STEPS)\n");

	for (int i = 0; i < LCOLS6; i++)
	{
		bool found = false;
		for (int j = 0; j < 7; j++)
		{
			if (std::string(element_output_headers[i].begin(), element_output_headers[i].end()) == K6SYM[j])
			{
				if (j == 0 && i != 0)
				{
					SimulationControl::exitOnError("INP-8C-1");
				}
				if (j == 4 && KTYPE[0] == 2)
				{
					SimulationControl::exitOnError("INP-8C-2");
				}
				J6COL[i] = j;
				found = true;
				break;
			}
		}
		if (!found)
			SimulationControl::exitOnError("INP-8C-3");
	}
	for (int j = 0; j < LCOLS6; j++)
	{
		_snprintf(buff, 1024, "             COLUMN %1d :  ", j+1);
		logLine.append(buff);
		logLine.append(VARNK6[J6COL[j]]);
		logLine.append("\n");
	}


	lstWriter->add_line(logLine);
	logLine.clear();


	NOBCYC = ITMAX + 1;
	if (NOBSN - 1 == 0)
		goto _999;

	NOBS = NOBSN - 1;
	NOBCYC = -1;
	{ std::string obs_nam,obs_sch,obs_fmt;
	double obs_x, obs_y, obs_z;
		for (std::vector<std::vector<char>> str : obsData)
		{
			obs_nam = std::string(str[0].begin(),str[0].end());
			obs_x = std::stod(std::string(str[1].begin(), str[1].end()));
			obs_y = std::stod(std::string(str[2].begin(), str[2].end()));
			obs_z = std::stod(std::string(str[3].begin(), str[3].end()));
			obs_sch = std::string(str[4].begin(), str[4].end());
			obs_fmt = std::string(str[5].begin(), str[5].end());
			obsContainer.push_back(obsPoint(obs_nam, obs_x, obs_y, obs_z, obs_sch, obs_fmt));
		}
	}

	// Check if Obs Schedule defined
	//!! TODO
	if (ISSTRA == 1){
		logLine.append("\n\n          .OBS AND .OBC FILES \n             -------------------\n\n");
		logLine.append(std::string(13, ' '));
		logLine.append("UNIT ASSIGNED TO ");
		logLine.append(InputFiles::instance()->getFilesForWriting()["OBS"]);
		logLine.append("\n\n             NOTE: BECAUSE FLOW AND TRANSPORT ARE STEADY-STATE, USER-DEFINED SCHEDULES ARE NOT IN EFFECT.\n");
		logLine.append("             STEADY-STATE OBSERVATIONS WILL BE WRITTEN TO THE APPROPRIATE OUTPUT FILES.");
	}
	else
	{
		logLine.append("\n\n          .OBS AND .OBC FILES \n             -------------------\n\n");
		logLine.append(std::string(13, ' '));
		logLine.append("SCHEDULE ");
		logLine.append(obsContainer[0].get_obs_sch());
		logLine.append(" , FORMAT OBS");
		logLine.append(", ASSIGNED TO ");
		logLine.append(InputFiles::instance()->getFilesForWriting()["OBS"]);
		logLine.append("\n");

	}

	logLine.append("\n\n\n\n           O B S E R V A T I O N   P O I N T S\n");
	if (KTYPE[0] == 3){
		logLine.append("\n\n             NAME");
		logLine.append(std::string(44, ' '));
		logLine.append("COORDINATES");
		logLine.append(std::string(43, ' '));
		logLine.append("SCHEDULE");
		logLine.append(std::string(4, ' '));
		logLine.append("FORMAT\n");
		logLine.append(std::string(13, ' '));
		logLine.append("----");
		logLine.append(std::string(44, ' '));
		logLine.append("-----------");
		logLine.append(std::string(43, ' '));
		logLine.append("--------");
		logLine.append(std::string(4, ' '));
		logLine.append("------\n");
	}

	if (NOBCYC != -1)
	{
		
	}
	else
	{
		for (obsPoint obs : obsContainer)
		{
			_snprintf(buff,1024,"             %10s %s ( %+14.7e, %+14.7e, %+14.7e )   %s  %s\n",obs.get_name().c_str(),std::string(35,'.').c_str(),obs.get_x(),obs.get_y(),obs.get_z(),obs.get_obs_sch().c_str(),obs.get_format().c_str());
			logLine.append(buff);
		}
	}

	lstWriter->add_line(logLine);
	logLine.clear();

	logLine.append("             .BCOF, .BCOS, .BCOP, AND .BCOU FILES\n");
	logLine.append("             ------------------------------------\n\n");
	_snprintf(buff, 1024, "                %4d   PRINTED OUTPUT CYCLE FOR FLUID SOURCES/SINK NODES TO .BCOF FILE (IN TIME STEPS)\n", NBCFPR);
	logLine.append(buff);
	_snprintf(buff, 1024, "                %4d   PRINTED OUTPUT CYCLE FOR SOLUTE SOURCES/SINK NODES TO .BCOS FILE (IN TIME STEPS)\n", NBCSPR);
	logLine.append(buff);
	_snprintf(buff, 1024, "                %4d   PRINTED OUTPUT CYCLE FOR SPECIFIED PRESSURE NODES TO .BCOP FILE (IN TIME STEPS)\n", NBCPPR);
	logLine.append(buff);
	_snprintf(buff, 1024, "                %4d   PRINTED OUTPUT CYCLE FOR SPECIFIED CONCENTRATION NODES TO .BCOU FILE (IN TIME STEPS)\n", NBCUPR);
	logLine.append(buff);

	if (CINACT == 'Y')
		logLine.append("\n             - PRINT INACTIVE BOUNDARY CONDITIONS\n");
	else
		logLine.append("\n             - CANCEL PRINT OF INACTIVE BOUNDARY CONDITIONS\n");
	lstWriter->add_line(logLine);
	logLine.clear();
	// ADSMOD check
	// !!TODO
	_999:

	if (ME <= 0)
	{
		CS = 0.0;
		CW = 1.0;
		SIGMAS = 0.0;
	} else
	{
		adsorption_string = { 'N', 'O', 'N', 'E' };
		CHI1 = 0.0;
		CHI2 = 0.0;
		PRODF1 = 0.0;
		PRODS1 = 0.0;
	}


	logLine.append("\n\n\n\n           C O N S T A N T   P R O P E R T I E S   O F   F L U I D   A N D   S O L I D   M A T R I X\n\n");
	_snprintf(buff, 1024, "           %+15.4e     COMPRESSIBILITY OF FLUID\n",COMPFL);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     COMPRESSIBILITY OF POROUS MATRIX\n\n",COMPMA);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     FLUID VISCOSITY\n\n",VISC0);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     DENSITY OF A SOLID GRAIN\n\n",RHOS);
	logLine.append(buff);

	logLine.append("            FLUID DENSITY, RHOW\n");
	logLine.append("            CALCULATED BY IPSSIM IN TERMS OF SOLUTE CONCENTRATION,U , AS:\n");
	logLine.append("            RHOW = RHOW0 + DRWDU*(U-URHOW0)\n\n");

	_snprintf(buff, 1024, "           %+15.4e     FLUID BASE DENSITY, RHOW0\n",RHOW0);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     COEFFICIENT OF DENSITY CHANGE WITH SOLUTE CONCENTRATION, DRWDU\n",DRWDU);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     SOLUTE CONCENTRATION, URHOW0, AT WHICH FLUID DENSITY IS AT BASE VALUE, RHOW0\n\n",URHOW0);
	logLine.append(buff);

	_snprintf(buff, 1024, "           %+15.4e     MOLECULAR DIFFUSIVITY OF SOLUTE IN FLUID\n",SIGMAW);
	logLine.append(buff);
	lstWriter->add_line(logLine);
	logLine.clear();

	std::string ads = std::string(adsorption_string.begin(),adsorption_string.end());
	if (ads == "NONE")
	{
		logLine.append("\n\n\n\n           A D S O R P T I O N   P A R A M E T E R S\n\n");
		logLine.append("                NON-SORBING SOLUTE\n");
		lstWriter->add_line(logLine);
		logLine.clear();
	}

	logLine.append("\n\n\n\n            P R O D U C T I O N   A N D   D E C A Y   O F   S P E C I E S   M A S S\n\n");
	logLine.append("             PRODUCTION RATE (+)\n");
	logLine.append("             DECAY RATE (-)\n");
	_snprintf(buff, 1024, "           %+15.4e     ZERO-ORDER RATE OF SOLUTE MASS PRODUCTION/DECAY IN FLUID\n", PRODF0);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     ZERO-ORDER RATE OF ADSORBATE MASS PRODUCTION/DECAY IN IMMOBILE PHASE\n", PRODS0);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     FIRST-ORDER RATE OF SOLUTE MASS PRODUCTION/DECAY IN FLUID\n", PRODF1);
	logLine.append(buff);
	_snprintf(buff, 1024, "           %+15.4e     FIRST-ORDER RATE OF ADSORBATE MASS PRODUCTION/DECAY IN IMMOBILE PHASE\n", PRODS1);
	logLine.append(buff);
	lstWriter->add_line(logLine);
	logLine.clear();

	if (KTYPE[0] == 3)
	{
		logLine.append("\n\n\n\n           C O O R D I N A T E   O R I E N T A T I O N   T O   G R A V I T Y\n\n");
		logLine.append("             COMPONENT OF GRAVITY VECTOR\n");
		logLine.append("             IN +X DIRECTION, GRAVX\n");
		_snprintf(buff, 1024, "           %+15.4e     GRAVX = -GRAV * D(ELEVATION)/DX\n\n", GRAVX);
		logLine.append(buff);
		logLine.append("             COMPONENT OF GRAVITY VECTOR\n");
		logLine.append("             IN +Y DIRECTION, GRAVY\n");
		_snprintf(buff, 1024, "           %+15.4e     GRAVY = -GRAV * D(ELEVATION)/DY\n\n", GRAVY);
		logLine.append(buff);
		logLine.append("             COMPONENT OF GRAVITY VECTOR\n");
		logLine.append("             IN +Z DIRECTION, GRAVZ\n");
		_snprintf(buff, 1024, "           %+15.4e     GRAVZ = -GRAV * D(ELEVATION)/DZ\n\n", GRAVZ);
		logLine.append(buff);
		lstWriter->add_line(logLine);
		logLine.clear();

		// print node information

		if (KNODAL == 0)
		{
			logLine.append("\n\n\n\n          N O D E   I N F O R M A T I O N\n\n");
			logLine.append("                PRINTOUT OF NODE COORDINATES, THICKNESSES AND POROSITIES CANCELLED.\n\n");
			logLine.append("                SCALE FACTORS :\n");
			logLine.append(std::string(33, ' '));
			_snprintf(buff, 1024, "%+15.4e     X-SCALE\n", SCALX);
			logLine.append(buff);
			logLine.append(std::string(33, ' '));
			_snprintf(buff, 1024, "%+15.4e     Y-SCALE\n", SCALY);
			logLine.append(buff);
			logLine.append(std::string(33, ' '));
			_snprintf(buff, 1024, "%+15.4e     Z-SCALE\n", SCALZ);
			logLine.append(buff);
			logLine.append(std::string(33, ' '));
			_snprintf(buff, 1024, "%+15.4e     POROSITY FACTOR\n", PORFAC);
			logLine.append(buff);
			lstWriter->add_line(logLine);
			logLine.clear();
		} else if (IUNSAT == 1 && KNODAL == 0 && NRTEST != 1)
		{
			logLine.append(std::string(33, ' '));
			logLine.append("MORE THAN ONE REGION OF UNSATURATED PROPERTIES HAS BEEN SPECIFIED AMONG THE NODES\n");
			lstWriter->add_line(logLine);
			logLine.clear();
			
		} else if (IUNSAT == 1 && KNODAL == 0 && NRTEST == 1)
		{
			logLine.append(std::string(33, ' '));
			logLine.append("ONLY ONE REGION OF UNSATURATED PROPERTIES HAS BEEN SPECIFIED AMONG THE NODES\n");
			lstWriter->add_line(logLine);
			logLine.clear();
		} else if (IUNSAT != 1 && KNODAL == 1)
		{
			logLine.append("\n\n           N O D E   I N F O R M A T I O N\n\n");
			logLine.append("           NODE       X                Y                 THICKNESS      POROSITY\n\n");
			for (int i = 0; i < NN; i++)
			{
				_snprintf(buff, 1024, "         %9d   %+14.5e   %+14.5e   %+14.5e      %8.5f\n", i + 1, node_x[i], node_y[i], node_z[i], node_por[i]);
			}
			lstWriter->add_line(logLine);
			logLine.clear();
		} else if (IUNSAT == 1 && KNODAL == 1)
		{
			logLine.append("\n\n           N O D E   I N F O R M A T I O N\n\n");
			logLine.append("           NODE   REGION       X                Y                 THICKNESS      POROSITY\n\n");
			for (int i = 0; i < NN; i++)
			{
				_snprintf(buff, 1024, "         %9d   %6d   %+14.5e   %+14.5e   %+14.5e      %8.5f\n", i + 1,node_nreg[i], node_x[i], node_y[i], node_z[i], node_por[i]);
			}
			lstWriter->add_line(logLine);
			logLine.clear();
		}
	} else
	{
		
	}
	
	if (KTYPE[0] == 3)
	{
		PMAXFA = std::stod(std::string(element_props[1].begin(), element_props[1].end()));
		PMIDFA = std::stod(std::string(element_props[2].begin(), element_props[2].end()));
		PMINFA = std::stod(std::string(element_props[3].begin(), element_props[3].end()));
		ANG1FA = std::stod(std::string(element_props[4].begin(), element_props[4].end()));
		ANG1FA = std::stod(std::string(element_props[5].begin(), element_props[5].end()));
		ANG1FA = std::stod(std::string(element_props[6].begin(), element_props[6].end()));
		ALMAXF = std::stod(std::string(element_props[7].begin(), element_props[7].end()));
		ALMIDF = std::stod(std::string(element_props[8].begin(), element_props[8].end()));
		ALMINF = std::stod(std::string(element_props[9].begin(), element_props[9].end()));
		ATMXF = std::stod(std::string(element_props[10].begin(), element_props[10].end()));
		ATMDF = std::stod(std::string(element_props[11].begin(), element_props[11].end()));
		ATMNF = std::stod(std::string(element_props[12].begin(), element_props[12].end()));
	} else
	{
		
	}
	
	//Create Elements
//	elementContainer.reserve(NE);
	allocate_element_arrays();
	Timer t;
	if (KTYPE[0] == 3){
		/*int el_num;
		int lreg;
		double pmax, pmid, pmin, ang1, ang2, ang3, almax, almid, almin, atmax, atmin, atmid;*/
		
		for (int j = 0; j < elementData.size();j++)
		{
			//str.push_back('\0');
	/*		el_num = std::stoi(strtok(elementData[j].data(), " "));
			lreg = std::stoi(strtok(NULL, " "));
			pmax = std::stod(strtok(NULL, " "))*PMAXFA;
			pmid = std::stod(strtok(NULL, " "))*PMIDFA;
			pmin = std::stod(strtok(NULL, " "))*PMINFA;
			ang1 = std::stod(strtok(NULL, " "))*ANG1FA;
			ang2 = std::stod(strtok(NULL, " "))*ANG2FA;
			ang3 = std::stod(strtok(NULL, " "))*ANG3FA;
			almax = std::stod(strtok(NULL, " "))*ALMAXF;
			almid = std::stod(strtok(NULL, " "))*ALMIDF;
			almin = std::stod(strtok(NULL, " "))*ALMINF;
			atmax = std::stod(strtok(NULL, " "))*ATMXF;
			atmid = std::stod(strtok(NULL, " "))*ATMDF;
			atmin = std::stod(strtok(NULL, " "))*ATMNF;*/
			el_num[j] = std::stoi(strtok(elementData[j].data(), " "));
			el_lreg[j] = std::stoi(strtok(NULL, " "));
			el_pmax[j] = std::stod(strtok(NULL, " "));
			el_pmid[j] = std::stod(strtok(NULL, " "));
			el_pmin[j] = std::stod(strtok(NULL, " "));
			el_ang1[j] = std::stod(strtok(NULL, " "));
			el_ang2[j] = std::stod(strtok(NULL, " "));
			el_ang3[j] = std::stod(strtok(NULL, " "));
			el_almax[j] = std::stod(strtok(NULL, " "));
			el_almid[j] = std::stod(strtok(NULL, " "));
			el_almin[j] = std::stod(strtok(NULL, " "));
			el_atmax[j] = std::stod(strtok(NULL, " "));
			el_atmid[j] = std::stod(strtok(NULL, " "));
			el_atmin[j] = std::stod(strtok(NULL, " "));
		
		//	elementContainer.push_back(Element(&el_num[j],&el_lreg[j],&el_pmax[j],&el_pmid[j],&el_pmin[j],&el_ang1[j],&el_ang2[j],&el_ang3[j],&el_almax[j],&el_almid[j],&el_almin[j],&el_atmax[j],&el_atmid[j],&el_atmin[j]));
		}
	}
	else
	{
		//2D Nodes
	}

	if (KTYPE[0] == 3){
		for (int j = 0; j < NE; j++)
			el_pmax[j] = el_pmax[j] * PMAXFA;
		for (int j = 0; j < NE; j++)
			el_pmid[j] = el_pmid[j] * PMIDFA;
		for (int j = 0; j < NE; j++)
			el_pmin[j] = el_pmin[j] * PMINFA;
		
		for (int j = 0; j < NE; j++)
			el_ang1[j] = el_ang1[j] * ANG1FA;
		for (int j = 0; j < NE; j++)
			el_ang2[j] = el_ang2[j] * ANG2FA;
		for (int j = 0; j < NE; j++)
			el_ang3[j] = el_ang3[j] * ANG3FA;

		for (int j = 0; j < NE; j++)
			el_atmax[j] = el_atmax[j] * ATMXF;
		for (int j = 0; j < NE; j++)
			el_atmid[j] = el_atmid[j] * ATMDF;
		for (int j = 0; j < NE; j++)
			el_atmin[j] = el_atmin[j] * ATMNF;

		for (int j = 0; j < NE; j++)
			el_almax[j] = el_almax[j] * ALMAXF;
		for (int j = 0; j < NE; j++)
			el_almid[j] = el_almid[j] * ALMIDF;
		for (int j = 0; j < NE; j++)
			el_almin[j] = el_almin[j] * ALMINF;

		for (int j = 0; j < NE; j++)
			el_pangl1[j] = el_ang1[j] * D2R;
		for (int j = 0; j < NE; j++)
			el_pangl2[j] = el_ang2[j] * D2R;
		for (int j = 0; j < NE; j++)
			el_pangl3[j] = el_ang3[j] * D2R;

		delete[] el_ang1;
		delete[] el_ang2;
		delete[] el_ang3;

		for (int j = 0; j < NE; j++)
		{
			std::vector<double> rotMat;
			ROTMAT(el_pangl1[j], el_pangl2[j], el_pangl3[j], rotMat);
			TENSYM(el_pmax[j], el_pmid[j], el_pmin[j], rotMat, &el_permxx[j], &el_permxy[j], &el_permxz[j], &el_permyx[j], &el_permyy[j], &el_permyz[j], &el_permzx[j], &el_permzy[j], &el_permzz[j]);
		}
		delete[] el_pmax;
		delete[] el_pmid;
		delete[] el_pmin;
	}


	std::cout << t << " seconds" << std::endl;

	// Element Information

	if (KTYPE[0] == 3)
	{
		if (KELMNT)
		{
			if (IUNSAT)
			{
				logLine.append("\n\n           E L E M E N T   I N F O R M A T I O N\n\n");
				logLine.append("           ELEMENT   REGION    MAXIMUM         MIDDLE          MINIMUM                  ANGLE1         ANGLE2         ANGLE3    LONGITUDINAL   LONGITUDINAL   LONGITUDINAL     TRANSVERSE     TRANSVERSE     TRANSVERSE\n");
				logLine.append("                               PERMEABILITY    PERMEABILITY    PERMEABILITY       (IN DEGREES)    (IN DEGREES)  (IN DEGREES)    DISPERSIVITY   DISPERSIVITY   DISPERSIVITY   DISPERSIVITY   DISPERSIVITY   DISPERSIVITY\n");
				logLine.append("					      																									  IN MAX-PERM    IN-MID PERM    IN MIN-PERM     IN MAX-PERM    IN-MID PERM    IN MIN-PERM\n");
				logLine.append("                                                                               (IN DEGREES)   DIRECTION      DIRECTION         DIRECTION      DIRECTION\n");
				for (int i = 0; i < NE; i++)
				{
					_snprintf(buff, 1024, "        %9d    %5d  %+14.5e  %+14.5e  %+14.5e       %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e\n", i + 1,el_lreg[i], el_pmax[i], el_pmid[i], el_pmin[i], el_ang1[i],
						el_ang2[i], el_ang3[i], el_almax[i], el_almid[i], el_almin[i], el_atmax[i], el_atmid[i], el_atmin[i]);
					logLine.append(buff);
				}
			} else
			{
				logLine.append("\n\n           E L E M E N T   I N F O R M A T I O N\n\n");
				logLine.append("           ELEMENT       MAXIMUM         MIDDLE          MINIMUM                  ANGLE1         ANGLE2         ANGLE3    LONGITUDINAL   LONGITUDINAL   LONGITUDINAL     TRANSVERSE     TRANSVERSE     TRANSVERSE\n");
				logLine.append("                         PERMEABILITY    PERMEABILITY    PERMEABILITY       (IN DEGREES)    (IN DEGREES)  (IN DEGREES)    DISPERSIVITY   DISPERSIVITY   DISPERSIVITY   DISPERSIVITY   DISPERSIVITY   DISPERSIVITY\n");
				logLine.append("																														  IN MAX-PERM    IN-MID PERM    IN MIN-PERM     IN MAX-PERM    IN-MID PERM    IN MIN-PERM\n");
				logLine.append("                                                                                                                          DIRECTION      DIRECTION      DIRECTION         DIRECTION      DIRECTION      DIRECTION\n");
				for (int i = 0; i < NE; i++)
				{
					_snprintf(buff, 1024, "        %9d  %+14.5e  %+14.5e  %+14.5e       %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e    %+11.4e\n", i + 1, el_pmax[i], el_pmid[i], el_pmin[i], el_ang1[i],
						el_ang2[i], el_ang3[i], el_almax[i], el_almid[i], el_almin[i], el_atmax[i], el_atmid[i], el_atmin[i]);
					logLine.append(buff);
				}
				lstWriter->add_line(logLine);
				logLine.clear();
			}
		} else
		{
			if (!IUNSAT){
				logLine.append("\n\n\n\n           E L E M E N T   I N F O R M A T I O N\n\n");
				logLine.append("               PRINTOUT OF ELEMENT PERMEABILITIES AND DISPERSIVITIES CANCELLED.\n\n");
				logLine.append("               SCALE FACTORS : \n");
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      MAXIMUM PERMEABILITY FACTOR\n", PMAXFA);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      MIDDLE PERMEABILITY FACTOR\n", PMIDFA);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      MINIMUM PERMEABILITY FACTOR\n", PMINFA);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      ANGLE1 FACTOR\n", ANG1FA);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      ANGLE2 FACTOR\n", ANG2FA);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      ANGLE3 FACTOR\n", ANG3FA);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      FACTOR FOR LONGITUDINAL DISPERSIVITY IN MAX-PERM DIRECTION\n", ALMAXF);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      FACTOR FOR LONGITUDINAL DISPERSIVITY IN MID-PERM DIRECTION\n", ALMIDF);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      FACTOR FOR LONGITUDINAL DISPERSIVITY IN MIN-PERM DIRECTION\n", ALMINF);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      FACTOR FOR TRANSVERSE DISPERSIVITY IN MAX-PERM DIRECTION\n", ATMXF);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      FACTOR FOR TRANSVERSE DISPERSIVITY IN MID-PERM DIRECTION\n", ATMDF);
				logLine.append(buff);
				logLine.append(std::string(33, ' '));
				_snprintf(buff, 1024, "%+15.4e      FACTOR FOR TRANSVERSE DISPERSIVITY IN MIN-PERM DIRECTION\n", ATMNF);
				logLine.append(buff);
				lstWriter->add_line(logLine);
				logLine.clear();
			}
			else{
				if (LRTEST != 1)
				{
					logLine.append("                                 MORE THAN ONE REGION OF UNSATURATED PROPERTIES HAS BEEN SPECIFIED AMONG THE ELEMENTS.\n");
					lstWriter->add_line(logLine);
					logLine.clear();
				}
				else
				{
					logLine.append("                                ONLY ONE REGION OF UNSATURATED PROPERTIES HAS BEEN SPECIFIED AMONG THE ELEMENTS.\n");
					lstWriter->add_line(logLine);
					logLine.clear();
				}
			}
		}
	} else
	{
		
	}

	NSOPI = NSOP - 1;
	NSOUI = NSOU - 1;
	IQSOPT = 1;
	IQSOUT = 1;
	if (NSOPI != 0)
	{
		logLine.append("\n\n\n\n           F L U I D   S O U R C E   D A T A\n\n\n\n");
		logLine.append("           **** NODES AT WHICH FLUID INFLOWS OR OUTFLOWS ARE SPECIFIED ****\n\n");
		logLine.append("                             DEFAULT FLUID    DEFAULT CONCENTRATION\n");
		logLine.append("                       INFLOW(+)/OUTFLOW(-)      OF INFLOWING FLUID\n");
		logLine.append("            NODE       (FLUID MASS/SECOND)  (MASS SOLUTE/MASS WATER)\n\n");

		int IQCP, IQCPA;
		double qinc, uinc;
		for (int j = 0; j < nsopData.size() - 1; j++)
		{
			IQCP = std::stoi(strtok(nsopData[j].data(), " "));
			IQCPA = abs(IQCP);
			if (IQCPA > NN)
			{
				SimulationControl::exitOnError("INP-17-1");
			}

			if (IQCP >0)
			{
				qinc = std::stod(strtok(NULL, " "));
				if (qinc > 0){
					uinc = std::stod(strtok(NULL, " "));
					_snprintf(buff, 1024, "       %9d      %+20.13e      %+20.13e\n", IQCP,qinc,uinc);
					logLine.append(buff);
				}
				else{
					uinc = 0.0;
					_snprintf(buff, 1024, "       %9d      %+20.13e\n", IQCP, qinc);
					logLine.append(buff);
				}
			} else
			{
				_snprintf(buff, 1024, "       %9d\n", IQCP);
				logLine.append(buff);
				qinc = 0.0;
				uinc = 0.0;
			}
			IQSOP.push_back(IQCP);
			if (IQCP < 0)
				IQSOPT = -1;
			node_qin[IQCPA - 1]=qinc;
			node_uin[IQCPA - 1]=uinc;
		}
		if (IQSOPT != -1)
		{
			logLine.append("\n           SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
			logLine.append("\n           DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");
		} else
		{
			logLine.append("\n\n            TIME-DEPENDENT FLUID SOURCE/SINK OR INFLOW CONCENTRATION");
			logLine.append("\n            SET IN FUNCTION BCTIME IS INDICATED BY NEGATIVE NODE NUMBER\n");
			logLine.append("\n           SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
			logLine.append("\n           DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");
		
		}

			lstWriter->add_line(logLine);
			logLine.clear();
	}

	if (NSOUI != 0)
	{
		logLine.append("\n\n\n\n\n\n\n\n           S O L U T E   S O U R C E   D A T A\n\n\n\n");
		logLine.append("\n\n\n\n           **** NODES AT WHICH SOURCES OR SINKS OF SOLUTE MASS ARE SPECIFIED ****\n\n");
		logLine.append("                            DEFAULT SOLUTE\n");
		logLine.append("                         SOURCE(+)/SINK(-)\n");
		logLine.append("            NODE      (SOLUTE MASS/SECOND)\n\n");
		int IQCU, IQCUA;
		double quinc;
		for (int j = 0; j < nsouData.size() - 1; j++)
		{
			IQCU = std::stoi(strtok(nsouData[j].data(), " "));
			IQCUA = abs(IQCU);
			if (IQCUA > NN)
				SimulationControl::exitOnError("INP-18-1");

			if (IQCU > 0)
			{
				quinc = std::stod(strtok(NULL, " "));
				_snprintf(buff, 1024, "       %9d      %+20.13e\n", IQCU, quinc);
				logLine.append(buff);
			} else
			{
				_snprintf(buff, 1024, "       %9d\n", IQCU);
				logLine.append(buff);
				quinc = 0;
			}
			IQSOU.push_back(IQCU);
			if (IQCU < 0)
				IQSOUT = -1;

			node_quin[IQCUA - 1]=quinc;
		}

		if (IQSOPT != -1)
		{
			logLine.append("\n          SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
			logLine.append("\n          DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");
		} else
		{
			logLine.append("\n\n            TIME-DEPENDENT SOLUTE SOURCE/SINK SET IN\n");
			logLine.append("            FUNCTION BCTIME IS INDICATED BY NEGATIVE NODE NUMBER.\n");
			logLine.append("\n          SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
			logLine.append("\n          DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");
		}
		lstWriter->add_line(logLine);
		logLine.clear();
	}

	if (NBCN - 1 >0)
	{
		IPBCT = 1;
		IUBCT = 1;
		if (NPBC != 0)
		{
			logLine.append("\n\n\n\n           S P E C I F I E D   P R E S S U R E   D A T A\n\n\n\n");
			logLine.append("           **** NODES AT WHICH PRESSURES ARE SPECIFIED ****\n");
			logLine.append("                (AS WELL AS SOLUTE CONCENTRATION OF ANY\n");
			logLine.append("                FLUID INFLOW WHICH MAY OCCUR AT THE POINT\n");
			logLine.append("                OF SPECIFIED PRESSURE)\n\n");
			logLine.append("            NODE      DEFAULT PRESSURE     DEFAULT CONCENTRATION\n\n");
			int IDUM, IDUMA;
			double pbc_, ubc_;
			GNUP1 = std::vector<double>(npbcData.size(), 0);
			for (int j = 0; j < npbcData.size() - 1; j++)
			{
				IDUM = std::stoi(strtok(npbcData[j].data(), " "));
				IDUMA = abs(IDUM);
				if (IDUMA > NN)
					SimulationControl::exitOnError("INP-19-1");
				IPBC.push_back(IDUM);
				if (IDUM > 0)
				{
					pbc_ = std::stod(strtok(NULL, " "));
					ubc_ = std::stod(strtok(NULL, " "));
					node_pbc[IDUMA - 1]=pbc_;
					node_ubc[IDUMA - 1]=ubc_;
					_snprintf(buff, 1024, "       %9d      %+20.13e      %+20.13e\n", IDUM,pbc_,ubc_);
					logLine.append(buff);
				} else if (IDUM < 0)
				{
					IPBCT = -1;
					_snprintf(buff, 1024, "       %9d\n", IDUM);
					logLine.append(buff);
				}
				GNUP1[j] = GNUP;
			}
			if (IPBCT != -1)
			{
				logLine.append("\n           SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
				logLine.append("\n           DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");
			} else
			{
				logLine.append("\n\n            TIME-DEPENDENT SPECIFIED PRESSURE OR INFLOW CONCENTRATION\n");
				logLine.append("            SET IN FUNCTION BCTIME IS INDICATED BY NEGATIVE NODE NUMBER\n");
				logLine.append("\n           SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
				logLine.append("\n           DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");
			}
			lstWriter->add_line(logLine);
			logLine.clear();
			QPLITR = std::vector<double>(npbcData.size(), 0);
			IBCPBC = 0;
		}
	
		if (NUBC != 0)
		{
			logLine.append("\n\n\n\n           S P E C I F I E D   C O N C E N T R A T I O N   D A T A\n\n\n\n");
			logLine.append("           **** NODES AT WHICH SOLUTE CONCENTRATIONS ARE SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID SOURCES****\n");
			logLine.append("            NODE     DEFAULT CONCENTRATION\n\n");
			int IDUM, IDUMA;
			double ubc_;
			GNUU1 = std::vector<double>(nubcData.size(), 0);
			for (int j = 0; j < nubcData.size() - 1; j++)
			{
				IDUM = std::stoi(strtok(nubcData[j].data(), " "));
				IDUMA = abs(IDUM);
				if (IDUMA > NN)
					SimulationControl::exitOnError("INP-20-1");
				IUBC.push_back(IDUM);
				if (IDUM > 0)
				{
					ubc_ = std::stod(strtok(NULL, " "));
					node_ubc[IDUMA - 1]=ubc_;
					_snprintf(buff, 1024, "       %9d       %+20.13e\n", IDUM,ubc_);
					logLine.append(buff);
				} else if (IDUM < 0)
				{
					_snprintf(buff, 1024, "       %9d\n", IDUM);
					logLine.append(buff);
					IUBCT = -1;
				}
				GNUU1[j] = GNUU;
			}
			if (IUBCT != -1)
			{
				logLine.append("\n           SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
				logLine.append("\n           DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");

				
			} else
			{
				logLine.append("\n\n            TIME-DEPENDENT SPECIFIED CONCENTRATIONS USER-PROGRAMMED IN\n");
				logLine.append("            SET IN FUNCTION BCTIME IS INDICATED BY NEGATIVE NODE NUMBER\n");
				logLine.append("\n           SPECIFICATIONS MADE IN (OPTIONAL) BCS INPUT FILES TAKE PRECEDENCE OVER THE");
				logLine.append("\n           DEFAULT VALUES LISTED ABOVE AND ANY VALUES SET IN FUNCTION BCTIME.\n");

			}
			
			IBCUBC = 0;
		}
	}


	
	




	if (!KINCID)
	{
		logLine.append("\n\n\n\n           M E S H   C O N N E C T I O N   D A T A\n\n");
		logLine.append("                PRINTOUT OF NODAL INCIDENCES CANCELLED.\n");
		lstWriter->add_line(logLine);
		logLine.clear();
	} else
	{
		logLine.append("\n\n\n\n           M E S H   C O N N E C T I O N   D A T A\n\n\n");
		logLine.append("           **** NODAL INCIDENCES ****\n");
		lstWriter->add_line(logLine);
		logLine.clear();
	}
	if (!strncmp(incidenceData[0].data(), "INCIDENCE", incidenceData[0].size()))
	{
		std::cout << "First Line of Data Set 22 must be 'INCIDENCE'.." << std::endl;
		SimulationControl::exitOnError("INP-22-1");
	}
	incidence_vector.reserve(NE*N48);
	{
		for (std::pair<int, std::vector<int>> inc_ : incidenceContainer)
		{
			_snprintf(buff, 1024, "           ELEMENT %9d     NODES AT :      CORNERS *****", inc_.first);
			logLine.append(buff);
			for (int nod_ : inc_.second){
				incidence_vector.push_back(nod_);
				_snprintf(buff, 1024, "%9d ", nod_);
				logLine.append(buff);
			}
			logLine.append("*****\n");
		}
		lstWriter->add_line(logLine);
		logLine.clear();
	}



	if (NOBCYC != -1)
	{
		std::cout << "NOBCYC != 1, IPSSIM works with 2.2 DataSet" << std::endl;
		SimulationControl::exitOnError();
	}

	if (KTYPE[0] == 3){
		for (int j = 0; j < obsContainer.size();j++)
		{
			int inout;
			double xsi_, eta_, zet_;
			for (int i = 1; i <= NE; i++)
			{
				inout = FINDL3(i, obsContainer[j],xsi_,eta_,zet_);
				if (inout == 1)
				{
					obsContainer[j].set_element(i);
					obsContainer[j].set_xsi(xsi_);
					obsContainer[j].set_zet(zet_);
					obsContainer[j].set_eta(eta_);
					break;
				}
			}
			if (inout != 1)
			{
				SimulationControl::exitOnError("INP-8D-3");
			}

		}
	} else
	{
		
	}

	if (strncmp(p_solver_string.data(),"DIRECT",p_solver_string.size()))
	{
		NDIMJA = NN + 1;
		Timer t;
		PTRSET();
		std::cout << "PTR Set done in " << t << " seconds." << std::endl;
	}

	BANWID();
	{ for (int i = 0; i <= ITMAX; i++){
		ITBCS = i;
		BCSTEP();
		}
	}
	{
		ITRST = 0;
		if (!strncmp(p_ics_string.data(), "'UNIFORM'", p_ics_string.size()))
		{
			double PUNI = p_ics[0];
			for (int i = 0; i < NN; i++)
			{
				node_pvec[i]=PUNI;
			}
		}
		else if (!strncmp(p_ics_string.data(), "'NONUNIFORM'", p_ics_string.size()))
		{
			for (int i = 0; i < NN; i++)
			{
				node_pvec[i]=p_ics[i];
			}
		} else
		{
			SimulationControl::exitOnError("ICS-2-1");
		}

		if (!strncmp(u_ics_string.data(), "'UNIFORM'", u_ics_string.size()))
		{
			double UUNI = u_ics[0];
			for (int i = 0; i < NN; i++)
			{
				node_uvec[i]=UUNI;
			}
		}
		else if (!strncmp(u_ics_string.data(), "'NONUNIFORM'", u_ics_string.size()))
		{
			for (int i = 0; i < NN; i++)
			{
				node_uvec[i]=u_ics[i];
			}
		}
		else
		{
			SimulationControl::exitOnError("ICS-3-1");
		}
		DELTP = DELT*1e16;
		DELTU = DELT*1e16;

		if (IPBCT < 0)
		{
			for (int j = 0; j < NPBC; j++)
			{
				if (IPBC[j] < 0)
					node_pbc[IPBC[j] - 1]=node_pvec[IPBC[j] - 1];
			}
		} else
		{
			for (int i = 0; i < NN; i++)
			{
				node_pm1[i]=node_pvec[i];
				node_um1[i] = node_uvec[i];
				node_um2[i] = node_uvec[i];
				node_rcit[i] = RHOW0 + DRWDU*(node_uvec[i] - URHOW0);
			}
		}
		for (int i = 0; i < NN; i++)
		{
			node_sw[i]=1.0;
			node_swt[i]=1.0;
			node_swb[i]=1.0;
			node_cs1[i]=CS;

		}

		if (IUNSAT)
		{
			// if unsaturated call
		}

		TSEC = TSTART;
	}

	std::cout << "" << std::endl;
}

void Storage::allocate_element_arrays()
{
	if (KTYPE[0] == 3){
		el_num = new int[NE];
		el_lreg = new int[NE];
		el_pmax = new double[NE];
		el_pmid = new double[NE];
		el_pmin = new double[NE];
		el_ang1 = new double[NE];
		el_ang2 = new double[NE];
		el_ang3 = new double[NE];
		el_pangl1 = new double[NE];
		el_pangl2 = new double[NE];
		el_pangl3 = new double[NE];
		el_almax = new double[NE];
		el_almid = new double[NE];
		el_almin = new double[NE];
		el_atmax = new double[NE];
		el_atmid = new double[NE];
		el_atmin = new double[NE];
		el_permxx = new double[NE];
		el_permxy = new double[NE];
		el_permxz = new double[NE];
		el_permyx = new double[NE];
		el_permyy = new double[NE];
		el_permyz = new double[NE];
		el_permzx = new double[NE];
		el_permzy = new double[NE];
		el_permzz = new double[NE];
		el_vmag = new double[NE]{};
		el_vang1 = new double[NE]{};
		el_vang2 = new double[NE]{};
		el_gxsi = std::vector<std::vector<double>>(NE,std::vector<double>(N48,0));
		el_geta = std::vector<std::vector<double>>(NE, std::vector<double>(N48, 0));
		el_gzet = std::vector<std::vector<double>>(NE, std::vector<double>(N48, 0));
		el_det = std::vector<std::vector<double>>(NE, std::vector<double>(N48, 0));
		elementContainer.reserve(NE);
	}
}
void Storage::allocate_node_arrays()
{
	if (KTYPE[0] == 3)
	{
		node_num = new int[NN];
		node_nreg = new int[NN];
		node_x = new double[NN];
		node_y = new double[NN];
		node_z = new double[NN];
		node_por = new double[NN];
		node_sop = new double[NN];
		node_qin = new double[NN]{};
		node_uin = new double[NN]{};
		node_quin = new double[NN]{};
		node_pbc = new double[NN]{};
		node_ubc = new double[NN]{};
		node_pvec = new double[NN]{};
		node_uvec = new double[NN]{};
		node_um1 = new double[NN]{};
		node_um2 = new double[NN]{};
		node_pm1 = new double[NN]{};
		node_rcit = new double[NN]{};
		node_sw = new double[NN]{};
		node_swt = new double[NN]{};
		node_swb = new double[NN]{};
		node_cnub = new double[NN]{};
		node_cnubm1 = new double[NN]{};
		node_dswdp = new double[NN]{};
		node_cs1 = new double[NN]{};
		node_cs2 = new double[NN]{};
		node_cs3 = new double[NN]{};
		node_sl = new double[NN]{};
		node_sr = new double[NN]{};
		node_dpdtitr = new double[NN]{};
		node_uiter = new double[NN]{};
		node_piter = new double[NN]{};
		node_pvel = new double[NN]{};
		node_rcitm1 = new double[NN]{};
		node_qinitr = new double[NN]{};
		node_p_rhs = new double[NN]{};
		node_u_rhs = new double[NN]{};
		node_p_solution = new double[NN]{};
		node_u_solution = new double[NN]{};
		node_vol = new double[NN]{};
		node_rho = new double[NN]{};
		node_relkb = new double[NN]{};
		node_relk = new double[NN]{};
		node_relkt = new double[NN]{};
		node_neighbors = std::vector<std::vector<int>>(NN, std::vector<int>());
	}
}

void Storage::de_allocate_node_arrays()
{ 
	if (KTYPE[0] == 3)
	{
		delete[] node_num;
		delete[] node_x;
		delete[] node_y;
		delete[] node_z;
		delete[] node_nreg;
		delete[] node_por;
		delete[] node_sop;
		delete[] node_qin;
		delete[] node_uin;
		delete[] node_quin;
		delete[] node_pbc;
		delete[] node_ubc;
		delete[] node_pvec;
		delete[] node_uvec;
		delete[] node_um1;
		delete[] node_um2;
		delete[] node_pm1;
		delete[] node_rho;
		delete[] node_relkt;
		delete[] node_relkb;
		delete[] node_relk;
	}
}


void Storage::de_allocate_element_arrays()
{

	if (KTYPE[0] == 3){
		delete[]el_num;
		delete[]el_lreg;
		delete[]el_almax;
		delete[]el_almid;
		delete[]el_almin;
		delete[]el_atmax;
		delete[]el_atmid;
		delete[]el_atmin;
		delete[]el_permxx;
		delete[]el_permxy;
		delete[]el_permxz;
		delete[]el_permyx;
		delete[]el_permyy;
		delete[]el_permyz;
		delete[]el_permzx;
		delete[]el_permzy;
		delete[]el_permzz;
		delete[]el_pangl1;
		delete[]el_pangl2;
		delete[]el_pangl3;
		delete[]el_vmag;
		delete[]el_vang1;
		delete[]el_vang2;
	
	}
}

void Storage::ROTMAT(double& a1, double& a2, double& a3, std::vector<double>& vec)
{
	
	double s1 = sin(a1);
	double s2 = sin(a2);
	double s3 = sin(a3);
	double c1 = cos(a1);
	double c2 = cos(a2);
	double c3 = cos(a3);
	vec.assign({ c1*c2, -1 * c1*s2*s3 - s1*c3, -1 * c1*s2*c3 + s1*s3, s1*c2, -1 * s1*s2*s3 + c1*c3, -1 * s1*s2*c3 - c1*s3, s2, c2*s3 - s1*c3, c2*c3 });

}


void Storage::ROTMAT(double * a1, double * a2, double * a3, double * vec)
{

	double s1 = sin(*a1);
	double s2 = sin(*a2);
	double s3 = sin(*a3);
	double c1 = cos(*a1);
	double c2 = cos(*a2);
	double c3 = cos(*a3);
	vec[0] = c1*c2;
	vec[1] = -1 * c1*s2*s3 - s1*c3;
	vec[2] = -1 * c1*s2*c3 + s1*s3;
	vec[3] = s1*c2;
	vec[4] = -1 * s1*s2*s3 + c1*c3;
	vec[5] = -1 * s1*s2*c3 - c1*s3;
	vec[6] = s2;
	vec[7] = c2*s3 - s1*c3;
	vec[8] = c2*c3;

}

void Storage::TENSYM(double& pmax,double& pmid,double& pmin,std::vector<double>& rotMat,double* permxx,double* permxy,double *permxz,double* permyx,double*permyy,double* permyz,double* permzx,double* permzy,double* permzz)
{
	*permxx = rotMat[0] * rotMat[0] * pmax +
		rotMat[1] * rotMat[1] * pmid +
		rotMat[2] * rotMat[2] *pmin;
	*permxy = rotMat[0] * rotMat[3] * pmax +
		rotMat[1] * rotMat[4] * pmid +
		rotMat[2] * rotMat[5] * pmin;
	*permxz = rotMat[0] * rotMat[6] * pmax +
		rotMat[1] * rotMat[7] * pmid +
		rotMat[2] * rotMat[8] * pmin;
	*permyy = rotMat[3] * rotMat[3] * pmax +
		rotMat[4] * rotMat[4] * pmid +
		rotMat[5] * rotMat[5] * pmin;
	*permyz = rotMat[3] * rotMat[6] * pmax +
		rotMat[4] * rotMat[7] * pmid +
		rotMat[5] * rotMat[8] * pmin;
	*permzz = rotMat[6] * rotMat[6] * pmax +
		rotMat[7] * rotMat[7] * pmid +
		rotMat[8] * rotMat[8] * pmin;
	*permyx = *permxy;
	*permzx = *permxz;
	*permzy = *permyz;
}


void Storage::TENSYM(double * pmax, double * pmid, double * pmin, double * rotMat, double* permxx, double* permxy, double *permxz, double* permyx, double*permyy, double* permyz, double* permzx, double* permzy, double* permzz)
{
	*permxx = rotMat[0] * rotMat[0] * (*pmax) +
		rotMat[1] * rotMat[1] * (*pmid) +
		rotMat[2] * rotMat[2] * (*pmin);
	*permxy = rotMat[0] * rotMat[3] * (*pmax) +
		rotMat[1] * rotMat[4] * (*pmid) +
		rotMat[2] * rotMat[5] * (*pmin);
	*permxz = rotMat[0] * rotMat[6] * (*pmax) +
		rotMat[1] * rotMat[7] * (*pmid) +
		rotMat[2] * rotMat[8] * (*pmin);
	*permyy = rotMat[3] * rotMat[3] * (*pmax) +
		rotMat[4] * rotMat[4] * (*pmid) +
		rotMat[5] * rotMat[5] * (*pmin);
	*permyz = rotMat[3] * rotMat[6] * (*pmax) +
		rotMat[4] * rotMat[7] * (*pmid) +
		rotMat[5] * rotMat[8] * (*pmin);
	*permzz = rotMat[6] * rotMat[6] * (*pmax) +
		rotMat[7] * rotMat[7] * (*pmid) +
		rotMat[8] * rotMat[8] * (*pmin);
	*permyx = *permxy;
	*permzx = *permxz;
	*permzy = *permyz;
}

void Storage::determine_tmax()
{
	if (ISSTRA == 0)
	{
		TMAX = schedule_list[time_steps_index]->get_max_time();
	}
	else
	{
		TMAX = TSTART;
	}
	TEMAX = TMAX - TSTART;
}

void Storage::check_restart()
{
	IT = ITRST;
	ITBCS = IT;
	DIT = (double)IT;
	if (IT == 0)
	{
		DELTLC = DELT;
	}
	else
	{
		std::cout << "Restart condition is not implemented in IPSSIM yet." << std::endl;
		SimulationControl::exitOnError();
	}
}

void Storage::set_flags()
{
	ONCEP = false;
	onceNOD = false;
	onceELE = false;
	SETBCS = true;
	IBCT = IQSOPT + IQSOUT + IPBCT + IUBCT;
	
}

void Storage::set_starting_time()
{
	TSECP0 = TSEC;
	TSECU0 = TSEC;
	TMIN = TSEC / 60;
	THOUR = TMIN / 60;
	TDAY = THOUR / 24;
	TWEEK = TDAY / 7;
	TMONTH = TDAY / 30.4375;
	TYEAR = TDAY / 365.25;
}

void Storage::output_initial_starting_if_transient()
{
	if (ISSTRA != 1)
	{
		if (KTYPE[0] == 3)
		{
			
		}
		else
		{
			
		}


		if (ISSFLO == 0)
		{
			outNOD();
		}
	}
}
void Storage::set_steady_state_switches()
{
	if (ISSFLO == 1)
	{
		ML = 1;
		NOUMAT = 0;
		ISSFLO = 2;
		ITER = 0;
		DLTPM1 = DELTP;
		DLTUM1 = DELTU;
		BDELP1 = 1.0;
		BDELP = 0.0;
		BDELU = 0.0;
		switch_set = true;
		if (ISSTRA != 0)
		{
			
		} else
		{
			TELAPS = TSEC - TSTART;
		}
	} else
	{
		switch_set = false;
	}
}

void Storage::simulation()
{
	Writer * smyWriter = Writer::instance("SMY");
	Writer * nodWriter = Writer::instance("NOD");
	Writer * eleWriter = Writer::instance("ELE");
	Writer * obsWriter = Writer::instance("OBS");
	std::string smyFile,nodFile,eleFile,obsFile;
	char buff[1024];
	std::string logLine;
	smyFile.append(InputFiles::instance()->getInputDirectory());
	smyFile.append(InputFiles::instance()->getFilesForWriting()["SMY"]);
	smyWriter->set_filename(smyFile);
	nodFile.append(InputFiles::instance()->getInputDirectory());
	nodFile.append(InputFiles::instance()->getFilesForWriting()["NOD"]);
	nodWriter->set_filename(nodFile);
	eleFile.append(InputFiles::instance()->getInputDirectory());
	eleFile.append(InputFiles::instance()->getFilesForWriting()["ELE"]);
	eleWriter->set_filename(eleFile);
	obsFile.append(InputFiles::instance()->getInputDirectory());
	obsFile.append(InputFiles::instance()->getFilesForWriting()["OBS"]);
	obsWriter->set_filename(obsFile);
	logLine.append("\n          ");
	logLine.append(std::string(53, '='));
	logLine.append("\n\n");
	logLine.append(std::string(26, ' '));
	logLine.append("I   P   S   S   I   M\n\n");
	logLine.append(std::string(31, ' '));
	logLine.append("Version v1.0 \n\n");
	logLine.append("          ");
	logLine.append(std::string(53, '='));
	logLine.append("\n");
	smyWriter->add_line(logLine);
	logLine.clear();


	if (switch_set)
		goto BEGIN_ITERATION;
	BEGIN_TIMESTEP:
	IT++;
	ITREL = IT - ITRST;
	ITBCS = IT;
	DIT = (double)IT;
	ITER = 0;
	ML = 0;
	NOUMAT = 0;
	if (ONCEP && ITREL > 2)
	{
		if (((IT - 1) % NPCYC != 0) && (IT % NPCYC != 0) && (!BCSFL[IT - 1]) && (!BCSFL[IT]) && ITRMAX == 1)
			NOUMAT = 1;
	}

	if (IT!=1 && ISSFLO==2)
	{
		if ((IT%NPCYC != 0) && (!BCSFL[IT]))
			ML = 2;
		if ((IT%NUCYC != 0) && (!BCSTR[IT]))
			ML = 1;
	}

	TSECM1 = TSEC;
	TSEC = schedule_list[time_steps_index]->get_time_at_step(IT);
	TMIN = TSEC / 60.0;
	THOUR = TMIN / 60;
	TDAY = THOUR / 24;
	TWEEK = TDAY / 7;
	TMONTH = TDAY / 30.4375;
	TYEAR = TDAY / 365.25;

	DELTM1 = DELT;
	DELT = TSEC - TSECM1;
	RELCHG = abs((DELT - DELTLC) / DELTLC);

	if (RELCHG > 1e-14)
	{
		DELTLC = DELT;
		NOUMAT = 0;
	}
	// ML=0 P and U ,, ML =1 P only, ML =2 U only

	if (ISSTRA != 0)
	{
		std::cout << "TIME STEP " << IT << " OF " << ITMAX << std::endl;
	} else
	{
		TELAPS = TSEC - TSTART;
		std::cout << "TIME STEP " << IT << " OF " << ITMAX << " ; ELAPSED TIME " << TELAPS << " OF " << TEMAX << " [s]" << std::endl;
	}

	if (ML == 0)
	{
		DLTPM1 = DELTP;
		DLTUM1 = DELTU;
		DELTP = TSEC - TSECP0;
		DELTU = TSEC - TSECU0;
		TSECP0 = TSEC;
		TSECU0 = TSEC;
	} else if (ML == 1)
	{
		DLTPM1 = DELTP;
		DELTP = TSEC - TSECP0;
		TSECP0 = TSEC;
	} else if (ML == 2)
	{
		DLTUM1 = DELTU;
		DELTU = TSEC - TSECU0;
		TSECU0 = TSEC;
	} else
	{
		std::cout << "Error in ML Value ." << std::endl;
		SimulationControl::exitOnError();
	}

	BDELP = (DELTP / DLTPM1)*0.5;
	BDELU = (DELTU / DLTUM1)*0.5;
	BDELP1 = BDELP + 1.0;
	BDELU1 = BDELU + 1.0;
BEGIN_ITERATION:
	ITER = ITER + 1;
	if (ITRMAX != 1)
		std::cout << " NON _LINEARITY ITERATION " << std::endl;

	if (ML == 0)
	{
		for (int j = 0; j < NN; j++)
		{
			node_dpdtitr[j] = (node_pvec[j] - node_pm1[j]) / DELTP;
			node_piter[j] = node_pvec[j];
			node_pvel[j] = node_pvec[j];
			node_uiter[j] = node_uvec[j];
			node_rcitm1[j] = node_rcit[j];
			node_rcit[j] = RHOW0 + DRWDU *(node_uiter[j] - URHOW0);
		}

		for (int i = 0; i < NPBC; i++)
		{
			int v = abs(IPBC[i]) - 1;
			QPLITR[i] = GNUP1[i] * (node_pbc[v] - node_piter[v]);
		}

		if (ITER <= 2)
		{
			for (int k = 0; k < NN; k++)
				node_qinitr[k] = node_qin[k];
		}

		if (ITER == 1)
		{
			for (int ii = 0; ii < NN; ii++)
			{
				node_piter[ii] = BDELP1 * node_pvec[ii] - BDELP*node_pm1[ii];
				node_uiter[ii] = BDELU1 * node_uvec[ii] - BDELU*node_um1[ii];

				node_dpdtitr[ii] = (node_pvec[ii] / node_pm1[ii]) / DLTPM1;
				node_pm1[ii] = node_pvec[ii];
				node_um2[ii] = node_um1[ii];
				node_um1[ii] = node_uvec[ii];
			}
		}
	}
	else if (ML == 1)
	{
		for (int i = 0; i < NN; i++)
		{
			node_pvel[i] = node_pvec[i];
			node_piter[i] = node_pvec[i];
		}
		if (ITER == 1)
		{
			for (int j = 0; j < NN; j++)
			{
				node_piter[j] = BDELP1*node_pvec[j] - BDELP*node_pm1[j];
				node_uiter[j] = node_uvec[j];
				node_rcitm1[j] = node_rcit[j];
				node_rcit[j] = RHOW0 + DRWDU*(node_uiter[j] - URHOW0);
				node_pm1[j] = node_pvec[j];
			}
		}
	}
	else if (ML == 2)
	{
		if (ITER == 1)
		{
			for (int i = 0; i < NN; i++)
				node_uiter[i] = BDELU1*node_uvec[i] - BDELU*node_um1[i];
			
		} else
		{
			for (int i = 0; i < NN; i++)
				node_uiter[i] = node_uvec[i];
		}

		if (NOUMAT == 1)
		{
			for (int i = 0; i < NN; i++){
				node_um2[i] = node_um1[i];
				node_um1[i] = node_uvec[i];
				node_cnubm1[i] = node_cnub[i];
			}
		}

		if (ITER == 1)
		{
			for (int i = 0; i < NN; i++){
				node_dpdtitr[i] = (node_pvec[i] - node_pm1[i]) / DELTP;
				node_qinitr[i] = node_qin[i];
				node_piter[i] = node_pvec[i];
				node_pvel[i] = node_pvec[i];
				node_rcitm1[i] = node_rcit[i];
			}

			for (int i = 0; i < NPBC; i++)
			{
				int v = abs(IPBC[i]) - 1;
				QPLITR[i] = GNUP1[i] * (node_pbc[v] - node_piter[v]);
			}

			for (int i = 0; i < NN; i++){
				node_um2[i] = node_um1[i];
				node_um1[i] = node_uvec[i];
				node_cnubm1[i] = node_cnub[i];
			}
		}
	}
	else
	{
		std::cout << "Error in ML Value ." << std::endl;
		SimulationControl::exitOnError();
	}

	// ZERO OUT ARRAYS
	if (ML <= 1)
	{
		init_a_val(PMAT, NELT, 0.0);
		init_a_val(node_p_rhs, NN, 0.0);
		init_a_val(node_vol, NN, 0.0);
		if (ML != 1)
		{
			if (NOUMAT <= 1)
			{
				init_a_val(UMAT, NELT, 0.0);
			}
			init_a_val(node_p_rhs, NN, 0.0);
		}
	} else
	{
		if (NOUMAT <= 1)
		{
			init_a_val(UMAT, NELT, 0.0);
		}
		init_a_val(node_p_rhs, NN, 0.0);
	}


	if (ITER == 1 && IBCT != 4)
	{
		BCTIME();
	}

	if ((ITER == 1) && bcs_defined)
	{
		BCSTEP();
	}


	if (ML != 1 && ME == -1 && NOUMAT == 0 && !strncmp(adsorption_string.data(), "NONE", 4))
		ADSORB();

	if (NOUMAT == 0)
	{
		if (KTYPE[0] == 3)
		{
			ELEMN3();
		} else
		{
			ELEMN2();
		}
	}

		NODAL();

		BC();

		std::string f_file = "p_vec";
		f_file.append(std::to_string(IT));
		f_file.append(".bin");
		
		std::ofstream outpvecbin("C:/Users/Mishac/Desktop/pvec_uvec/" + f_file, std::ios::binary);
		for (int i = 0; i < NN; i++)
			outpvecbin.write(reinterpret_cast < const char*>(&node_p_rhs[0] + i), sizeof(double));
		outpvecbin.close();

		f_file = "u_vec";
		f_file.append(std::to_string(IT));
		f_file.append(".bin");
		std::ofstream outuvecbin("C:/Users/Mishac/Desktop/pvec_uvec/" + f_file, std::ios::binary);
		for (int i = 0; i < NN; i++)
			outuvecbin.write(reinterpret_cast < const char*>(&node_u_rhs[0] + i), sizeof(double));
		outuvecbin.close();

		f_file = "u_MAT";
		f_file.append(std::to_string(IT));
		f_file.append(".bin");
		std::ofstream outumatbin("C:/Users/Mishac/Desktop/pvec_uvec/" + f_file, std::ios::binary);
		for (int i = 0; i < NELT; i++)
			outumatbin.write(reinterpret_cast < const char*>(&UMAT[0] + i), sizeof(double));
		outumatbin.close();

		f_file = "p_MAT";
		f_file.append(std::to_string(IT));
		f_file.append(".bin");
		std::ofstream outpmatbin("C:/Users/Mishac/Desktop/pvec_uvec/" + f_file, std::ios::binary);
		for (int i = 0; i < NELT; i++)
			outpmatbin.write(reinterpret_cast < const char*>(&PMAT[0] + i), sizeof(double));
		outpmatbin.close();


		
		double pnorm = DNRM2(NN, node_p_rhs, 1);
		if (pnorm == 0)
		{
			for (int i = 0; i < NN; i++)
				node_p_solution[i] = 0.0;
			std::cout << " P solution inferred from Matrix equation. No solver called." << std::endl;
		} else
		{
			CompCol_Mat_double A;
			A.newsize(NN, NN, NELT);
			for (int i = 0; i < NELT; i++)
			{
				A.val(i) = PMAT[i];
				A.row_ind(i) = IA[i];
			}
			for (int i = 0; i < NN + 1; i++)
				A.col_ptr(i) = JA[i];

			VECTOR_double b, x(A.dim(1), 0.0);
			b.newsize(NN);
			for (int i = 0; i < NN; i++)
				b(i) = node_p_rhs[i];
			int restart = 32, result, it = 2000;
			MATRIX_double H(restart + 1, restart, 0.0);
			CompCol_ILUPreconditioner_double M(A);
			double tol = 1.e-13;

			result = GMRES(A, x, b, M, H, restart, it, tol);

			if (!result)
				ONCEP = true;
			for (int i = 0; i < NN; i++)
				node_p_solution[i] = x(i);

			std::cout << " GMRES flag " << result << std::endl;
			std::cout << "iterations performed : " << it << std::endl;
			std::cout << "tolerance achieved : " << tol << std::endl;
			std::cout << "done " << std::endl;
			// if but we will try gmres
			// solve for p;
			// convert to row compressed
			/*
			re_orient_matrix(NN + 1, NELT, PMAT, JA, IA, new_MAT, row_jumper, col_indices);
			// set vienna cl rhs vector
			viennacl::vector<double> vcl_rhs = viennacl::scalar_vector<double>(NN, 0.0);
			for (int i = 0; i < NN; i++)
				vcl_rhs[i] = node_p_rhs[i];
			// Set Matrix;
			viennacl::compressed_matrix<double> A;
			A.set(row_jumper, col_indices, new_MAT, NN, NN,NELT);

			//As initial guess we take a vector consisting of all 0.9s, except for the first entry, which we set to zero:
			
			viennacl::vector<double> init_guess = viennacl::scalar_vector<double>(A.size2(), double(0.9));
			init_guess[0] = 0;

			// Set up the monitor data, holding the system matrix, the right hand side, and the initial guess:
		
			monitor_user_data<viennacl::compressed_matrix<double>, viennacl::vector<double> > my_monitor_data(A, vcl_rhs, init_guess);


			for (int i = 0; i < NN + 1; i++)
			{
				row_jumper[i] = 0;
			}
			for (int i = 0; i < NELT; i++)
			{
				col_indices[i] = 0;
				new_MAT[i] = 0;
			}
			// Preconditioner
			viennacl::linalg::ilu0_tag my_tag;
			viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double>> ilu0_precond(A, my_tag);
			viennacl::linalg::gmres_tag my_gmres_tag(1e-13, 2000,40);
			viennacl::linalg::gmres_solver<viennacl::vector<double> > my_gmres_solver(my_gmres_tag);
			my_gmres_solver.set_monitor(my_custom_monitor<viennacl::vector<double>, double, viennacl::compressed_matrix<double> >, &my_monitor_data);
			my_gmres_solver.set_initial_guess(init_guess);
			viennacl::vector<double> vcl_results(A.size2());
			vcl_results= my_gmres_solver(A, vcl_rhs, ilu0_precond);
			for (int i = 0; i < NN; i++)
				node_p_rhs[i] = vcl_results[i];
			int ITRS = my_gmres_solver.tag().iters();
			double ERR = my_gmres_solver.tag().error();*/
		}

		if (ISSFLO != 0)
		{
			for (int i = 0; i < NN; i++)
				node_pm1[i] = node_p_solution[i];
		}

		if (ML != 1)
		{
			//if (NOUMAT)
			double unorm = DNRM2(NN, node_u_rhs, 1);
			if (unorm == 0)
			{
				for (int i = 0; i < NN; i++)
					node_u_solution[i] = 0.0;
				std::cout << " U Solution inferred from matrix equation" << std::endl;
			} else
			{
				CompCol_Mat_double A;
				A.newsize(NN, NN, NELT);
				for (int i = 0; i < NELT; i++)
				{
					A.val(i) = UMAT[i];
					A.row_ind(i) = IA[i];
				}
				for (int i = 0; i < NN + 1; i++)
					A.col_ptr(i) = JA[i];

				VECTOR_double b, x(A.dim(1), 0.0);
				b.newsize(NN);
				for (int i = 0; i < NN; i++)
					b(i) = node_u_rhs[i];
				int restart = 32,result,it=1600;
				MATRIX_double H(restart + 1, restart, 0.0);
				CompCol_ILUPreconditioner_double M(A);
				double tol = 1.e-13;

				result = GMRES(A, x, b, M, H, restart, it, tol);
				
				for (int i = 0; i < NN; i++)
					node_u_solution[i] = x(i);

				std::cout << " GMRES flag " << result << std::endl;
				std::cout << "iterations performed : " << it << std::endl;
				std::cout << "tolerance achieved : " << tol << std::endl;
				std::cout << "done " << std::endl;
				//// Solve for U;
				//re_orient_matrix(NN + 1, NELT, UMAT, JA, IA, new_MAT, row_jumper, col_indices);
				//// set vienna cl rhs vector
				//viennacl::vector<double> vcl_rhs = viennacl::scalar_vector<double>(NN, 0.0);
				//for (int i = 0; i < NN; i++)
				//	vcl_rhs[i] = node_u_rhs[i];
				//// Set Matrix;
				//viennacl::compressed_matrix<double> A;
				//A.set(row_jumper, col_indices, new_MAT, NN, NN,NELT);

				////As initial guess we take a vector consisting of all 0.9s, except for the first entry, which we set to zero:
				//
				//viennacl::vector<double> init_guess = viennacl::scalar_vector<double>(A.size2(), double(0.9));
				//init_guess[0] = 0;

				//// Set up the monitor data, holding the system matrix, the right hand side, and the initial guess:
				//
				//monitor_user_data<viennacl::compressed_matrix<double>, viennacl::vector<double> > my_monitor_data(A, vcl_rhs, init_guess);


				//for (int i = 0; i < NN + 1; i++)
				//{
				//	row_jumper[i] = 0;
				//}
				//for (int i = 0; i < NELT; i++)
				//{
				//	col_indices[i] = 0;
				//	new_MAT[i] = 0;
				//}
				//// Preconditioner
				//viennacl::linalg::ilu0_tag my_tag;
				//viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double>> ilu0(A, my_tag);
				//viennacl::linalg::gmres_tag my_gmres_tag(1e-13, 1600,40);
				//viennacl::linalg::gmres_solver<viennacl::vector<double> > my_gmres_solver(my_gmres_tag);
				//my_gmres_solver.set_monitor(my_custom_monitor<viennacl::vector<double>, double, viennacl::compressed_matrix<double> >, &my_monitor_data);
				//my_gmres_solver.set_initial_guess(init_guess);
				//viennacl::vector<double> vcl_results(A.size2());
				//vcl_results = my_gmres_solver(A, vcl_rhs, ilu0);
				//for (int i = 0; i < NN; i++)
				//	node_u_rhs[i] = vcl_results[i];
				//int ITRS = my_gmres_solver.tag().iters();
				//double ERR = my_gmres_solver.tag().error();
			}
		}

		for (int i = 0; i < NN; i++)
		{
			node_cnubm1[i] = node_cnub[i];
			double ueff = max(node_u_solution[i], 1e-10);
			node_cnub[i] = node_cnubm1[i] + (-0.5 * PRODF1*(node_rho[i] * ueff / SMWH))*DELTU;

			if (IUNSAT - 2 == 0)
			{
				if (node_p_solution[i] < 0)
				{
					std::cout << "UNSAT is not implemented " << std::endl;
					SimulationControl::exitOnError("UNSAT P");
				} else
				{
					node_sw[i] = 1.0;
					node_relk[i] = 1.0;
					BUBSAT(node_swb[i], node_relkb[i], node_p_solution[i], node_cnub[i], node_relkt[i], node_swt[i], node_sw[i],node_relk[i]);

				}

			} else
			{
				node_sw[i] = 1.0;
				node_relk[i] = 1.0;
				BUBSAT(node_swb[i], node_relkb[i], node_p_solution[i], node_cnub[i], node_relkt[i], node_swt[i],node_sw[i], node_relk[i]);
			}
		}
		for (int i = 0; i < NN; i++)
		{
			node_pvec[i] = node_p_solution[i];
			node_uvec[i] = node_u_solution[i];
		}
		outELE();
		goto BEGIN_TIMESTEP;
}


void Storage::BCTIME()
{
	int NSOUI = NSOU - 1;
	int NSOPI = NSOP - 1;

	if (IPBCT < 0)
	{
		for (int i = 0; i < NPBC; i++)
		{
			int v = IPBC[i];
			if (v<0)
			{
				IBCPBC[i] = -1;
			}
		}
	}

	if (IUBCT < 0)
	{
		for (int i = 0; i < NUBC; i++)
		{
			int v = IPBC[i+NPBC];
			int g = IUBC[v];
			if (g < 0)
			{
				IBCUBC[v] = -1;
			}
		}
	}

	if (IQSOPT < 0)
	{
		for (int i = 0; i < NSOPI; i++)
		{
			int v = IQSOP[i];
			if (v < 0)
			{
				IBCSOP[i] = -1;
			}
		}
	}

	if (IQSOUT < 0)
	{
		for (int i = 0; i < NSOUI; i++)
		{
			int v = IQSOP[i];
			if (v < 0)
			{
				IBCSOU[i] = -1;
			}
		}
	}
}

void Storage::BCSTEP()
{
	if (ITBCS != 0)
	{
		BCSFL[ITBCS] = false;
		BCSTR[ITBCS] = false;
	}

	bool found = false;
	Bcs * bcs = nullptr;
	for (int i = 0; i < bcsContainer.size(); i++)
	{
		if (ITBCS == bcsContainer[i]->getTimeStep())
		{
			found = true;
			bcs = bcsContainer[i];
			break;
		}

	}
	if (found){
		int NBCN1, NSOP1, NSOU1, NSOPI, NSOUI, NSOPI1, NSOUI1;
		bool usefl, anyfl, anytr, setfl, settr;
		NBCN1 = bcs->getNumberOfPBC() + bcs->getNumberOfUBC() + 1;
		NSOP1 = bcs->getNumberOfQINC() + 1;
		NSOU1 = bcs->getNumberOfQUINC() + 1;
		NSOPI = NSOP - 1;
		NSOUI = NSOU - 1;
		NSOPI1 = NSOP1 - 1;
		NSOUI1 = NSOU1 - 1;

		usefl = ((ISSFLO != 0) && (ITBCS == 0)) || ((ISSFLO == 0) && (ITBCS != 0));
		anyfl = (NSOPI1 + bcs->getNumberOfPBC()) > 0;
		anytr = (NSOUI1 + bcs->getNumberOfUBC()) > 0;
		BCSFL[ITBCS] = usefl && anyfl;
		BCSTR[ITBCS] = anytr;
		setfl = SETBCS && BCSFL[ITBCS];
		settr = SETBCS && BCSTR[ITBCS];

		if ((bcs->getNumberOfQINC() + bcs->getNumberOfQUINC()) > 0)
		{
			/*double * qin1 = new double[NN];
			double * uin1 = new double[NN];
			double * quin1 = new double[NN];
			double * iqsop1 = new double[NSOP1];
			double * iqsou1 = new double[NSOU1];
			
			delete[] qin1;
			delete[] uin1;
			delete[] quin1;
			delete[] iqsop1;
			delete[] iqsou1;*/
			for (int i = 0; i < NSOPI1; i++)
			{
				int IQP = -1;
				int nod_ = bcs->getNodes()[i];
				for (int j = 0; j < NSOPI; j++)
				{
					int nod__ = IQSOP[j];
					if (abs(nod_) == abs(nod__))
					{
						IQP = j;
						break;
					}
				}
				if (IQP == -1)
					SimulationControl::exitOnError("BCS-3-2");

				if (setfl)
				{
					if (nod_ > 0)
					{
						node_qin[nod_ - 1] = bcs->getQINC()[i];
						if (node_qin[nod_ - 1] > 0.0)
							node_uin[nod_ - 1] = bcs->getUINC()[i];
						IBCSOP[IQP] = 1;
					} else
					{
						node_qin[-nod_ - 1] = 0.0;
						IBCSOP[IQP] = 2;
					}
				}
			}

			for (int i = 0; i < NSOUI1; i++)
			{
				int IQU = -1;
				int nod_ = bcs->getNodes()[i];
				for (int j = 0; j < NSOUI; j++)
				{
					int nod__ = IQSOU[j];
					if (abs(nod_) == abs(nod__))
					{
						IQU = j;
						break;
					}
				}
				if (IQU == -1)
					SimulationControl::exitOnError("BCS-4-2");

				if (settr)
				{
					if (nod_ > 0)
					{
						node_quin[nod_ - 1] = bcs->getQUINC()[i];
						IBCSOU[IQU] = 1;
					}
					else
					{
						node_quin[-nod_ - 1] = 0.0;
						IBCSOU[IQU] = 2;
					}
				}
			}

		}

		if (NBCN1 - 1 != 0)
		{
		/*	double * IPBC1 = new double[NBCN1];
			double * PBC1 = new double[NBCN1];
			double * IUBC1 = new double[NBCN1];
			double * UBC1 = new double[NBCN1];

			delete[] IPBC1;
			delete[] PBC1;
			delete[] IUBC1;
			delete[] UBC1;*/
			for (int i = 0; i < bcs->getNumberOfPBC(); i++)
			{
				int IP = -1;
				int nod_ = bcs->getNodes()[i];
				for (int j = 0; j < NPBC; j++)
				{
					int nod__ = IPBC[j];
					if (abs(nod_) == abs(nod__))
					{
						IP = j;
						break;
					}
				}
				if (IP == -1)
					SimulationControl::exitOnError("BCS-5-2");

				if (setfl)
				{
					if (nod_ > 0)
					{
						node_pbc[nod_ - 1] = bcs->getPBC()[i];
						node_ubc[nod_ - 1] = bcs->getUBC()[i];
						GNUP1[IP] = GNUP;
						IBCPBC[IP] = 1;
					}
					else
					{
						GNUP1[IP] = 0.0;
						IBCPBC[IP] = 2;
					}
				}
			}

			for (int i = 0; i < bcs->getNumberOfUBC(); i++)
			{
				int IP = -1;
				int nod_ = bcs->getNodes()[i];
				for (int j = 0; j < NPBC; j++)
				{
					int nod__ = IUBC[j];
					if (abs(nod_) == abs(nod__))
					{
						IP = j;
						break;
					}
				}
				if (IP == -1)
					SimulationControl::exitOnError("BCS-5-2");

				if (settr)
				{
					if (nod_ > 0)
					{
						node_ubc[nod_ - 1] = bcs->getUBC()[i];
						GNUU1[IP] = GNUU;
						IBCUBC[IP] = 1;
					}
					else
					{
						GNUU1[IP] = 0.0;
						IBCUBC[IP] = 2;
					}
				}
			}
		}
	}





}

void Storage::ADSORB()
{
	
}

void Storage::ELEMN3()
{
	int ivcalc, jvcalc, kvcalc;
	ivcalc = jvcalc = kvcalc = 0;
	if ((ML != 2) && (ITER == 1))
		ivcalc = 1;
	if (IT == 1) ivcalc = 1;
	if ((KVEL == 1))
		jvcalc = 1;

	kvcalc = ivcalc + jvcalc;

	if (INTIM)
	{
		INTIM = false;
		for (int i = 0; i < NE; i++)
		{
			double  DET[8];
			for (int j = 0; j < N48; j++)
			{
				double xloc = GXLOC[j];
				double yloc = GYLOC[j];
				double zloc = GZLOC[j];
				double CJ[9]={0,0,0,0,0,0,0,0,0};
				BASIS3_Simple(i, xloc, yloc, zloc, DET[j], CJ);
				el_gxsi[i][j] = CJ[0] * GRAVX + CJ[1] * GRAVY + CJ[2] * GRAVZ;
				el_geta[i][j] = CJ[3] * GRAVX + CJ[4] * GRAVY + CJ[5] * GRAVZ;
				el_gzet[i][j] = CJ[6] * GRAVX + CJ[7] * GRAVY + CJ[8] * GRAVZ;
				if (DET[j] <= 0)
				{
					ISTOP = ISTOP + 1;
					std::cout << "Determinant of the jacobian at node " << (incidenceContainer[i]).second[j] << " in element "
						<< i + 1 << " is negative or zero " << DET[j] << std::endl;
				}
			}
		}
	}

	if (ISTOP != 0)
	{
		SimulationControl::exitOnError("INP-14B,22-1");
	}

	if (IUNSAT != 0)
		IUNSAT = 2;
	// Main Element Loop
	double XIX, YIY, ZIZ;
	int kgx;
	double F[8][8];
	double W[8][8];
	double CJ[9];
	double DET[8];
	double DWDXG[8][8];
	double DWDYG[8][8];
	double DWDZG[8][8];
	double DFDXG[8][8];
	double DFDYG[8][8];
	double DFDZG[8][8];
	double swbg[8];
	double relkbg[8];
	double swtg[8];
	double viscg[8];
	double rhog[8];
	double relktg[8];
	double rgxg[8], rgyg[8], rgzg[8];
	double vole[8], dflowe[8], bflowe[8][8];
	double vxg[8], vyg[8], vzg[8], vgmag[8];
	double porg[8];
	double RXXG[8], RXYG[8], RXZG[8], RYXG[8], RYYG[8], RYZG[8], RZXG[8], RZYG[8], RZZG[8];
	double EXG[8], EYG[8], EZG[8];
	double bxxg[8], bxyg[8], bxzg[8];
	double byxg[8], byyg[8], byzg[8];
	double bzxg[8], bzyg[8], bzzg[8];
	double BTRANE[8][8];
	double DTRANE[8][8];
	double rddfjx, rddfjy, rddfjz;
	double xloc, yloc, zloc;
	double axsum, aysum, azsum;
	double SWTEST;
	double dxxg, dxyg, dxzg, dyxg, dyyg, dyzg, dzxg, dzyg, dzzg;
	double eswg, rhocwg, esrcg;
	double bddfjx, bddfjy, bddfjz, eddfj;
	double rxxgd, rxygd, rxzgd, ryxgd, ryygd, ryzgd, rzxgd, rzygd, rzzgd;
	double rdrx, rdry, rdrz;

	double BXXGD,BXYGD,BXZGD,BYXGD,BYYGD, BYZGD,BZXGD,BZYGD,BZZGD;
	double EXGD, EYGD, EZGD;
	for (int l = 0; l < NE; l++)
	{
		
		XIX = YIY = ZIZ = -1.0;
		kgx = 0;

		
		for (int izl = 0; izl < 2; izl++){
			for (int iyl = 0; iyl < 2; iyl++){
				for (int ixl = 0; ixl < 2; ixl++)
				{
					xloc = XIX * GLOC;
					yloc = YIY * GLOC;
					zloc = ZIZ * GLOC;
					BASIS3(1, l, xloc, yloc, zloc, F[kgx], W[kgx], DET[kgx], CJ, DFDXG[kgx], DFDYG[kgx], DFDZG[kgx], DWDXG[kgx], DWDYG[kgx], DWDZG[kgx], swbg[kgx], relkbg[kgx], vxg[kgx], vyg[kgx], vzg[kgx], vgmag[kgx], swtg[kgx], relktg[kgx], viscg[kgx], rhog[kgx], rgxg[kgx], rgyg[kgx], rgzg[kgx], porg[kgx]);
					XIX = -XIX;
					kgx++;
				}
				YIY = -YIY;
			}
			ZIZ = -ZIZ;
		}
		
		//calculate velocity at element centroid
		if (kvcalc - 2 == 0)
		{
			
			axsum = aysum = azsum = 0.0;
			for (int i = 0; i < 8; i++)
			{
				axsum = axsum + vxg[i];
				aysum = aysum + vyg[i];
				azsum = azsum + vzg[i];
			}
			el_vmag[l] = sqrt(axsum*axsum + aysum + aysum + azsum + azsum);
			if (el_vmag[l] != 0.0)
			{
				el_vang2[l] = asin(azsum / el_vmag[l]) *57.2957795130823;
				el_vmag[l] = el_vmag[l] * 0.125;
				el_vang1[l] = atan2(aysum, axsum)*57.2957795130823;
			}
			else
			{
				el_vang1[l] = 0.0;
				el_vang2[l] = 0.0;
			}

		}


		// calculate parameters for fluid mass balance at gauss points

		if (ML == 2)
			goto u_only;

		SWTEST = 0.0;
		
		for (int i = 0; i < 8; i++)
		{
			SWTEST = SWTEST + swtg[i];
			double ROMG = rhog[i] * relktg[i] / viscg[i];
			RXXG[i] = el_permxx[l] * ROMG;
			RXYG[i] = el_permxy[l] * ROMG;
			RXZG[i] = el_permxz[l] * ROMG;
			RYXG[i] = el_permyx[l] * ROMG;
			RYYG[i] = el_permyy[l] * ROMG;
			RYZG[i] = el_permyz[l] * ROMG;
			RZXG[i] = el_permzx[l] * ROMG;
			RZYG[i] = el_permzy[l] * ROMG;
			RZZG[i] = el_permzz[l] * ROMG;
		}

		// integrate fluid mass balance in an unsaturated element using asymetric weighting functions
		if (UP <= 1e-6)
			goto symmetric;

		if (SWTEST - 3.999 >= 0)
			goto symmetric;

		for (int i = 0; i < 8; i++)
		{
			vole[i] = 0.0;
			dflowe[i] = 0.0;
			for (int j = 0; j < 8; j++)
			{
				bflowe[i][j] = 0.0;
			}
		}

		for (int kg = 0; kg < 8; kg++)
		{
			 rxxgd = RXXG[kg] * DET[kg];
			 rxygd = RXYG[kg] * DET[kg];
			 rxzgd = RXZG[kg] * DET[kg];
			 ryxgd = RYXG[kg] * DET[kg];
			 ryygd = RYYG[kg] * DET[kg];
			 ryzgd = RYZG[kg] * DET[kg];
			 rzxgd = RZXG[kg] * DET[kg];
			 rzygd = RZYG[kg] * DET[kg];
			 rzzgd = RZZG[kg] * DET[kg];
			 rdrx = rxxgd*rgxg[kg] + rxygd*rgyg[kg] + rxzgd*rgzg[kg];
			 rdry = ryxgd*rgxg[kg] + ryygd*rgyg[kg] + ryzgd*rgzg[kg];
			 rdrz = rzxgd*rgxg[kg] + rzygd*rgyg[kg] + rzzgd*rgzg[kg];
			for (int i = 0; i < 8; i++)
			{
				vole[i] = vole[i] + F[kg][i] * DET[kg];
				dflowe[i] = dflowe[i] + rdrx*DWDXG[kg][i] + rdry*DWDYG[kg][i] + rdrz*DWDZG[kg][i];
			}

			for (int j = 0; j < 8; j++)
			{
			
				rddfjx = rxxgd*DFDXG[kg][j] + rxygd*DFDYG[kg][j] + rxzgd*DFDZG[kg][j];
				rddfjy = ryxgd*DFDXG[kg][j] + ryygd*DFDYG[kg][j] + ryzgd*DFDZG[kg][j];
			    rddfjz = rzxgd*DFDXG[kg][j] + rzygd*DFDYG[kg][j] + rzzgd*DFDZG[kg][j];
				for (int p = 0; p < 8; p++){
					bflowe[j][p] = bflowe[j][p] + DWDXG[kg][p] * rddfjx + DWDYG[kg][p] * rddfjy + DWDZG[kg][p] * rddfjz;
				}
			}
		}
		goto check;

	symmetric:

		for (int i = 0; i < 8; i++)
		{
			vole[i] = 0.0;
			dflowe[i] = 0.0;
			for (int j = 0; j < 8; j++)
			{
				bflowe[i][j] = 0.0;
			}
		}
		for (int kg = 0; kg < 8; kg++)
		{
			 rxxgd = RXXG[kg] * DET[kg];
			 rxygd = RXYG[kg] * DET[kg];
			 rxzgd = RXZG[kg] * DET[kg];
			 ryxgd = RYXG[kg] * DET[kg];
			 ryygd = RYYG[kg] * DET[kg];
			 ryzgd = RYZG[kg] * DET[kg];
			 rzxgd = RZXG[kg] * DET[kg];
			 rzygd = RZYG[kg] * DET[kg];
			 rzzgd = RZZG[kg] * DET[kg];
			 rdrx = rxxgd*rgxg[kg] + rxygd*rgyg[kg] + rxzgd*rgzg[kg];
			 rdry = ryxgd*rgxg[kg] + ryygd*rgyg[kg] + ryzgd*rgzg[kg];
			 rdrz = rzxgd*rgxg[kg] + rzygd*rgyg[kg] + rzzgd*rgzg[kg];
			for (int i = 0; i < 8; i++)
			{
			
				vole[i] = vole[i] + F[kg][i] * DET[kg];
				dflowe[i] = dflowe[i] + rdrx*DFDXG[kg][i] + rdry*DFDYG[kg][i] + rdrz*DFDZG[kg][i];
			}

			for (int j = 0; j < 8; j++)
			{
			
				rddfjx = rxxgd*DFDXG[kg][j] + rxygd*DFDYG[kg][j] + rxzgd*DFDZG[kg][j];
				rddfjy = ryxgd*DFDXG[kg][j] + ryygd*DFDYG[kg][j] + ryzgd*DFDZG[kg][j];
				rddfjz = rzxgd*DFDXG[kg][j] + rzygd*DFDYG[kg][j] + rzzgd*DFDZG[kg][j];
				for (int p = 0; p < 8; p++)
				{
			
					bflowe[j][p] = bflowe[j][p] + DFDXG[kg][p] * rddfjx + DFDYG[kg][p] * rddfjy + DFDZG[kg][p] * rddfjz;
				}
					
			}
		}

	check:
		if (ML == 1)
			goto send;
	u_only:
		if (NOUMAT == 1)
			goto send;

		// calculate parameters for energy balanca or solute mass balance at gauss points
		
		
	
		//add difusion and dispersion terms to total dispersion tensor
		for (int kg = 0; kg < 8; kg++)
		{
			eswg = porg[kg] * swtg[kg];
			rhocwg = rhog[kg] * CW;
			esrcg = eswg * rhocwg;
			if (vgmag[kg] <= 0)
			{
				EXG[kg] = 0.0;
				EYG[kg] = 0.0;
				EZG[kg] = 0.0;
				dxxg = dxyg = dxzg = dyxg = dyyg = dyzg = dzxg = dzyg = dzzg = 0.0;

			}
			else
			{
				EXG[kg] = esrcg*vxg[kg];
				EYG[kg] = esrcg*vyg[kg];
				EZG[kg] = esrcg*vzg[kg];
				//Timer t;
				DISPR3(&vxg[kg], &vyg[kg], &vzg[kg], &vgmag[kg], &el_pangl1[l], &el_pangl2[l], &el_pangl3[l], &el_almax[l], &el_almid[l], &el_almin[l],
					&el_atmax[l], &el_atmid[l], &el_atmin[l], &dxxg, &dxyg, &dxzg, &dyxg, &dyyg, &dyzg, &dzxg, &dzyg, &dzzg);
				//std::cout << t << std::endl;
			}
			double ESE;
			if (ME == 1)
				ESE = eswg*SIGMAW + (1.0 - porg[kg])*SIGMAS;
			else
				ESE = esrcg*SIGMAW + (1.0 - porg[kg])*rhocwg*SIGMAS;

			bxxg[kg] = esrcg*dxxg + ESE;
			bxyg[kg] = esrcg*dxyg;
			bxzg[kg] = esrcg*dxzg;
			byxg[kg] = esrcg*dyxg;
			byyg[kg] = esrcg*dyyg + ESE;
			byzg[kg] = esrcg*dyzg;
			bzxg[kg] = esrcg*dzxg;
			bzyg[kg] = esrcg*dzyg;
			bzzg[kg] = esrcg*dzzg + ESE;

		}
	
		for (int i = 0; i < 8; i++){
			for (int j = 0; j < 8; j++)
			{
				BTRANE[i][j] = 0.0;
				DTRANE[i][j] = 0.0;
			}
		}

		for (int kg = 0; kg < 8; kg++)
		{
			 BXXGD = bxxg[kg] * DET[kg];
			 BXYGD = bxyg[kg] * DET[kg];
			 BXZGD = bxzg[kg] * DET[kg];
			 BYXGD = byxg[kg] * DET[kg];
			 BYYGD = byyg[kg] * DET[kg];
			 BYZGD = byzg[kg] * DET[kg];
			 BZXGD = bzxg[kg] * DET[kg];
			 BZYGD = bzyg[kg] * DET[kg];
			 BZZGD = bzzg[kg] * DET[kg];
			 EXGD = EXG[kg] * DET[kg];
			 EYGD = EYG[kg] * DET[kg];
			 EZGD = EZG[kg] * DET[kg];
			for (int j = 0; j < 8; j++)
			{
				
				bddfjx = BXXGD*DFDXG[kg][j] + BXYGD*DFDYG[kg][j] + BXZGD*DFDZG[kg][j];
				bddfjy = BYXGD*DFDXG[kg][j] + BYYGD*DFDYG[kg][j] + BYZGD*DFDZG[kg][j];
				bddfjz = BZXGD*DFDXG[kg][j] + BZYGD*DFDYG[kg][j] + BZZGD*DFDZG[kg][j];
				eddfj = EXGD*DFDXG[kg][j] + EYGD*DFDYG[kg][j] + EZGD*DFDZG[kg][j];
				for (int i = 0; i < 8; i++)
				{
					BTRANE[j][i] = BTRANE[j][i] + DFDXG[kg][i] * bddfjx + DFDYG[kg][i] * bddfjy + DFDZG[kg][i] * bddfjz;
					DTRANE[j][i] = DTRANE[j][i] + eddfj*W[kg][i];
				}
			}
		}

		send:
			if (KSOLVP == 0)
			{
				
			} else
			{
				GLOCOL(l,vole,bflowe,dflowe,BTRANE,DTRANE);
			}
	}
}

void Storage::ELEMN2()
{
	
}

void Storage::NODAL()
{
	int JMID, IMID;
	if (IUNSAT != 0)
		IUNSAT = 1;

	if (KSOLVP == 0)
		JMID = NBHALF;
	else
		JMID = 0;

	if (NOUMAT < 1)
	{
		for (int i = 0; i < NN; i++)
		{
			if (IUNSAT - 1 == 0)
			{
				if (node_piter[i] < 0)
				{
					std::cout << "UNSATURATED is not implemented" << std::endl;
					SimulationControl::exitOnError("UNSATURATED");
				}
				else
				{
					node_sw[i] = 1.0;
					node_dswdp[i] = 0.0;
				}
			}
		}


		for (int i = 0; i < NN; i++)
		{
			BUBSAT(node_swb[i], node_relkb[i], node_piter[i], node_cnub[i], node_relkt[i], node_swt[i], node_sw[i], node_relk[i]);
		}

		for (int i = 0; i < NN; i++)
		{
			node_rho[i] = RHOW0 + DRWDU*(node_uiter[i] - URHOW0);
		}
	}


		for (int i = 0; i < NN; i++)
		{
			if (KSOLVP == 0)
				IMID = i;
			else
				IMID = JA[i];

			double SWRHON = node_swt[i] * node_rho[i];
			if (ML != 2)
			{
				double afln = (1 - ISSFLO / 2) * (SWRHON * node_sop[i] + node_por[i] * node_rho[i] * node_swb[i] * (node_dswdp[i] + node_sw[i] * (1.0 - node_swb[i]) / (PSTAR + node_piter[i])))*node_vol[i] / (TMAX - TSTART);//DELTP;
				double cfln = node_por[i] * node_swt[i] * DRWDU*node_vol[i];
				double dudt = (1 - ISSFLO / 2)*(node_um1[i] - node_um2[i]) / DLTUM1;
				cfln = cfln * dudt - (node_sw[i] * GCONST*TEMP*node_por[i] * node_rho[i] * ((node_swb[i] * node_swb[i]) / (PSTAR + node_piter[i]))*(-0.5*PRODF1*(node_rho[i] * node_uiter[i] / SMWH)))*node_vol[i];
				PMAT[IMID] = PMAT[IMID] + afln;
				node_p_rhs[i] = node_p_rhs[i] - cfln + afln * node_pm1[i] + node_qin[i];
			}
			if (ML == 1)
				continue;

			double EPRS = (1.0 - node_por[i])*RHOS;
			double atrn = (1 - ISSTRA)*(node_por[i] * SWRHON*CW + EPRS*node_cs1[i])*node_vol[i] / DELTU;
			double gtrn = node_por[i] * SWRHON*PRODF1*node_vol[i];
			double gsv = EPRS*PRODS1*node_vol[i];
			double gsltrn = gsv * node_sl[i];
			double gsrtrn = gsv * node_sr[i];
			double etrn = (node_por[i] * SWRHON*PRODF0 + EPRS*PRODS0)*node_vol[i];
			double qur = 0.0;
			double qul = 0.0;

			if (node_qinitr[i]<1)
			{
				qul = -CW*node_qinitr[i];
				qur = -qul*node_uin[i];
			}
			if (NOUMAT != 1)
				UMAT[IMID] = UMAT[IMID] + atrn - gtrn - gsltrn - qul;
			node_u_rhs[i] = node_u_rhs[i] + atrn*node_um1[i] + etrn + gsrtrn + qur + node_quin[i];
		}

	
}

void Storage::BC()
{
	int JMID, IMID;
	if (KSOLVP == 0)
		JMID = NBHALF;
	else
		JMID = 1;

		if (NPBC != 0)
		{
			for (int ip = 0; ip < NPBC; ip++)
			{
				int i = abs(IPBC[ip]);
				i = i - 1;
				if (KSOLVP == 0)
					IMID = i;
				else
					IMID = JA[i];

				if (ML < 2)
				{
					double gpinl = -GNUP1[ip];
					double gpinr = GNUP1[ip] * node_pbc[ip];
					PMAT[IMID] = PMAT[IMID] - gpinl;
					node_p_rhs[i] = node_p_rhs[i] + gpinr;
				}
				if (ML == 1)
					continue;

				double gur = 0.0;
				double gul = 0.0;
				if (QPLITR[ip] > 0)
				{
					gul = -CW*QPLITR[ip];
					gur = -gul*node_ubc[ip];
				}
				if (NOUMAT != 1)
					UMAT[IMID] = UMAT[IMID] - gul;
				node_u_rhs[i] = node_u_rhs[i] + gur;
			}
		}

		if (ML == 1)
			return;

		if (NUBC != 0)
		{
			for (int ip = 0; ip < NUBC; ip++)
			{
				int i = abs(IUBC[ip]);
				i = i - 1;
				if (KSOLVU == 0)
					IMID = i;
				else
					IMID = JA[i];

				if (NOUMAT != 1)
				{
					double GUINL = -GNUU1[ip];
					UMAT[IMID] = UMAT[IMID] - GUINL;
				}
				double guinr = GNUU1[ip] * node_ubc[ip];
				node_u_rhs[i] = node_u_rhs[i] + guinr;
			}
		}
}

void Storage::BASIS3_Simple( int L, double XLOC, double YLOC, double ZLOC, double& DET, double CJ[])
{
	double XIIX[8] = { -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0 };
	double YIIY[8] = { -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0 };
	double ZIIZ[8] = { -1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0 };

	double XF[2] = { 1.0 - XLOC, 1.0 + XLOC };
	double YF[2] = { 1.0 - YLOC, 1.0 + YLOC };
	double ZF[2] = { 1.0 - ZLOC, 1.0 + ZLOC };

	double FX[8] = { XF[0], XF[1], XF[1], XF[0], XF[0], XF[1], XF[1], XF[0] };
	double FY[8] = { YF[0], YF[0], YF[1], YF[1], YF[0], YF[0], YF[1], YF[1] };
	double FZ[8] = { ZF[0], ZF[0], ZF[0], ZF[0], ZF[1], ZF[1], ZF[1], ZF[1] };
	double DFDXL[8]{};
	double DFDYL[8]{};
	double DFDZL[8]{};
	
	for (int i = 0; i < 8; i++)
	{
		DFDXL[i] = XIIX[i] * 0.125 * FY[i] * FZ[i];
		DFDYL[i] = YIIY[i] * 0.125 * FX[i] * FZ[i];
		DFDZL[i] = ZIIZ[i] * 0.125 * FX[i] * FY[i];
	}
	for (int i = 0; i < 9; i++)
		CJ[i] = 0.0;

	for (int il = 0; il < 8; il++)
	{
		int ii = L * 8 + il;
		int i = incidence_vector[ii];
		i = i - 1;
		CJ[0] = CJ[0] + DFDXL[il] * node_x[i];
		CJ[1] = CJ[1] + DFDXL[il] * node_y[i];
		CJ[2] = CJ[2] + DFDXL[il] * node_z[i];


		CJ[3] = CJ[3] + DFDYL[il] * node_x[i];
		CJ[4] = CJ[4] + DFDYL[il] * node_y[i];
		CJ[5] = CJ[5] + DFDYL[il] * node_z[i];

		CJ[6] = CJ[6] + DFDZL[il] * node_x[i];
		CJ[7] = CJ[7] + DFDZL[il] * node_y[i];
		CJ[8] = CJ[8] + DFDZL[il] * node_z[i];
	}

	DET = CJ[0] * (CJ[4] * CJ[8] - CJ[7] * CJ[5]) - CJ[3] * (CJ[1] * CJ[8] - CJ[7] * CJ[2]) +
		CJ[6] * (CJ[1] * CJ[5] - CJ[4] * CJ[2]);
}
void Storage::BASIS3(int ICALL, int L, double XLOC, double YLOC, double ZLOC, double F[],double W[], double& DET, double CJ[],
	double DFDXG[],double DFDYG[],double DFDZG[],double DWDXG[],double DWDYG[],double DWDZG[],double& swbg,double& relkbg,
	double &vxg,double & vyg,double&vzg,double&vgmag,double& swtg,double&relktg,double &viscg,double& rhog,double&rgxg,double&rgyg,double&rgzg,double& porg)
{
	double XIIX[8] = { -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0 };
	double YIIY[8] = { -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0 };
	double ZIIZ[8] = { -1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0 };

	double XF[2] = { 1.0 - XLOC, 1.0 + XLOC };
	double YF[2] = { 1.0 - YLOC, 1.0 + YLOC };
	double ZF[2] = { 1.0 - ZLOC, 1.0 + ZLOC };

	double FX[8] = { XF[0], XF[1], XF[1], XF[0], XF[0], XF[1], XF[1], XF[0] };
	double FY[8] = { YF[0], YF[0], YF[1], YF[1], YF[0], YF[0], YF[1], YF[1] };
	double FZ[8] = { ZF[0], ZF[0], ZF[0], ZF[0], ZF[1], ZF[1], ZF[1], ZF[1] };
	double DFDXL[8]{};
	double DFDYL[8]{};
	double DFDZL[8]{};
	for (int i = 0; i < 8; i++)
		F[i] = 0.125 *FX[i] * FY[i] * FZ[i];

	for (int i = 0; i < 8; i++)
	{
		DFDXL[i] = XIIX[i] * 0.125 * FY[i] * FZ[i];
		DFDYL[i] = YIIY[i] * 0.125 * FX[i] * FZ[i];
		DFDZL[i] = ZIIZ[i] * 0.125 * FX[i] * FY[i];
	}

	//std::vector<double> CJ(9, 0);
	for (int i = 0; i < 9; i++)
		CJ[i] = 0.0;

	for (int il = 0; il < 8; il++)
	{
		int ii = L * 8 + il;
		int i = incidence_vector[ii];
		i = i - 1;
		CJ[0] = CJ[0] + DFDXL[il] * node_x[i];
		CJ[1] = CJ[1] + DFDXL[il] * node_y[i];
		CJ[2] = CJ[2] + DFDXL[il] * node_z[i];


		CJ[3] = CJ[3] + DFDYL[il] * node_x[i];
		CJ[4] = CJ[4] + DFDYL[il] * node_y[i];
		CJ[5] = CJ[5] + DFDYL[il] * node_z[i];

		CJ[6] = CJ[6] + DFDZL[il] * node_x[i];
		CJ[7] = CJ[7] + DFDZL[il] * node_y[i];
		CJ[8] = CJ[8] + DFDZL[il] * node_z[i];
	}

	DET = CJ[0] * (CJ[4] * CJ[8] - CJ[7] * CJ[5]) - CJ[3] * (CJ[1] * CJ[8] - CJ[7] * CJ[2]) +
		CJ[6] * (CJ[1] * CJ[5] - CJ[4] * CJ[2]);
	//cj = CJ;
	if (ICALL == 0){
		return;
	}

	double ODET = 1.0 / DET;
	double CIJ[9]={0,0,0,0,0,0,0,0,0};

	CIJ[0] = +ODET * (CJ[4] * CJ[8] - CJ[7] * CJ[5]);
	CIJ[1] = -ODET * (CJ[1] * CJ[8] - CJ[7] * CJ[2]);
	CIJ[2] = +ODET * (CJ[1] * CJ[5] - CJ[4] * CJ[2]);

	CIJ[3] = -ODET * (CJ[3] * CJ[8] - CJ[6] * CJ[5]);
	CIJ[4] = +ODET * (CJ[0] * CJ[8] - CJ[6] * CJ[2]);
	CIJ[5] = -ODET * (CJ[0] * CJ[5] - CJ[3] * CJ[2]);

	CIJ[6] = +ODET * (CJ[3] * CJ[7] - CJ[6] * CJ[4]);
	CIJ[7] = -ODET * (CJ[0] * CJ[7] - CJ[6] * CJ[1]);
	CIJ[8] = +ODET * (CJ[0] * CJ[4] - CJ[3] * CJ[1]);

	// CALCULATE DERIVATIVES WRT TO GLOBAL COORDINATES
	for (int i = 0; i < 8; i++)
	{
		DFDXG[i] = CIJ[0] * DFDXL[i] + CIJ[1] * DFDYL[i] + CIJ[2] * DFDZL[i];
		DFDYG[i] = CIJ[3] * DFDXL[i] + CIJ[4] * DFDYL[i] + CIJ[5] * DFDZL[i];
		DFDZG[i] = CIJ[6] * DFDXL[i] + CIJ[7] * DFDYL[i] + CIJ[8] * DFDZL[i];
	}

	// Calculate consistent components of (RHO*GRAV) term in local coordinates
	double RGXL, RGYL, RGZL, RGXLM1, RGYLM1, RGZLM1;
	RGXL = RGYL = RGZL = RGXLM1= RGYLM1 = RGZLM1 = 0;
	for (int il = 0; il < 8; il++)
	{
		int ii = L * 8 + il;
		int i = incidence_vector[ii];
		i = i - 1;
		double ADFDXL = abs(DFDXL[il]);
		double ADFDYL = abs(DFDYL[il]);
		double ADFDZL = abs(DFDZL[il]);
		RGXL = RGXL + node_rcit[i] * el_gxsi[L][il] * ADFDXL;
		RGYL = RGYL + node_rcit[i] * el_geta[L][il] * ADFDYL;
		RGZL = RGZL + node_rcit[i] * el_gzet[L][il] * ADFDZL;
		RGXLM1 = RGXLM1 + node_rcitm1[i] * el_gxsi[L][il] * ADFDXL;
		RGYLM1 = RGYLM1 + node_rcitm1[i] * el_geta[L][il] * ADFDYL;
		RGZLM1 = RGZLM1 + node_rcitm1[i] * el_gzet[L][il] * ADFDZL;
	}

	// transform consistent components to global coordinates
	rgxg = CIJ[0] * RGXL + CIJ[1] * RGYL + CIJ[2] * RGZL;
	rgyg = CIJ[3] * RGXL + CIJ[4] * RGYL + CIJ[5] * RGZL;
	rgzg = CIJ[6] * RGXL + CIJ[7] * RGYL + CIJ[8] * RGZL;

	double RGXGM1 = CIJ[0] * RGXLM1 + CIJ[1] * RGYLM1 + CIJ[2] * RGZLM1;
	double RGYGM1 = CIJ[3] * RGXLM1 + CIJ[4] * RGYLM1 + CIJ[5] * RGZLM1;
	double RGZGM1 = CIJ[6] * RGXLM1 + CIJ[7] * RGYLM1 + CIJ[8] * RGZLM1;

	// calculate Parameter values at this location
	double piterg, uiterg, dpdxg, dpdyg, dpdzg, cnubg;
	piterg = uiterg = dpdxg = dpdyg = dpdzg = porg = cnubg = 0.0;
	for (int il = 0; il < 8; il++)
	{
		int ii = L * 8 + il;
		int i = incidence_vector[ii];
		i = i - 1;
		dpdxg = dpdxg + node_pvel[i] * DFDXG[il];
		dpdyg = dpdyg + node_pvel[i] * DFDYG[il];
		dpdzg = dpdzg + node_pvel[i] * DFDZG[il];
		porg = porg + node_por[i] * F[il];
		piterg = piterg + node_piter[i] * F[il];
		uiterg = uiterg + node_uiter[i] * F[il];
		cnubg = cnubg + node_cnub[i] * F[il];
	}

	// Set values for density and viscosity
	rhog = RHOW0 + DRWDU*(uiterg - URHOW0);
	viscg = 0;
	if (ME == 1)
	{
		viscg = VISC0*239.4e-7*(pow(10.0, (248.37 / (uiterg + 133.5))));
	}
	else
	{
		viscg = VISC0;
	}

	// Set unsaturated flow parameters swg and relkg
	double relkg, swg;
	if (IUNSAT - 2 == 0)
	{
		if (piterg < 0)
		{
			std::cout << " unsaturated flow not implemented yet .." << std::endl;
			SimulationControl::exitOnError();
		}
		else
		{
			swg = 1.0;
			relkg = 1.0;
		}
	}
	else
	{
		swg = 1.0;
		relkg = 1.0;
	}
	BUBSAT(swbg,relkbg,piterg,cnubg,relktg,swtg,swg,relkg);

	// Calculate consistent fluid velocities wrt global coordinates
	double denom = 1.0 / (porg*swtg*viscg);
	double pgx = dpdxg - RGXGM1;
	double pgy = dpdyg - RGYGM1;
	double pgz = dpdzg - RGZGM1;

	if (dpdxg != 0)
		if (abs(pgx / dpdxg) - 1e-10 <= 0)
			pgx = 0.0;

	if (dpdyg != 0)
		if (abs(pgy / dpdyg) - 1e-10 <= 0)
			pgy = 0.0;

	if (dpdzg != 0)
		if (abs(pgz / dpdzg) - 1e-10 <= 0)
			pgz = 0.0;


	vxg = -denom*relktg*(el_permxx[L] * pgx + el_permxy[L] * pgy + el_permxz[L] * pgz);
	vyg = -denom*relktg*(el_permyx[L] * pgx + el_permyy[L] * pgy + el_permyz[L] * pgz);
	vzg = -denom*relktg*(el_permzx[L] * pgx + el_permzy[L] * pgy + el_permzz[L] * pgz);
	vgmag = sqrt(vxg*vxg + vyg*vyg + vzg*vzg);

	// calculate asymmetric weighting functions

	if (UP > 1.0e-6 && NOUMAT == 0)
		goto fv;

		for (int i = 0; i < 8; i++)
		{
			W[i] = F[i];
			DWDXG[i] = DFDXG[i];
			DWDYG[i] = DFDYG[i];
			DWDZG[i] = DFDZG[i];
		}

		return;
	
	fv:
	//calculate local fluid velocities
	double vxl = CIJ[0] * vxg + CIJ[1] * vyg + CIJ[2] * vzg;
	double vyl = CIJ[3] * vxg + CIJ[4] * vyg + CIJ[5] * vzg;
	double vzl = CIJ[6] * vxg + CIJ[7] * vyg + CIJ[8] * vzg;
	double vlmag = sqrt(vxl*vxl + vyl*vyl + vzl*vzl);
	double aa, bb, gg,xixi,yiyi,zizi;
	aa = bb = gg = 0.0;
	if (vlmag > 0)
	{
		aa = UP * vxl / vlmag;
		bb = UP * vyl / vlmag;
		gg = UP * vzl / vlmag;
	}
	xixi = 0.75*aa*XF[0] * XF[1];
	yiyi = 0.75*bb*YF[0] * YF[1];
	zizi = 0.75*gg*ZF[0] * ZF[1];
	double  AFX[8];
	double  AFY[8];
	double  AFZ[8];
	for (int i = 0; i < 8; i++)
	{
		AFX[i] = 0.5*FX[i] + XIIX[i] * xixi;
		AFY[i] = 0.5*FY[i] + YIIY[i] * yiyi;
		AFZ[i] = 0.5*FZ[i] + ZIIZ[i] * zizi;
	}

	// Calculate asymmetric weighting funcs
	for (int i = 0; i < 8; i++)
		W[i] = AFX[i] * AFY[i] * AFZ[i];


	double thaax = 0.5 - 1.5*aa*XLOC;
	double thbby = 0.5 - 1.5*bb*YLOC;
	double thggz = 0.5 - 1.5*gg*ZLOC;
	double XDW[8];
	double YDW[8];
	double ZDW[8];

	for (int i = 0; i < 8; i++)
	{
		XDW[i] = XIIX[i] * thaax;
		YDW[i] = YIIY[i] * thbby;
		ZDW[i] = ZIIZ[i] * thggz;
	}

	// calculate derivatives wrt local
	double DWDXL[8];
	double DWDYL[8];
	double DWDZL[8];

	for (int i = 0; i < 8; i++)
	{
		DWDXL[i] = XDW[i] * AFY[i] * AFZ[i];
		DWDYL[i] = YDW[i] * AFX[i] * AFZ[i];
		DWDZL[i] = ZDW[i] * AFX[i] * AFY[i];
	}


	// calculate derviatives wrt global;
	for (int i = 0; i < 8; i++)
	{
		DWDXG[i] = CIJ[0] * DWDXL[i] + CIJ[1] * DWDYL[i] + CIJ[2] * DWDZL[i];
		DWDYG[i] = CIJ[3] * DWDXL[i] + CIJ[4] * DWDYL[i] + CIJ[5] * DWDZL[i];
		DWDZG[i] = CIJ[6] * DWDXL[i] + CIJ[7] * DWDYL[i] + CIJ[8] * DWDZL[i];
	}

}



void Storage::BUBSAT(double& SWB,double& RELKB,double PRES,double CNUB,double & RELKT,double& SWT,double SW,double RELK)
{
	SWB = (PSTAR + PRES) / ((PSTAR + PRES) + CNUB*GCONST*TEMP);
	RELKB = pow(SWB, 2);
	SWT = SWB*SW;
	RELKT = RELKB*RELK;
}

void Storage::DISPR3(double * vx,double * vy,double * vz,double * vmag,double * ang1,double * ang2,double * ang3,
	double * ALMAX,double * ALMID,double * ALMIN,double * ATMAX,double * ATMID,double * ATMIN,
	double* dxx,double * dxy,double* dxz, double * dyx, double * dyy, double * dyz,double* dzx, double * dzy, double * dzz)
{
	//VX,VY,VZ,VMAG,ANG1,ANG2,ANG3,ALMAX,ALMID,ALMIN,  DISPR3.........800
	//  ATMAX, ATMID, ATMIN, DXX, DXY, DXZ, DYX, DYY, DYZ, DZX, DZY, DZZ
	double toliso = 1e-7;
	double tolvrt = 1e-7;
	double tolcir = 9.999999e-1;
	double dt1, dt2;
	double unx, uny, unz, wnz, wny, wnx;
	double untxx, untyy, untzz, wntxx, wntyy, wntzz;
	double vnx = *vx / *vmag;
	double vny = *vy / *vmag;
	double vnz = *vz / *vmag;
	double vnvec[3] = { 0, 0, 0 };
	double unvec[3] = { 0, 0, 0 };
	double wnvec[3] = { 0, 0, 0 };
	double vec[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	//std::vector<double> vnvec(3,0);
	//std::vector<double> unvec(3, 0);
	//std::vector<double> wnvec(3, 0);
	//std::vector<double> vec(9,0);
	bool liso = false;
	bool tiso = false;
	double AL[3] = { *ALMAX, *ALMID, *ALMIN };
	//std::vector<double> AL({ALMAX,ALMID,ALMIN});
	//double almxvl = *std::max_element(AL.begin(), AL.end());
	double almxvl = AL[0];
	double tmpal, tmpat;
	for (int i = 1; i < 3; i++)
	{
		tmpal = AL[i];
		if (tmpal > almxvl)
			almxvl = tmpal;
	}
	//double almnvl = *std::min_element(AL.begin(), AL.end());
	double almnvl = AL[0];
	for (int i = 1; i < 3; i++)
	{
		tmpal = AL[i];
		if (tmpal < almnvl)
			almnvl = tmpal;
	}
	//std::vector<double> AT({ ATMAX, ATMID, ATMIN });
	double AT[3] = { *ATMAX, *ATMID, *ATMIN };
	//double atmxvl = *std::max_element(AT.begin(), AT.end());
	double atmxvl = AT[0];
	for (int i = 1; i < 3; i++)
	{
		tmpat = AT[i];
		if (tmpat > atmxvl)
			atmxvl = tmpat;
	}
	//double atmnvl = *std::min_element(AT.begin(), AT.end());
	double atmnvl = AT[0];
	for (int i = 1; i < 3; i++)
	{
		tmpat = AT[i];
		if (tmpat < atmnvl)
			atmnvl = tmpat;
	}
	double at1, at2;
	double DL;
	int J[3];
	double VN[3];
	double UN[3], WN[3];
	//std::vector<double> rotmat;
	if (almxvl == 0)
		liso = true;
	else
		liso = ((almxvl - almnvl) / almxvl < toliso);

	
	if (liso)
	{
		DL = *ALMAX* *vmag;
	} else
	{
		
		ROTMAT(ang1, ang2, ang3, vec);
		ROTATE(&vec[0], &vnx, &vny, &vnz, &vnvec[0]);
		vnx = vnvec[0];
		vny = vnvec[1];
		vnz = vnvec[2];

		DL = *vmag / (vnvec[0] * vnvec[0] / *ALMAX + vnvec[1] * vnvec[1] / *ALMID + vnvec[2] * vnvec[2] / *ALMIN);
	}


	if (atmxvl == 0.0)
		tiso = true;
	else
		tiso = ((atmxvl - atmnvl) / atmxvl < toliso);

	if (tiso)
	{
		double term = 1.0 - vnz*vnz;
		if (term < tolvrt)
		{
			unx = wny = 1.0;
			uny = unz = wnx = wnz = 0.0;
		} else
		{
			double termh = sqrt(term);
			unx = -vny / termh;
			uny = vnx / termh;
			unz = 0.0;
			wnx = -vnz*uny;
			wny = vnz*unx;
			wnz = termh;
		}
		at1 = *ATMAX;
		at2 = at1;
	} else
	{
		if (liso)
		{
			ROTMAT(ang1, ang2, ang3, &vec[0]);
			ROTATE(&vec[0], &vnx, &vny, &vnz, &vnvec[0]);
			vnx = vnvec[0];
			vny = vnvec[1];
			vnz = vnvec[2];
		}

	
		J[0] = std::distance(AT, std::max_element(AT, AT+ sizeof(AT)/sizeof(double)));
		J[2] = std::distance(AT, std::min_element(AT, AT + sizeof(AT) / sizeof(double)));
		J[1] = 3 - J[0] - J[2];
		
		VN[0] = vnvec[0];
		VN[1] = vnvec[1];
		VN[2] = vnvec[2];
		double vntxx = VN[J[0]];
		double vntyy = VN[J[1]];
		double vntzz = VN[J[2]];
		double a2 = AT[J[0]];
		double b2 = AT[J[1]];
		double c2 = AT[J[2]];

		double a2b2 = a2*b2;
		double a2c2 = a2*c2;
		double b2c2 = b2*c2;

		double cos2av = (a2c2 - b2c2) / (a2b2 - b2c2);
		double sin2av = 1.0 - cos2av;
		double cosav = sqrt(cos2av);
		double sinav = sqrt(sin2av);
		double term1 = cosav*vntxx;
		double term2 = sinav*vntzz;
		double oa1v = term1 + term2;
		double oa2v = term1 - term2;

		if (max(abs(oa1v), abs(oa2v)) > tolcir)
		{
			untxx = -vntzz;
			untyy = 0.0;
			untzz = vntxx;
			wntxx = 0.0;
			wntyy = 1.0;
			wntzz = 0.0;
			at1 = b2;
			at2 = b2;
		} else
		{
			double rvj1mg = 1.0 / sqrt(1.0 - oa1v*oa1v);
			double rvj2mg = 1.0 / sqrt(1.0 - oa2v*oa2v);
			double rsum = rvj1mg + rvj2mg;
			double rdif = rvj1mg - rvj2mg;
			double oauxx = cosav*rsum;
			double oauzz = sinav * rdif;
			double oawxx = cosav * rdif;
			double oawzz = sinav * rsum;
			double oauv = oauxx*vntxx + oauzz * vntzz;
			double oawv = oawxx * vntxx + oawzz * vntzz;
			double oauoau = oauxx*oauxx + oauzz*oauzz;
			double oawoaw = oawxx * oawxx + oawzz * oawzz;
			double umterm = oauoau - oauv*oauv;
			double wmterm = oawoaw - oawv * oawv;

			if (umterm > wmterm)
			{
				double rumagh = 1.0 / sqrt(umterm);
				untxx = (oauxx - oauv * vntxx) * rumagh;
				untyy = -oauv*vntyy*rumagh;
				untzz = (oauzz - oauv*vntzz)*rumagh;
				wntxx = untyy*vntzz - untzz*vntyy;
				wntyy = untzz*vntxx - untxx*vntzz;
				wntzz = untxx * vntyy - untyy*vntxx;
			} else
			{
				double rwmagh = 1.0 / sqrt(wmterm);
				wntxx = (oawxx - oawv*vntxx) * rwmagh;
				wntyy = -oawv*vntyy*rwmagh;
				wntzz = (oawzz - oawv*vntzz) * rwmagh;
				untxx = wntyy*vntzz - wntzz * vntyy;
				untyy = wntzz*vntxx - wntxx * vntzz;
				untzz = wntxx*vntyy - wntyy*vntxx;

			}
			double a2b2c2 = a2b2*c2;
			double den1 = b2c2 * untxx*untxx + a2c2 * untyy*untyy + a2b2*untzz*untzz;
			double den2 = b2c2 * wntxx*wntxx + a2c2 * wntyy*wntyy + a2b2*wntzz*wntzz;
			at1 = a2b2c2 / den1;
			at2 = a2b2c2 / den2;
		}

		
		UN[J[0]] = untxx;
		UN[J[1]] = untyy;
		UN[J[2]] = untzz;
		double unxx = UN[0];
		double unyy = UN[1];
		double unzz = UN[2];
		WN[J[0]] = wntxx;
		WN[J[1]] = wntyy;
		WN[J[2]] = wntzz;
		double wnxx = WN[0];
		double wnyy = WN[1];
		double wnzz = WN[2];
		ROTATE(&vec[0], &unxx, &unyy, &unzz, &unvec[0]);
		ROTATE(&vec[0], &wnxx, &wnyy, &wnzz, &wnvec[0]);
		unx = unvec[0];
		uny = unvec[1];
		unz = unvec[2];
		wnx = wnvec[0];
		wny = wnvec[1];
		wnz = wnvec[2];
	}
	dt1 = at1* *vmag;
	dt2 = at2* *vmag;

	double rMat[9] = { vnx, unx, wnx, vny, uny, wny, vnz, unz, wnz };
	//rotmat.assign({ vnx, unx, wnx, vny, uny, wny, vnz, unz, wnz });
	TENSYM(&DL, &dt1, &dt2,&rMat[0], dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz);
}

void Storage::ROTATE(std::vector<double>& vec, double& v1, double& v2, double& v3, std::vector<double>& out_vec)
{
	out_vec[0] = vec[0] * v1 + vec[1] * v2 + vec[2] * v3;
	out_vec[1] = vec[3] * v1 + vec[4] * v2 + vec[5] * v3;
	out_vec[2] = vec[6] * v1 + vec[7] * v2 + vec[8] * v3;
}

void Storage::ROTATE(double * vec, double * v1, double* v2, double* v3, double * out_vec)
{
	out_vec[0] = vec[0] * *v1 + vec[1] * *v2 + vec[2] * *v3;
	out_vec[1] = vec[3] * *v1 + vec[4] * *v2 + vec[5] * *v3;
	out_vec[2] = vec[6] * *v1 + vec[7] * *v2 + vec[8] * *v3;
}

void Storage::GLOCOL(int L, double vole[], double bflowe[8][8], double dflowe[], double btrane[8][8], double dtrane[8][8])
{
	int n1 = L*N48;
	int n8 = n1 + N48;

	if (ML == 2)
		goto u_only;
	int ie = 0;
	int m = -1;
	for (int i = n1; i < n8; i++)
	{
		
		int ib = incidence_vector[i];
		ib = ib - 1;

		node_vol[ib] = node_vol[ib] + vole[ie];
		node_p_rhs[ib] = node_p_rhs[ib] + dflowe[ie];
		int je = 0;
		for (int j = n1; j < n8; j++)
		{
			int jb = incidence_vector[j];
			jb = jb - 1;
			int mbeg = JA[jb];
			int mend = JA[jb + 1];
			bool found = false;
			for (int mm = mbeg; mm < mend; mm++)
			{
				if (ib == IA[mm])
				{
					m = mm;
					found = true;
					break;
				}
			}
			if (found)
			{
				if (m != -1)
					PMAT[m] = PMAT[m] + bflowe[ie][je];
				else
					SimulationControl::exitOnError(" m negative");
			}
			
			je++;
		}
		
	

		ie++;
	}

	if (ML == 1)
		return;

	u_only:
	if (NOUMAT != 1)
	{
		ie = 0;
		for (int i = n1; i < n8; i++)
		{
			int ib = incidence_vector[i];
			ib = ib - 1;
			int je = 0;
			for (int j = n1; j < n8; j++)
			{
				int jb = incidence_vector[j];
				jb = jb - 1;
				int mbeg = JA[jb];
				int mend = JA[jb + 1];
				bool found = false;
				for (int mm = mbeg; mm < mend; mm++)
				{
					if (ib == IA[mm])
					{
						m = mm;
						found = true;
						break;
					}
				}
				if (found)
				{
					if (m != -1)
						UMAT[m] = UMAT[m] + dtrane[ie][je] + btrane[ie][je];
					else
						SimulationControl::exitOnError(" m negative");
				}
				je++;
			}
			
			ie++;
		}
	}
}

double Storage::DNRM2(int N, double * X, int INCX)
{
	double NORM, SSQ, SCALE;
	double ABSXI;
	NORM = 0.0;
	if (N < 1 || INCX < 1)
	{
		NORM = 0.0;
	}
	else if (N == 1)
	{
		NORM = abs(X[0]);
	}
	else
	{
		SCALE = 0.0;
		SSQ = 1.0;

		for (int it = 1; it < (1 + (N - 1)*INCX); it += INCX)
		{
			if (X[it] != 0.0)
			{
				ABSXI = abs(X[it]);
				if (SCALE < ABSXI)
				{
					SSQ = 1.0 + SSQ * pow(SCALE / ABSXI, 2);
					SCALE = ABSXI;
				}
				else
				{
					SSQ = SSQ + pow(ABSXI / SCALE, 2);
				}
			}
		}
		NORM = SCALE * sqrt(SSQ);
	}
	return NORM;
}

void Storage::re_orient_matrix(int jmper_size, int vals_size, double vals[], std::vector<int>&jmper, std::vector<int>& indices, double * new_vals, int * new_jmper, int * new_indices)
{
	int * rr = new int[vals_size];
	int * nn = new int[jmper_size];
	for (int k = 0, i = 0; i < jmper_size-1; i++)
		for (int j = 0; j < jmper[i + 1] - jmper[i]; j++)
			rr[k++] = i;

	for (int i = 0; i < vals_size; i++)
		new_jmper[indices[i] + 1]++;

	for (int i = 1; i <= jmper_size - 1; i++)
		new_jmper[i] += new_jmper[i - 1];

	memcpy(nn, new_jmper, sizeof(int)*jmper_size);

	for (int i = 0; i < vals_size; i++) {
		int x = nn[indices[i]]++;
		new_vals[x] = vals[i];
		new_indices[x] = rr[i];
	}
	delete[]rr;
	delete[]nn;
}


void Storage::banner(){
	std::string logLine = "";
	Writer * logWriter = Writer::instance("LST");
	std::string lstFile;
	lstFile.append(InputFiles::instance()->getInputDirectory());
	lstFile.append(InputFiles::instance()->getFilesForWriting()["LST"]);
	logWriter->set_filename(lstFile);
	for (int i = 0; i < 133; i++){
		logLine.append("*");
	}
	logLine.append("\n\n\n");
	for (int i = 0; i < 2; i++){
		logLine.append(logLine);
	}

	logLine.append("\n\n\n\n\n");
	logLine.append(std::string(37, ' ') + "II   PPPPPPP     SSS            SSS     II   MMMM      MMMM\n");
	logLine.append(std::string(37, ' ') + "II   PP    PP  SS   SS        SS   SS   II   MM  MM  MM  MM\n");
	logLine.append(std::string(37, ' ') + "II   PP    PP  SS             SS        II   MM    MM    MM\n");
	logLine.append(std::string(37, ' ') + "II   PPPPPPP    SSSSS   ====   SSSSS    II   MM    MM    MM\n");
	logLine.append(std::string(37, ' ') + "II   PP             SS             SS   II   MM          MM\n");
	logLine.append(std::string(37, ' ') + "II   PP        SS   SS        SS   SS   II   MM          MM\n");
	logLine.append(std::string(37, ' ') + "II   PP          SSS            SSS     II   MM          MM\n");
	logLine.append("\n\n\n\n\n\n");
	logLine.append(setLSTLine("N O R T H E A S T E R N   U N I V E R S I T Y"));
	logLine.append("\n\n\n");
	logLine.append(setLSTLine("INDUCED PARTIAL SATURATION SIMULATION"));
	logLine.append("\n\n");
	logLine.append(setLSTLine("-IPSSIM VERSION v1.0-"));
	logLine.append("\n\n\n\n\n\n");
	logWriter->add_line(logLine);
	logLine.clear();
	for (int j = 0; j < 3; j++){
		for (int i = 0; i < 133; i++){
			logLine.append("*");
		}
		logLine.append("\n\n\n");
	}
	for (int i = 0; i < 133; i++){
		logLine.append("*");
	}
	logLine.append("\n\n");

	logWriter->add_line(logLine);
}

std::string Storage::setLSTLine(std::string s){
	int sLen = s.length();
	std::string tempStr1 = "";
	std::string tempStr2 = "";
	std::string tempStr = "";
	int p1, p2;
	p1 = p2 = 0;
	if ((133 - sLen) % 2 == 0){
		p1 = p2 = (133 - sLen) / 2;
	}
	else{
		p1 = (132 - sLen) / 2;
		p2 = p1 + 1;
	}
	tempStr1 = std::string(p1, ' ');
	tempStr2 = std::string(p2, ' ');
	tempStr.append(tempStr1 + s + tempStr2);
	return tempStr;

}


void Storage::outLST()
{
	Writer * lstWriter = Writer::instance("LST");
	char buff[1024];
	std::string logLine;
	if (KTYPE[0] == 3){ //3D
		if (ITREL > 0 || ISSFLO == 2 || ISSTRA == 1)
		{

		} else
		{
			logLine.append("\n\n\n\n");
			logLine.append(std::string(11, ' '));
			logLine.append("I N I T I A L   C O N D I T I O N S\n");
			logLine.append(std::string(11, ' '));
			logLine.append("___________________________________\n");
			if (IREAD == -1)
			{
				logLine.append("\n\n");
				logLine.append(std::string(11, ' '));
				logLine.append("INITIAL CONDITIONS RETRIVED FROM A RESTART FILE (WARM START)\n");
				logLine.append(std::string(11, ' '));
				logLine.append("THAT WAS SAVED AT THE END OF TIME STEP");
				logLine.append(std::to_string(ITRST));
				logLine.append("OF THE ORIGINAL SIMULATION.\n");
			}
			if (IT == 0 && ISSFLO == 2)
			{
				if (KPANDS == 1)
				{
					logLine.append("\n\n\n\n           S T E A D Y - S T A T E   P R E S S U R E\n\n  ");
					logLine.append("      NODE                      NODE                      NODE                      NODE                      NODE                \n");
					int ctr = 0;
					for (int j = 0; j < NN; j++)
					{
						_snprintf(buff, 1024, "   %9d %15.8e", j + 1, node_swt[j]);
						logLine.append(buff);
						ctr++;
						if (ctr == 5 || j == NN - 1){
							logLine.append("\n");
							ctr = 0;
						}
					}
				}
				return;
			} else
			{
				if (ISSTRA == 1)
				{
					if (KCORT == 1)
					{
						if (ME <= 0)
						{
							logLine.append("\n\n\n\n           S T E A D Y - S T A T E   C O N C E N T R A T I O N\n\n  ");
							logLine.append("      NODE                      NODE                      NODE                      NODE                      NODE                \n");
							int ctr = 0;
							for (int j = 0; j < NN; j++)
							{
								_snprintf(buff, 1024, "   %9d %15.8e", j + 1, node_uvec[j]);
								logLine.append(buff);
								ctr++;
								if (ctr == 5 || j == NN - 1){
									logLine.append("\n");
									ctr = 0;
								}
							}
						}

					}
				} else
				{
					logLine.append("\n\n\n");
					logLine.append(std::string(11, ' '));
					_snprintf(buff, 1024, "TIME INCREMENT : %+15.4e SECONDS\n\n", DELT);
					logLine.append(buff);
					logLine.append(std::string(11, ' '));
					_snprintf(buff, 1024, "TIME AT END : %+15.4e SECONDS\n", TSEC);
					logLine.append(buff);
					logLine.append(std::string(11, ' '));
					_snprintf(buff, 1024, "OF STEP : %+15.4e MINUTES\n", TMIN);
					logLine.append(buff);
					logLine.append(std::string(11, ' '));
					_snprintf(buff, 1024, "OF STEP : %+15.4e MINUTES\n", TMIN);
					logLine.append(buff);
				}
			}

		}
	} else
	{ // 2D
		
	}

}
void Storage::outELE()
{
	std::string logLine = "";
	Writer * logWriter = Writer::instance("ELE");
	char buff[1024];
	int TS; // Time Step Information
	int JT; // Time Step Value
	int LCP,LCV;
	char CPVX, CPVY, CPVZ;
	double DELTK;
	std::vector<double> TT;
	std::vector<int> ITT;
	std::vector<int> ISVEL;
	LCOLPR = element_output_every;
	if (onceELE == false)
	{
		std::string eleFile;
		eleFile.append(InputFiles::instance()->getInputDirectory());
		eleFile.append(InputFiles::instance()->getFilesForWriting()["ELE"]);
		logWriter->set_filename(eleFile);
		KT = 0;
		for (int jt = 1; jt <= ITMAX; jt++)
		{
			if (jt%LCOLPR == 0 || jt == ITRST + 1)
				KT = KT + 1;
		}

		if (ITMAX > 1 && ITMAX%LCOLPR != 0)
			KT = KT + 1;

		KTMAX = KT;

		TS = TSTART;
		JT = 0;
		KT = 0;
		CPVX = CPVY = CPVZ = 'N';
		for (int i = 0; i < 7; i++){
			if (J6COL[i] == 4) CPVX = 'Y';
			if (J6COL[i] == 5) CPVY = 'Y';
			if (J6COL[i] == 6) CPVZ = 'Y';
		}
		LCP = 0;
		for (int i = 0; i <= ITMAX; i++){
			TS = schedule_list[time_steps_index]->get_time_at_step(i);
			JT = schedule_list[time_steps_index]->get_step_time()[i].first;
			LCV = LCP;
			if ((JT%NPCYC == 0) || (BCSFL[JT] || JT == 1))
				LCP = JT;
			if ((JT%LCOLPR == 0) || (JT == (ITRST + 1))){
				KT = KT + 1;
				TT.push_back(TS);
				ITT.push_back(JT);
				ISVEL.push_back(LCV);
				if (JT != 1 || ISSFLO == 2)
					if (JT %NUCYC != 0 && !BCSTR[JT])
						ISVEL[KT] = 0;
			}
		}
		if (ISSTRA == 1)
			TT[KT] = TSTART;

		if (ITMAX > 1 && ITMAX%LCOLPR != 0)
		{
			KT = KT + 1;
			TT.push_back(TS);
			ITT.push_back(ITMAX);
			ISVEL.push_back(LCV);
		}
		if (ISSFLO != 0)
		{
			KTMAX = 1;
			ISVEL[0] = 0;
		}

		if (IREAD == 1)
			KTPRN = KTMAX;
		else
			KTPRN = 0;

		for (int i = 0; i < KTMAX; i++)
		{
			if (ITT[i] > ITRST)
				KTPRN = KTPRN + 1;
		}

		logLine.append("## " + titles[0] + "\n");
		logLine.append("## " + titles[1] + "\n");
		logLine.append("## \n");
		std::string CTYPE2;
		if (KTYPE[1] > 1){

			if (KTYPE[1] == 3) { CTYPE2 = "BLOCKWISE MESH"; }
			else{ CTYPE2 = "REGULAR MESH"; }

			if (KTYPE[0] == 3){ // 3D
				_snprintf(buff, sizeof(buff), "## %1d-D,", KTYPE[0]);
				logLine.append(buff + CTYPE2);
				_snprintf(buff, sizeof(buff), "  (%9d)*(%9d)*(%9d) = %9d Nodes ( %9d Elems)\n", NN1, NN2, NN3, NN, NE);
				logLine.append(buff);
				logLine.append("## \n");
			}
			else{ // 2D
				_snprintf(buff, sizeof(buff), "## %1d-D,", KTYPE[0]);
				logLine.append(buff + CTYPE2);
				_snprintf(buff, sizeof(buff), "  (%9d)*(%9d) = %9d Nodes ( %9d Elems)\n", NN1, NN2, NN, NE);
				logLine.append(buff);
				logLine.append("## \n");

			}

		}
		else if (KTYPE[1] == 1){
			_snprintf(buff, sizeof(buff), "## %1d-D, LAYERED MESH [", KTYPE[0]);
			logLine.append(buff);
			logLine.append(LAYSTR + "]");
			_snprintf(buff, sizeof(buff), "       (%9d)*(%9d) = %9d Nodes ( %9d Elems)\n", NLAYS, NNLAY, NN, NE);
			logLine.append(buff);
			logLine.append("## \n");
		}
		else{
			_snprintf(buff, sizeof(buff), "## % 1d - D, IRREGULAR MESH", KTYPE[0]);
			logLine.append(buff);
			logLine.append(std::string(40, ' '));
			_snprintf(buff, sizeof(buff), "%9d Nodes ( %9d Elems)\n", NN, NE);
			logLine.append(buff);
			logLine.append("## \n");
		}

		logLine.append("## " + std::string(92, '=') + "\n## VELOCITY RESULTS" + std::string(48, ' '));
		_snprintf(buff, 1024, "%9d Time steps printed\n", KTPRN);
		logLine.append(buff);
		logLine.append("## " + std::string(92, '=') + "\n## \n");
		logLine.append("##     Time steps" + std::string(22, ' ') + "[Printed? / Time step on which V is based]\n");
		logLine.append("##    in this file      Time (sec)          VX             VY             VZ\n");
		logLine.append("##  " + std::string(14, '-') + "   " + std::string(13, '-') + "    " + std::string(12, '-') + "   " + std::string(12, '-') + "   " + std::string(12, '-')+"\n");
		for (int i = 1; i <= KTMAX; i++){
			if (ITT[i] >= ITRST){
				_snprintf(buff, sizeof(buff), "##        %8d    %+13.6e      %c %8d     %c %8d     %c %8d\n", ITT[i], TT[i],CPVX, ISVEL[i],CPVY,ISVEL[i],CPVZ,ISVEL[i]);
				logLine.append(buff);
			}
		}
		onceELE = true;
	}

	if ((ISSFLO == 2) && (IT > 1))
		return;

	if (IT == 1 && ISSTRA == 1)
	{
		DURN = 0.0;
		TOUT = TSTART;
	}
	else
	{
		DURN = DELT;
		TOUT = TSEC;
	}





	logLine.append("## \n## " + std::string(92, '=') + "\n ## TIME STEP ");
	_snprintf(buff, 1024, "%8d                       Duration: %11.4e sec      Time: %11.4e sec \n ## ", IT, DURN, TOUT);
	logLine.append(buff);
	logLine.append(std::string(92, '=') +"\n");


	logLine.append("##");
	for (std::string a : LCOL){
		int fil = 15 - a.length();
		logLine.append(std::string(fil, ' ') + a + "  ");
	}

	logLine.append("\n");
	double centrx, centry, centrz;
	for (int i = 0; i < NE; i++){
		centrx = centry = centrz = 0;
		for (int v : incidenceContainer[i].second)
		{
			centrx = centrx + node_x[v - 1];
			centry = centry + node_y[v - 1];
			centrz = centrz + node_z[v - 1];
		}
		centrx = centrx * 1 / N48;
		centry = centry * 1 / N48;
		centrz = centrz * 1 / N48;
		double va1 = 0.017453292*el_vang1[i];
		int ll = min(i + 1, NEX);
		double va2 = 0.017453292*el_vang2[ll - 1] * (KTYPE[0]-2);
		double cva = cos(va2);

		double vectrx = el_vmag[i] * cos(va1)*cva;
		double vectry = el_vmag[i] * sin(va1)*cva;
		double vectrz = el_vmag[i] * sin(va2);
		_snprintf(buff, sizeof(buff), "  %+14.7e  %+14.7e  %+14.7e  %+14.7e  %+14.7e  %+14.7e\n", centrx,centry,centrz,vectrx,vectry,vectrz);
		logLine.append(buff);
	}


	logWriter->add_line(logLine);






}
void Storage::outNOD()
{
	std::string logLine = "";
	Writer * logWriter = Writer::instance("NOD");
	char buff[512];
	
	int JT; // Time Step Value
	int TS; // Time Step Information
	double DELTK; // time sstep increment
	std::vector<double> TT;
	std::vector<int> ITT;
	std::vector<int> ISHORP, ISTORC, ISSATU;
	int LCHORP, LCTORC;
	char CPSATU, CPTORC, CPHORP;
	NCOLPR = node_output_every;
	std::string header = "\nNode              X              Y              Z       Pressure  Concentration     Saturation\n";
	if (onceNOD == false){

		std::string nodFile;
		nodFile.append(InputFiles::instance()->getInputDirectory());
		nodFile.append(InputFiles::instance()->getFilesForWriting()["NOD"]);
		logWriter->set_filename(nodFile);
		if (ISSTRA != 1){
			KT = 1;
		}
		else{
			KT = 0;
		}
		for (int i = 1; i <= ITMAX; i++){
			if (((i%NCOLPR) == 0) || (i == ITRST) || ((i == (ITRST + 1)) && ((ISSTRA != 0) || (NCOLPR >0)))){
				KT = KT + 1;
			}
		}
		if (ITMAX > 1 && (ITMAX%NCOLPR != 0))
			KT = KT + 1;

		KTMAX = KT;
		TS = TSTART;
		JT = 0;
		KT = 0;
		DELTK = DELT;



		// pressure conc sat print y or no
		// KPANDS pressure and sat
		if (KPANDS == 1){
			CPHORP = CPSATU = 'Y';
		}
		else{ CPHORP = CPSATU = 'N'; }
		if (KCORT == 1){
			CPTORC = 'Y';

		}
		else{ CPTORC = 'N'; }
		//KCORD concentration

		if (ISSTRA != 1){
			KT = KT + 1;
			TT.push_back(TS);
			ITT.push_back(JT);
			ISHORP.push_back(0);
			ISTORC.push_back(0);
			ISSATU.push_back(0);
		}

		for (int i = 0; i <= ITMAX; i++){
			TS = schedule_list[time_steps_index]->get_time_at_step(i);
			JT = schedule_list[time_steps_index]->get_step_time()[i].first;

			if ((JT%NPCYC == 0) || (BCSFL[JT] || JT == 1))
				LCHORP = JT;
			if ((JT%NUCYC == 0) || (BCSTR[JT] || JT == 1))
				LCTORC = JT;
			if ((JT%NCOLPR == 0) || (JT == ITRST) || ((JT == (ITRST + 1)) && ((ISSTRA != 0) || NCOLPR >0))){
				KT = KT + 1;
				TT.push_back(TS);
				ITT.push_back(JT);
				ISHORP.push_back(LCHORP);
				ISTORC.push_back(LCTORC);
				ISSATU.push_back(LCHORP);
			}
		}

		if (ISSTRA == 1){
			TT.push_back(TSTART);
			ITT.push_back(0);
		}

		if (ITMAX > 1 && (ITMAX%NCOLPR != 0)){
			KT = KT + 1;
			TS = schedule_list[time_steps_index]->get_time_at_step(ITMAX);
			TT.push_back(TS);
			ITT.push_back(ITMAX);
			if ((ITMAX%NPCYC == 0) || (BCSFL[ITMAX]))
				LCHORP = ITMAX;
			if ((ITMAX%NPCYC == 0) || (BCSTR[ITMAX]))
				LCTORC = ITMAX;
			ISHORP.push_back(LCHORP);
			ISTORC.push_back(LCTORC);
			ISSATU.push_back(LCHORP);
		}

		if (ISSFLO != 0){
			for (int i = 0; i < KTMAX; i++){
				ISHORP.push_back(0);
				ISSATU.push_back(0);
			}
		}

		if (IREAD == 1){ KTPRN = KTMAX; }
		else{
			KTPRN = 0;
			for (int i = 1; i <= KTMAX; i++){
				if (ITT[i] > ITRST)
					KTPRN = KTPRN + 1;
			}
		}

		logLine.append("## " + titles[0] + "\n");
		logLine.append("## " + titles[1]+ "\n");
		logLine.append("## \n");

		std::string CTYPE2;
		if (KTYPE[1] > 1){

			if (KTYPE[1] == 3) { CTYPE2 = "BLOCKWISE MESH"; }
			else{ CTYPE2 = "REGULAR MESH"; }

			if (KTYPE[0] == 3){ // 3D
				_snprintf(buff, sizeof(buff), "## %1d-D,", KTYPE[0]);
				logLine.append(buff + CTYPE2);
				_snprintf(buff, sizeof(buff), "  (%9d)*(%9d)*(%9d) = %9d Nodes ( %9d Elems)\n", NN1, NN2, NN3, NN, NE);
				logLine.append(buff);
				logLine.append("## \n");
			}
			else{ // 2D
				_snprintf(buff, sizeof(buff), "## %1d-D,", KTYPE[0]);
				logLine.append(buff + CTYPE2);
				_snprintf(buff, sizeof(buff), "  (%9d)*(%9d) = %9d Nodes ( %9d Elems)\n", NN1, NN2, NN, NE);
				logLine.append(buff);
				logLine.append("## \n");

			}

		}
		else if (KTYPE[1] == 1){
			_snprintf(buff, sizeof(buff), "## %1d-D, LAYERED MESH [", KTYPE[0]);
			logLine.append(buff);
			logLine.append(LAYSTR + "]");
			_snprintf(buff, sizeof(buff), "       (%9d)*(%9d) = %9d Nodes ( %9d Elems)\n", NLAYS, NNLAY, NN, NE);
			logLine.append(buff);
			logLine.append("## \n");
		}
		else{
			_snprintf(buff, sizeof(buff), "## % 1d - D, IRREGULAR MESH", KTYPE[0]);
			logLine.append(buff);
			logLine.append(std::string(40, ' '));
			_snprintf(buff, sizeof(buff), "%9d Nodes ( %9d Elems)\n", NN, NE);
			logLine.append(buff);
			logLine.append("## \n");
		}


		logLine.append("## " + std::string(92, '=') + "\n## NODEWISE RESULTS" + std::string(48, ' '));
		_snprintf(buff, sizeof(buff), "%9d Time steps printed\n", KTPRN);
		logLine.append(buff);
		logLine.append("## " + std::string(92, '=') + "\n## \n");
		logLine.append("##    Time steps" + std::string(23, ' ') + "[Printed? / Latest time step computed]\n");
		logLine.append("##    in this file      Time (sec)         Press          Conc           Sat\n");
		logLine.append("##   " + std::string(14, '-') + "   " + std::string(13, '-') + "    " + std::string(12, '-') + "   " + std::string(12, '-') + "   " + std::string(12, '-') + "\n");



		for (int i = 1; i <= KTMAX; i++){
			if (ITT[i] >= ITRST){
				_snprintf(buff, sizeof(buff), "##        %8d    %+13.6e      %c %8d     %c %8d     %c %8d\n", ITT[i], TT[i], CPHORP, ISHORP[i], CPTORC, ISTORC[i], CPSATU, ISSATU[i]);
				logLine.append(buff);
			}
		}



		onceNOD = true;
	}
	if ((IT == 0) || ((IT == 1) && (ISSTRA == 1))){
		DURN = 0.0;
		TOUT = TSTART;
		ITOUT = 0;
	}
	else{
		DURN = DELT;
		TOUT = TSEC;
		ITOUT = IT;
	}

	logLine.append("## \n");
	logLine.append("## " + std::string(98, '=') + "\n");
	_snprintf(buff, sizeof(buff), "## TIME STEP %8d", ITOUT);
	logLine.append(buff + std::string(26, ' '));
	_snprintf(buff, sizeof(buff), "Duration: %+11.4e sec      Time: %+11.4e sec\n", DURN, TOUT);
	logLine.append(buff);
	logLine.append("## " + std::string(98, '=') + "\n");




	logLine.append("##");
	for (std::string a : NCOL){
		int fil = 15 - a.length();
		logLine.append(std::string(fil, ' ') + a + "  ");
	}

	logLine.append("\n");

	for (int i = 0; i < NN; i++){

		_snprintf(buff, sizeof(buff), "  %+14.7e  %+14.7e  %+14.7e  %+14.7e  %+14.7e  %+14.7e\n", node_x[i], node_y[i], node_z[i], node_pvec[i], node_uvec[i], node_swt[i]);
		logLine.append(buff);
	}


	logWriter->add_line(logLine);



}
void Storage::outOBS()
{

}

/*	std::ofstream outbin("p_rhs.bin", std::ios::binary);
for (int i = 0; i < NN; i++)
outbin.write(reinterpret_cast<const char*>(&node_p_rhs[i]), sizeof(double));
outbin.close();*/


/*
*

for (int i = 0; i < 16621000;i++)
d[i] = i+0.123456789012345;
std::ofstream outbin("C:/Users/demiryurek.a/Desktop/pvec_uvec/p_rhs.bin", std::ios::binary);
for (int i = 0; i < 16621000; i++){


outbin.write(reinterpret_cast<const char *>(&d[0]+i), sizeof(double));
}
outbin.close();

std::streampos size;
char * memblock;
std::ifstream inbin("C:/Users/demiryurek.a/Desktop/pvec_uvec/p_rhs.bin", std::ios::in | std::ios::binary | std::ios::ate);

if (inbin.is_open())
{
size = inbin.tellg();
memblock = new char[size];
inbin.seekg(0, std::ios::beg);
inbin.read(memblock, size);
inbin.close();

std::cout << "the entire file content is in memory";
double * double_values = (double*)memblock;

delete[] memblock;
}
else
std::cout << "Unable to open file";

*/