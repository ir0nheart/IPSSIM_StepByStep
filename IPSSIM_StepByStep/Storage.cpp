#include "stdafx.h"
#include "Storage.h"
#include "SimulationControl.h"
#include "obsPoint.h"
#include "Timer.h"

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
	x1 = nodeContainer[incidence_vector[m0]-1].get_x();
	x2 = nodeContainer[incidence_vector[m0+1]-1].get_x();
	x3 = nodeContainer[incidence_vector[m0+2]-1].get_x();
	x4 = nodeContainer[incidence_vector[m0+3]-1].get_x();
	x5 = nodeContainer[incidence_vector[m0+4]-1].get_x();
	x6 = nodeContainer[incidence_vector[m0+5]-1].get_x();
	x7 = nodeContainer[incidence_vector[m0+6]-1].get_x();
	x8 = nodeContainer[incidence_vector[m0+7]-1].get_x();
	y1 = nodeContainer[incidence_vector[m0]-1].get_y();
	y2 = nodeContainer[incidence_vector[m0 + 1]-1].get_y();
	y3 = nodeContainer[incidence_vector[m0 + 2]-1].get_y();
	y4 = nodeContainer[incidence_vector[m0 + 3]-1].get_y();
	y5 = nodeContainer[incidence_vector[m0 + 4]-1].get_y();
	y6 = nodeContainer[incidence_vector[m0 + 5]-1].get_y();
	y7 = nodeContainer[incidence_vector[m0 + 6]-1].get_y();
	y8 = nodeContainer[incidence_vector[m0 + 7]-1].get_y();
	z1 = nodeContainer[incidence_vector[m0]-1].get_z();
	z2 = nodeContainer[incidence_vector[m0 + 1]-1].get_z();
	z3 = nodeContainer[incidence_vector[m0 + 2]-1].get_z();
	z4 = nodeContainer[incidence_vector[m0 + 3]-1].get_z();
	z5 = nodeContainer[incidence_vector[m0 + 4]-1].get_z();
	z6 = nodeContainer[incidence_vector[m0 + 5]-1].get_z();
	z7 = nodeContainer[incidence_vector[m0 + 6]-1].get_z();
	z8 = nodeContainer[incidence_vector[m0 + 7]-1].get_z();


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
				if (neighbor_ == nod_)
				{
					continue;
				} else
				{
					nodeContainer[nod_ - 1].add_neighbor(neighbor_);
				}
			}
		}
		
	}
	NELT = 0;
	for (int i = 0; i < nodeContainer.size(); i++)
	{
		nodeContainer[i].make_neighbors_unique();
		NELT += nodeContainer[i].get_neighbor_size();
	}
	NELT = NELT + NN;

	IA.reserve(NELT);
	JA.reserve(NDIMJA);

	for (int i = 0; i < nodeContainer.size(); i++)
	{
		JA.push_back(IA.size());
		IA.push_back(i + 1);
		for (int k = 0; k < nodeContainer[i].get_neighbor_size(); k++)
		{
			IA.push_back(nodeContainer[i].get_neighbor(k));
		}
	}
	JA.push_back(NELT);
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
	// Check Titles
	if (titles.size() != 2)
	{
		std::cout << "Titles are not properly defined." << std::endl;
		SimulationControl::exitOnError();
	}

	if (!strcmp(sutra_string.data() , "SUTRA"))
	{
		SimulationControl::exitOnError("INP-2A-1");
	}

	if (!strncmp(version_string.data(),"VERSION",7))
	{
		if (version_num_string.data() == "2D3D.1")
		{
			version_num_string = std::vector<char>{'2','.','0'};
		} else if (!strcmp(version_num_string.data(),"2.0") &&
			!strcmp(version_num_string.data(),"2.1") &&
			!strcmp(version_num_string.data(), "2.2"))
		{
			SimulationControl::exitOnError("INP-2A-4");
		}
	} else
	{
		version_num_string = std::vector<char>{'2', '.', '0'};
	}

	if (strcmp(simulation_type_string.data(), "SOLUTE"))
	{
		ME = -1;
		std::cout << "SOLUTE TRANSPORT " << std::endl;
	} else if (strcmp(simulation_type_string.data(), "ENERGY"))
	{
		ME = +1;
		std::cout << "ENERGY TRANSPORT " << std::endl;
	} else
	{
		std::cout << "Simulation TYPE is not Solute or Energy Transport " << std::endl;
		SimulationControl::exitOnError("INP-2A-2");
	}
	
	if (!strncmp(mesh_dim_string.data(), "2D",2))
	{
		KTYPE.push_back(2);
	} else if (!strncmp(mesh_dim_string.data(), "3D",2))
	{
		KTYPE.push_back(3);
	} else
	{
		std::cout << "Error in Mesh type definition." << std::endl;
		std::cout << "\t Mesh type should be either '2D' or '3D'" << std::endl;
		SimulationControl::exitOnError("INP-2B-1");
	}

	// This program runs only REGULAR and IRREGULAR MESH
	if (!strncmp(mesh_type_string.data(), "REGULAR",mesh_type_string.size()))
	{
		KTYPE.push_back(2);
		if (KTYPE[0] == 2)
		{
			NN1 = std::stoi(std::string(nn1_string.begin(), nn1_string.end()));
			NN2 = std::stoi(std::string(nn2_string.begin(), nn2_string.end()));
			NN3 = 1;
		} else if (KTYPE[0] == 3)
		{
			NN1 = std::stoi(std::string(nn1_string.begin(), nn1_string.end()));
			NN2 = std::stoi(std::string(nn2_string.begin(), nn2_string.end()));
			NN3 = std::stoi(std::string(nn3_string.begin(), nn3_string.end()));
		}
		if (NN1 < 2 || NN2 < 2 || (KTYPE[0] == 3 && NN3 < 2))
		{
			SimulationControl::exitOnError("INP-2B-3");
		}

	} else if (!strncmp(mesh_type_string.data(), "IRREGULAR",mesh_type_string.size()))
	{
		KTYPE.push_back(0);
	} else
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
		} else
		{
			NE123 = (NN1 - 1)*(NN2 - 1);
		}
		if (NE123 != NE)
		{
			SimulationControl::exitOnError("INP-2B,3-2");
		}
	}

	if (!strncmp(simulation_condition_string.data(), "UNSATURATED",simulation_condition_string.size()))
		IUNSAT = +1;
	else if (!strncmp(simulation_condition_string.data(), "SATURATED", simulation_condition_string.size()))
		IUNSAT = 0;
	else
		SimulationControl::exitOnError("INP-4-1");


	if (!strncmp(flow_type_string.data(), "TRANSIENT",flow_type_string.size()))
		ISSFLO = 0;
	else if (!strncmp(flow_type_string.data(), "STEADY", flow_type_string.size()))
		ISSFLO = +1;
	else
		SimulationControl::exitOnError("INP-4-2");

	if (!strncmp(transport_type_string.data(), "TRANSIENT", transport_type_string.size()))
		ISSTRA = 0;
	else if (!strncmp(transport_type_string.data(), "STEADY", transport_type_string.size()))
		ISSTRA = +1;
	else
		SimulationControl::exitOnError("INP-4-3");

	if (!strncmp(simulation_start_string.data(), "COLD", simulation_start_string.size()))
		IREAD = +1;
	else if (!strncmp(simulation_start_string.data(), "WARM", simulation_start_string.size()))
		IREAD = -1;
	else
		SimulationControl::exitOnError("INP-4-4");

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
		goto _l846;
	}

	NSCH = NSCH + NSCHAU;

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
			std::cout << "bitti" << std::endl;
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

	// Create Nodes
	nodeContainer.reserve(NN);
	if (KTYPE[0] == 3){
		int node_num;
		int nreg;
		double x, y, z, por;
		for (std::vector<char> str : nodeData)
		{
			str.push_back('\0');
			node_num = std::stoi(strtok(str.data(), " "));
			nreg = std::stoi(strtok(NULL, " "));
			x = std::stod(strtok(NULL, " ")) * SCALX;
			y = std::stod(strtok(NULL, " "))*SCALY;
			z = std::stod(strtok(NULL, " "))*SCALZ;
			por = std::stod(strtok(NULL, " "))*PORFAC;
			nodeContainer.push_back(Node(node_num, nreg, x, y, z, por));
		}
	} else
	{
		//2D Nodes
	}




	// define incidence
	incidenceContainer.reserve(NE);
		int el_num;
		for (int j = 1; j < incidenceData.size();j++)
		{
			std::vector<int> nodes;
			el_num = std::stoi(strtok(incidenceData[j].data(), " "));
			for (int k = 0; k < N48; k++)
			{
				nodes.push_back(std::stoi(strtok(NULL, " ")));
			}
			incidenceContainer.push_back(std::pair<int, std::vector<int>>(el_num, nodes));
		}

	

	// define iqsop,iqsou,ipbc,iubc,iidpbc,iidubc,iidsop,iidsou,ibcpbc,ibcubc,ibcsop,ibcsou
	 // !!TODO
	BCSTR.reserve(ITMAX);
	BCSFL.reserve(ITMAX);

	if (simulation_output_controls[0][0] == 'Y')
		KNODAL = +1;
	else
		KNODAL = 0;

	if (simulation_output_controls[1][0] == 'Y')
		KELMNT = +1;
	else
		KELMNT = 0;

	if (simulation_output_controls[2][0] == 'Y')
		KINCID = +1;
	else
		KINCID = 0;

	if (simulation_output_controls[3][0] == 'Y')
		KPANDS = +1;
	else
		KPANDS = 0;

	if (simulation_output_controls[4][0] == 'Y')
		KVEL = +1;
	else
		KVEL = 0;

	if (simulation_output_controls[5][0] == 'Y')
		KCORT = +1;
	else
		KCORT = 0;

	if (simulation_output_controls[6][0] == 'Y')
		KBUDG = +1;
	else
		KBUDG = 0;

	if (simulation_output_controls[7][0] == 'Y')
		KSCRN = +1;
	else
		KSCRN = 0;

	if (simulation_output_controls[8][0] == 'Y')
		KPAUSE = +1;
	else
		KPAUSE = 0;


	// Set NODAL OUTPUT HEADERS
	// !!TODO

	//Set ELE OUTPUT HEADERS
	// !!TODO

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

	for (int j = 0; j < nodeContainer.size();j++)
	{
		nodeContainer[j].set_sop(COMPMA, COMPFL);
		nodeContainer[j].set_qin(0.0);
		nodeContainer[j].set_uin(0.0);
		nodeContainer[j].set_quin(0.0);
		nodeContainer[j].reserve_neighbors(N48*(N48-1));
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
	elementContainer.reserve(NE);

	if (KTYPE[0] == 3){
		int el_num;
		int lreg;
		double pmax, pmid, pmin, ang1, ang2, ang3, almax, almid, almin, atmax, atmin, atmid;
		for (int j = 0; j < elementData.size();j++)
		{
			//str.push_back('\0');
			el_num = std::stoi(strtok(elementData[j].data(), " "));
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
			atmin = std::stod(strtok(NULL, " "))*ATMNF;
			elementContainer.push_back(Element(el_num, lreg, pmax, pmid, pmin, ang1, ang2, ang3, almax, almid, almin, atmax, atmid, atmin));
		}
	}
	else
	{
		//2D Nodes
	}

	NSOPI = NSOP - 1;
	NSOUI = NSOU - 1;
	IQSOPT = 1;
	IQSOUT = 1;
	if (NSOPI != 0)
	{
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
				if (qinc > 0)
					uinc = std::stod(strtok(NULL, " "));
				else
					uinc = 0.0;
			} else
			{
				qinc = 0.0;
				uinc = 0.0;
			}
			IQSOP.push_back(IQCP);
			if (IQCP < 0)
				IQSOPT = -1;
			nodeContainer[IQCPA - 1].set_qin(qinc);
			nodeContainer[IQCPA - 1].set_uin(uinc);
		}

	}

	if (NSOUI != 0)
	{
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
			} else
			{
				quinc = 0;
			}
			IQSOU.push_back(IQCU);
			if (IQCU < 0)
				IQSOUT = -1;

			nodeContainer[IQCUA - 1].set_quin(quinc);
		}
	}

	if (NBCN - 1 >0)
	{
		IPBCT = 1;
		IUBCT = 1;
		if (NPBC != 0)
		{
			int IDUM, IDUMA;
			double pbc_, ubc_;
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
					nodeContainer[IDUMA - 1].set_pbc(pbc_);
					nodeContainer[IDUMA - 1].set_ubc(ubc_);
				} else if (IDUM < 0)
				{
					IPBCT = -1;
				}
			}
			GNUP1 = GNUP;
			IBCPBC = 0;
		}
	
		if (NUBC != 0)
		{
			int IDUM, IDUMA;
			double ubc_;
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
					nodeContainer[IDUMA - 1].set_ubc(ubc_);
				} else if (IDUM < 0)
				{
					IUBCT = -1;
				}
			}
			GNUU1 = GNUU;
			IBCUBC = 0;
		}
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
			for (int nod_ : inc_.second)
				incidence_vector.push_back(nod_);
		}
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


	{
		ITRST = 0;
		if (!strncmp(p_ics_string.data(), "'UNIFORM'", p_ics_string.size()))
		{
			double PUNI = p_ics[0];
			for (int i = 0; i < nodeContainer.size(); i++)
			{
				nodeContainer[i].set_pvec(PUNI);
			}
		}
		else if (!strncmp(p_ics_string.data(), "'NONUNIFORM'", p_ics_string.size()))
		{
			for (int i = 0; i < nodeContainer.size(); i++)
			{
				nodeContainer[i].set_pvec(p_ics[i]);
			}
		} else
		{
			SimulationControl::exitOnError("ICS-2-1");
		}

		if (!strncmp(u_ics_string.data(), "'UNIFORM'", u_ics_string.size()))
		{
			double UUNI = u_ics[0];
			for (int i = 0; i < nodeContainer.size(); i++)
			{
				nodeContainer[i].set_uvec(UUNI);
			}
		}
		else if (!strncmp(u_ics_string.data(), "'NONUNIFORM'", u_ics_string.size()))
		{
			for (int i = 0; i < nodeContainer.size(); i++)
			{
				nodeContainer[i].set_uvec(u_ics[i]);
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
					nodeContainer[IPBC[j] - 1].set_pbc(nodeContainer[IPBC[j] - 1].get_pvec());
			}
		} else
		{
			for (int i = 0; i < NN; i++)
			{
				nodeContainer[i].set_pm1(nodeContainer[i].get_pvec());
				nodeContainer[i].set_um1(nodeContainer[i].get_uvec());
				nodeContainer[i].set_um2(nodeContainer[i].get_uvec());
				nodeContainer[i].set_rcit(RHOW0 + DRWDU*(nodeContainer[i].get_uvec() - URHOW0));
			}
		}
		for (int i = 0; i < NN; i++)
		{
			nodeContainer[i].set_sw(1.0);
			nodeContainer[i].set_swt(1.0);
			nodeContainer[i].set_swb(1.0);
			nodeContainer[i].set_cnub(0.0);
			nodeContainer[i].set_dswdp(0.0);
			nodeContainer[i].set_cs1(CS);
			nodeContainer[i].set_sl(0.0);
			nodeContainer[i].set_sr(0.0);
			nodeContainer[i].set_dpdtitr(0.0);
		}

		if (IUNSAT)
		{
			// if unsaturated call
		}

		TSEC = TSTART;
	}

	std::cout << "" << std::endl;
}
