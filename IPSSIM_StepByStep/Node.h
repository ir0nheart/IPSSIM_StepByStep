#ifndef NODE_H
#define NODE_H
#pragma once
class Node
{
public:
	Node(int& num, int& reg, double& xx, double& yy, double& zz, double& ppor);
	void set_sop(double compma, double compfl){ sop = (1.0 - por) * compma + por*compfl; }
	void change_sop(double sop_){ sop = sop_; }
	void set_quin(double val){ quin = val; }
	void set_uin(double val){ uin = val; }
	void set_qin(double val){ qin = val; }
	void set_pbc(double val){ pbc = val; }
	void set_ubc(double val){ ubc = val; }
	void set_pvec(double val){ pvec = val; }
	void set_uvec(double val){ uvec = val; }
	void set_pm1(double val){ pm1 = val; }
	void set_um1(double val){ um1 = val; }
	void set_um2(double val){ um2 = val; }
	void set_rcit(double val){ rcit = val; }
	void set_sw(double val){ sw = val; }
	void set_swt(double val){ swt = val; }
	void set_swb(double val){ swb = val; }
	void set_cnub(double val){ cnub = val; }
	void set_dswdp(double val){ dswdp = val; }
	void set_cs1(double val){ cs1 = val; }
	void set_sl(double val){ sl = val; }
	void set_sr(double val){ sr = val; }
	void set_dpdtitr(double val){ dpdtitr = val; }
	
	double get_pvec() const{ return pvec; }
	double get_uvec() const{ return uvec; }
	double get_x() const{ return x; }
	double get_y() const{ return y; }
	double get_z() const{ return z; }
	void reserve_neighbors(int val){ neighbors.reserve(val); }
	void add_neighbor(int val){ neighbors.push_back(val); }
	void make_neighbors_unique();
	int get_neighbor_size()const { return neighbors.size(); }
	int get_neighbor(int val) const { return neighbors[val]; }
	~Node();
private:
	int node_num;
	double pvec;
	double uvec;
	double piter;
	double uiter;
	double pm1;
	double um1;
	double dpdtitr;
	double um2;
	double pvel;
	double sl;
	double sr;
	double x;
	double y;
	double z;
	double vol;
	double por;
	double cs1;
	double cs2;
	double cs3;
	double sw;
	double dswdp;
	double rho;
	double sop;
	double qin;
	double uin;
	double quin;
	double qinitr;
	double rcit;
	double rcitm1;
	double swt;
	double cnub;
	double cnumb1;
	double swb;
	double relk;
	double relkb;
	double relkt;
	double effstr;
	double runod;
	double totstr;
	double pbc;
	double ubc;
	int nreg;
	std::vector<int> neighbors;
};
#endif
