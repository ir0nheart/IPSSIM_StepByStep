#ifndef ELEMENT_H
#define ELEMENT_H
#pragma once
class Element
{
public:
	Element(int &el_num, int &lreg, double &pmax, double &pmid, double &pmin, double &ang1, double &ang2, double &ang3, double &almax,
		double &almid,double &almin,double &atmax,double &atmid ,double &atmin);
	void ROTMAT();
	void TENSYM();
	void reserve_GXSI(int val){ GXSI.reserve(val); }
	void reserve_GETA(int val){ GETA.reserve(val); }
	void reserve_GZET(int val){ GZET.reserve(val); }
	double get_GXLOC(int val) const{ return GXLOC[val]; }
	double get_GYLOC(int val) const{ return GYLOC[val]; }
	double get_GZLOC(int val) const{ return GZLOC[val]; }

	void set_jacobi_det(int ind,double val){DET[ind] = val; }
	void set_gxsi_geta_gzet(int ind, std::vector<double>& jacobi, double GRAVX, double GRAVY, double GRAVZ);
	
	~Element();
private:
	double GLOC;
	int INTIM;
	int ISTOP;
	std::vector<double> GXLOC;
	std::vector<double> GYLOC;
	std::vector<double> GZLOC;
	std::vector<std::vector<double>> rotationMatrix;
	double D2R;
	double almax;
	double almin;
	double atmax;
	double atmin;
	double vmag;
	double vang1;
	double vang2;
	double permxx;
	double permxy;
	double permyx;
	double permyy;
	double pangl1;
	double almid;
	double atmid;
	double permxz;
	double permyz;
	double permzz;
	double permzx;
	double permzy;
	double pangl2;
	double pangl3;
	double permmax;
	double permmid;
	double permmin;
	double angle1;
	double angle2;
	double angle3;
	std::vector<double> GXSI;
	std::vector<double> GETA;
	std::vector<double> GZET;
	std::vector<double> DET; // row based {c11,c12,c13,c21,c22,c23,c31,c32,c33}
	
	int lreg;
	int el_num;
};
#endif
