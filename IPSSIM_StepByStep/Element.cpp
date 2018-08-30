#include "stdafx.h"
#include "Element.h"



Element::Element(int& el_num, int& lreg, double& pmax, double& pmid, double& pmin, double &ang1, double &ang2, double &ang3,
	double& almax, double& almid,double& almin, double& atmax, double& atmid, double& atmin):el_num(el_num),lreg(lreg),permmax(pmax),permmid(pmid),permmin(pmin),angle1(ang1),
	angle2(ang2), angle3(ang3), almax(almax), atmax(atmax),almid(almid),almin(almin),atmid(atmid),atmin(atmin)
{
	D2R = 0.01745329252;
	pangl1 = D2R * angle1;
	pangl2 = D2R * angle2;
	pangl3 = D2R * angle3;
	ROTMAT();
	GLOC = 0.577350269189626;
	INTIM = ISTOP = 0;
	GXSI.assign({ 0, 0, 0, 0, 0, 0, 0, 0 });
	GETA.assign({ 0, 0, 0, 0, 0, 0, 0, 0 });
	GZET.assign({ 0, 0, 0, 0, 0, 0, 0, 0 });
	DET.assign({ 0, 0, 0, 0, 0, 0, 0, 0 });
	GXLOC.assign({-1,1,1,-1,-1,1,1,-1});
	GYLOC.assign({-1,-1,1,1,-1,-1,1,1});
	GZLOC.assign({-1,-1,-1,-1,1,1,1,1});
}

void Element::TENSYM()
{
	permxx = rotationMatrix[0][0] * rotationMatrix[0][0] * permmax +
		rotationMatrix[0][1] * rotationMatrix[0][1] * permmid +
		rotationMatrix[0][2] * rotationMatrix[0][2] * permmin;
	permxy = rotationMatrix[0][0] * rotationMatrix[1][0] * permmax +
		rotationMatrix[0][1] * rotationMatrix[1][1] * permmid +
		rotationMatrix[0][2] * rotationMatrix[1][2] * permmin;
	permxz = rotationMatrix[0][0] * rotationMatrix[2][0] * permmax +
		rotationMatrix[0][1] * rotationMatrix[2][1] * permmid +
		rotationMatrix[0][2] * rotationMatrix[2][2] * permmin;
	permyy = rotationMatrix[1][0] * rotationMatrix[1][0] * permmax +
		rotationMatrix[1][1] * rotationMatrix[1][1] * permmid +
		rotationMatrix[1][2] * rotationMatrix[1][2] * permmin;
	permyz = rotationMatrix[1][0] * rotationMatrix[2][0] * permmax +
		rotationMatrix[1][1] * rotationMatrix[2][1] * permmid +
		rotationMatrix[1][2] * rotationMatrix[2][2] * permmin;
	permzz = rotationMatrix[2][0] * rotationMatrix[2][0] * permmax +
		rotationMatrix[2][1] * rotationMatrix[2][1] * permmid +
		rotationMatrix[2][2] * rotationMatrix[2][2] * permmin;
	permyx = permxy;
	permzx = permxz;
	permzy = permyz;






}
void Element::ROTMAT()
{
	double s1, s2, s3, c1, c2, c3;
	s1 = s2 = s3 = c1 = c2 = c3 = 0;
	s1 = sin(pangl1);
	s2 = sin(pangl2);
	s3 = sin(pangl3);
	c1 = cos(pangl1);
	c2 = cos(pangl2);
	c3 = cos(pangl3);

	rotationMatrix.push_back(std::vector<double>{c1*c2, -1 * c1*s2*s3 - s1*c3, -1 * c1*s2*c3 + s1*s3});
	rotationMatrix.push_back(std::vector<double>{s1*c2, -1 *s1*s2*s3 + c1*c3, -1 * s1*s2*c3 -c1*s3});
	rotationMatrix.push_back(std::vector<double>{s2, c2*s3 - s1*c3,c2*c3});
}

void Element::set_gxsi_geta_gzet(int i,std::vector<double>& jacobi,double GRAVX, double GRAVY, double GRAVZ)
{

		GXSI[i] = jacobi[0] * GRAVX + jacobi[1] * GRAVY + jacobi[2] * GRAVZ;
		GETA[i] = jacobi[3] * GRAVX + jacobi[4] * GRAVY + jacobi[5] * GRAVZ;
		GZET[i] = jacobi[6] * GRAVX + jacobi[7] * GRAVY + jacobi[8] * GRAVZ;
	
}

Element::~Element()
{
}
