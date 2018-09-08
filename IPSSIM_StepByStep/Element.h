#ifndef ELEMENT_H
#define ELEMENT_H
#pragma once
#define GLOC 0.577350269189626
#define D2R 0.01745329252

class Element
{
public:
	Element();
	static double get_GXLOC(int val) { return GXLOC[val]; }
	static double get_GYLOC(int val) { return GYLOC[val]; }
	static double get_GZLOC(int val) { return GZLOC[val]; }
	int get_el_num() const{ return *el_num; }
	~Element();
private:
	static double GXLOC[8];
	static double GYLOC[8];
	static double GZLOC[8];
	double * almax;
	double * almin;
	double * atmax;
	double * atmin;
	double * almid;
	double * atmid;
	double vmag;
	double vang1;
	double vang2;
	double * angle1;
	double * angle2;
	double * angle3;
	double * permxx;
	double * permxy;
	double * permyx;
	double * permyy;
	double * permmax;
	double * permmid;
	double * permmin;
	double * permxz;
	double * permyz;
	double * permzz;
	double * permzx;
	double * permzy;
	double * pangl1;
	double * pangl2;
	double * pangl3;	
	int * lreg;
	int * el_num;
};
#endif
