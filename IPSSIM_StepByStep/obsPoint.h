#ifndef OBS_POINT_H
#define OBS_POINT_H
#pragma once
class obsPoint
{
public:
	obsPoint(std::string _name,double _x,double _y, double _z,std::string _schedule,std::string _formt);
	double get_x() const{ return x; }
	double get_y() const{ return y; }
	double get_z() const{ return z; }
	int get_element() const{ return in_element; }
	void set_element(int val){ in_element = val; }
	void set_xsi(double val){ xsi = val; }
	void set_eta(double val){ eta = val; }
	void set_zet(double val){ zet = val; }
	std::string get_obs_sch()const { return schedule; }
	std::string get_name() const{ return name; }
	std::string get_format() const{ return format; }
	~obsPoint();
private:
	std::string name;
	std::string schedule;
	std::string format;
	int in_element;
	double x;
	double y;
	double z;
	double xsi;
	double eta;
	double zet;
};
#endif
