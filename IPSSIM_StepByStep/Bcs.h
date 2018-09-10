#ifndef BCS_H
#define BCS_H
#pragma once
class Bcs
{
public:
	Bcs();
	virtual ~Bcs();
	void setTimeStep(int ts);
	void setScheduleName(std::string scheduleName);
	void setBCSID(std::string BCSID);
	void setNumberOfQINC(int val);
	void setNumberOfUINC(int val);
	void setNumberOfPBC(int val);
	void setNumberOfUBC(int val);
	void addNode(int val);
	void addQINC(double val);
	void addQUINC(double val);
	void addUINC(double val);
	void addPBC(double val);
	void addUBC(double val);
	void addIsQINC(bool val);
	void addIsQUINC(bool val);
	void addIsUINC(bool val);
	void addIsPBC(bool val);
	void addIsUBC(bool val);

	int getTimeStep();
	std::string getScheduleName();
	std::string getBCSID();
	int getNumberOfQINC();
	int getNumberOfQUINC();
	int getNumberOfPBC();
	int getNumberOfUBC();
	std::vector<int>getNodes();
	std::vector<double>getQINC();
	std::vector<double>getQUINC();
	std::vector<double>getUINC();
	std::vector<double>getPBC();
	std::vector<double>getUBC();
	std::vector<bool>getIsQINC();
	std::vector<bool>getIsQUINC();
	std::vector<bool>getIsUINC();
	std::vector<bool>getIsPBC();
	std::vector<bool>getIsUBC();
private:
	int timeStep;
	std::string scheduleName;
	std::string BCSID;
	int numberOfQINC;
	int numberOfQUINC;
	int numberOfPBC;
	int numberOfUBC;

	std::vector<int> node;
	std::vector<double> QINC;
	std::vector<double> QUINC;
	std::vector<double> UINC;
	std::vector<double> PBC;
	std::vector<double> UBC;
	std::vector<bool> isQINC;
	std::vector<bool> isQUINC;
	std::vector<bool> isUINC;
	std::vector<bool> isPBC;
	std::vector<bool> isUBC;
};
#endif
