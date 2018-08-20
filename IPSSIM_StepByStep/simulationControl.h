#ifndef SIMULATION_CONTROL_H
#define SIMULATION_CONTROL_H
#pragma once
class SimulationControl
{
public:
	virtual void run();
	static SimulationControl * instance();
	static void wConsole(const char* s, WORD color);
	static void exitOnError();
private:
	static SimulationControl * m_pInstance;
	std::string inputDirectory;
	SimulationControl();
	~SimulationControl();
};
#endif

