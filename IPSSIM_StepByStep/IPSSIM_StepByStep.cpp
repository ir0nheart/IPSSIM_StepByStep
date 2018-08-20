// IPSSIM_StepByStep.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SimulationControl.h"


void main(int argc, _TCHAR* argv[])
{
	std::thread simControl = std::thread(&SimulationControl::run, SimulationControl::instance());
	simControl.join();
	std::cout << "Simulation Ended." << std::endl;
}

