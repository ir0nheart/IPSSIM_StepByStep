// IPSSIM_StepByStep.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SimulationControl.h"
#include "Writer.h"


void main(int argc, _TCHAR* argv[])
{

	std::thread simControl = std::thread(&SimulationControl::run, SimulationControl::instance());
	std::thread writer = std::thread(&Writer::run, Writer::instance("LST"));
	std::thread writer2 = std::thread(&Writer::run, Writer::instance("NOD"));
	simControl.join();
	writer.join();
	writer2.join();
	std::cout << "Simulation Ended." << std::endl;
}

