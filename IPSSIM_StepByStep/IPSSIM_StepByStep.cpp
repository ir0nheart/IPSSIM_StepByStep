// IPSSIM_StepByStep.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "SimulationControl.h"
#include "Writer.h"


void main(int argc, _TCHAR* argv[])
{


	std::thread simControl = std::thread(&SimulationControl::run, SimulationControl::instance());
	std::thread lstwriter = std::thread(&Writer::run, Writer::instance("LST"));
	std::thread nodwriter = std::thread(&Writer::run, Writer::instance("NOD"));
	std::thread smywriter = std::thread(&Writer::run, Writer::instance("SMY"));
	std::thread elewriter = std::thread(&Writer::run, Writer::instance("ELE"));
	std::thread obswriter = std::thread(&Writer::run, Writer::instance("OBS"));

	simControl.join();
	lstwriter.join();
	nodwriter.join();
	smywriter.join();
	elewriter.join();
	obswriter.join();
	std::cout << "Simulation Ended." << std::endl;
}

