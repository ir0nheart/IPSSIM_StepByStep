#include "stdafx.h"
#include "SimulationControl.h"
#include "Timer.h"
#include "InputFiles.h"
#include <string>
#include "Parser.h"
#include "Storage.h"
#include "Writer.h"


SimulationControl * SimulationControl::m_pInstance = nullptr;

enum Colors{
	BLACK = 0,
	NAVY = 1,
	DARK_GREEN = 2,
	DARK_CYAN = 3,
	DARK_RED = 4,
	DARK_PINK = 5,
	DARK_YELLOW = 6,
	LIGHT_GRAY = 7,
	DARK_GRAY = 8,
	BRIGHT_BLUE = 9,
	BRIGHT_GREEN = 10,
	BRIGHT_CYAN = 11,
	BRIGHT_RED = 12,
	BRIGHT_PINK = 13,
	BRIGHT_YELLOW = 14,
	BRIGHT_WHITE = 15
};

SimulationControl * SimulationControl::instance()
{
	if (m_pInstance == nullptr)
		m_pInstance = new SimulationControl;
	return m_pInstance;
}


SimulationControl::SimulationControl()
{
	std::cout << "Simulation Control is Created Now " << std::endl;
	char result[MAX_PATH];
	GetModuleFileNameA(0, result, MAX_PATH);
	std::string strPath(result);
	strPath = strPath.substr(0, strPath.find_last_of("\\"));
	strPath.append("\\");
	inputDirectory = strPath;
	std::cout << "Path is " << inputDirectory.c_str() << std::endl;
}


SimulationControl::~SimulationControl()
{
}

void SimulationControl::exitOnError()
{
	wConsole("Terminating Program in 10 seconds", BRIGHT_RED);
	std::string msg;
	for (int i = 10; i >0; i--){
		msg.append(std::to_string(i));
		msg.append(" seconds to terminate the program.");
		wConsole(msg.c_str(), BRIGHT_WHITE);
		msg.clear();
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
	wConsole("BYEE..", BRIGHT_RED);
	exit(1);
}
void SimulationControl::exitOnError(std::string errcod)
{
	std::string msg;
	wConsole(errcod.c_str(), BRIGHT_CYAN);
	wConsole("Terminating Program in 10 seconds", BRIGHT_RED);
	msg.clear();
	for (int i = 10; i >0; i--){
		msg.append(std::to_string(i));
		msg.append(" seconds to terminate the program.");
		wConsole(msg.c_str(), BRIGHT_WHITE);
		msg.clear();
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
	wConsole("BYEE..", BRIGHT_RED);
	exit(1);
}
void SimulationControl::wConsole(const char* s, WORD color)
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
	std::cout << s << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);

}
void SimulationControl::run()
{
//	double * p = new double[100]{0};

	InputFiles * iFiles = InputFiles::instance(inputDirectory);
	wConsole("I am running on another thread", BRIGHT_YELLOW);
	iFiles->checkInputFiles();
	iFiles->readPropsINP();
	iFiles->printAllFileInformation();

	Storage * store = Storage::instance();

	Parser * fileParser = Parser::instance(inputDirectory,store);
	fileParser->mapFile(iFiles->getFilesForReading()["INP"]);
	fileParser->findDataSetPositionsInMap();

	Timer t;
	fileParser->extractDataSets();
	std::cout << "All Data Sets are extracted in " << t << "seconds."<< std::endl;
	fileParser->unmapFile();
	fileParser->mapFile(iFiles->getFilesForReading()["ICS"]);
	fileParser->extractICS();
	fileParser->unmapFile();
	store->check_data_sets();

	// HERE SUTRA STARTS
	std::cout << "FAST Check " << std::endl;
	
}