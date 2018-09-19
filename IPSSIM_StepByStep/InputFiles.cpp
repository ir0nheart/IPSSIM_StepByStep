#include "stdafx.h"
#include "InputFiles.h"
#include "simulationControl.h"
#include "Storage.h"
#include "Writer.h"

InputFiles* InputFiles::m_pInstance = nullptr;

InputFiles* InputFiles::instance()
{
	if (m_pInstance == nullptr){
		std::cout << "instance() must be called with input directory first." << std::endl;
		SimulationControl::exitOnError();
	}
	return m_pInstance;
}
InputFiles* InputFiles::instance(std::string inputDirectory)
{
	if (m_pInstance == nullptr)
		m_pInstance = new InputFiles(inputDirectory);
	return m_pInstance;
}
std::unordered_map<std::string, std::string> InputFiles::getFilesForReading()
{
	return filesToRead;
}

std::unordered_map<std::string, std::string> InputFiles::getFilesForWriting()
{
	return filesToWrite;
}

InputFiles::InputFiles(std::string inputDirectory) :input_directory(inputDirectory)
{
	std::string fname;
	fname.append(inputDirectory);
	fname.append("SUTRA.FIL");
	std::ifstream simulationFileList;

	simulationFileList.open(fname.c_str());
	if (simulationFileList.fail())
	{
		std::cerr << "Error : " << strerror(errno) << "\n " << std::endl;
		std::cout << "PLEASE MAKE SURE SUTRA.FIL IS IN YOUR WORKING DIRECTORY" << std::endl;
		SimulationControl::exitOnError();
	}

	std::string mline;
	while (std::getline(simulationFileList, mline)){  // This argument gets a line from file and saves in mline
		// Example Line of SUTRA.FIL
		// FILETYPE MODULENO(FORTRAN REQUIRES) 'FILENAME'
		// INP    50 'GT_SUM_GT_LW.inp'

		//convert mline String to char vector and push terminating character
		std::vector<char> lineBuffer(mline.begin(), mline.end());
		lineBuffer.push_back('\0');

		std::string filetype, filename; // two strings for storing required information filetype and filename, in c++ we don't need module number
		filetype = std::string(strtok(lineBuffer.data(), " "));  // parse the line and save filetype 
		std::stoi(strtok(NULL, " ")); // eliminate module number from line
		filename = std::string(strtok(NULL, " ")); // readfilename and save it. note that filename given like 'filename'

		std::string delimiter = "'"; // for parsing filename we use a different delimiter
		size_t pos = 0;
		std::string actFileName;
		while ((pos = filename.find(delimiter)) != std::string::npos) {
			actFileName = filename.substr(0, pos);
			filename.erase(0, pos + delimiter.length());
		}

		inputFileMap[filetype] = actFileName; // Put Input file into map, relate actual file name to file extension
	}
	setFilesForReading();
	setFilesForWriting();
	printAllFileInformation();
}

void InputFiles::setFilesForReading()
{
	if (inputFileMap["BCS"] != "")
	{
		filesToRead["BCS"] = inputFileMap["BCS"];
		Storage::instance()->set_bcs_defined(true);
	}
	else
	{
		std::cout << "WARNING!!" << std::endl;
		std::cout << "BCS (Boundary Conditions) input file for simulation is not defined in SUTRA.FIL file." << std::endl;
		Storage::instance()->set_bcs_defined(false);
	}

	if (inputFileMap["ICS"] != "")
	{
		filesToRead["ICS"] = inputFileMap["ICS"];
	}
	else
	{
		std::cout << "ERROR!!" << std::endl;
		std::cout << "ICS (Initial Conditions) input file for simulation is not defined in SUTRA.FIL file." << std::endl;
		SimulationControl::exitOnError();
	}

	if (inputFileMap["INP"] != "")
	{
		filesToRead["INP"] = inputFileMap["INP"];
	}
	else
	{
		std::cout << "ERROR!!" << std::endl;
		std::cout << "INP input file for simulation is not defined in SUTRA.FIL file." << std::endl;
		SimulationControl::exitOnError();
	}
}

void InputFiles::setFilesForWriting()
{
	filesToWrite["LST"] = inputFileMap["LST"];
	filesToWrite["RST"] = inputFileMap["RST"];
	filesToWrite["NOD"] = inputFileMap["NOD"];
	filesToWrite["ELE"] = inputFileMap["ELE"];
	filesToWrite["OBS"] = inputFileMap["OBS"];
	filesToWrite["SMY"] = inputFileMap["SMY"];
}

void InputFiles::printAllFileInformation()
{
	SimulationControl::wConsole("List of Files for Reading..", 14);
	
	std::unordered_map<std::string, std::string>::iterator it = filesToRead.begin();
	int i = 1;
	while (it != filesToRead.end())
	{
		char buff[1024];
		_snprintf(buff, 1024, "\t %d", i);
		SimulationControl::wConsolex(buff,13);
		SimulationControl::wConsolex("\t -- \t", 15);
		SimulationControl::wConsole(it->second.c_str(), 11);
		//std::cout << "\t" << i << " -- " << it->second << std::endl;
		it++;
		i++;
	}

	SimulationControl::wConsoleSpacer();
	SimulationControl::wConsole("List of Files for Writing..", 14);

	i = 1;
	it = filesToWrite.begin();
	while (it != filesToWrite.end())
	{
	//	std::cout << "\t" << i << " -- " << it->second << std::endl;
		char buff[1024];
		_snprintf(buff, 1024, "\t %d", i);
		SimulationControl::wConsolex(buff, 13);
		SimulationControl::wConsolex("\t -- \t", 15);
		SimulationControl::wConsole(it->second.c_str(), 11);
		it++;
		i++;
	}
	SimulationControl::wConsoleSpacer();
}


//Check if all required input files exists in working Directory
void InputFiles::checkInputFiles()
{
	SimulationControl::wConsole("Checking if input files exist...",14);
	std::unordered_map<std::string, std::string>::iterator it = filesToRead.begin();
	while (it != filesToRead.end()){
		std::string str;
		str.append(input_directory);
		str.append(it->second);
		std::ifstream fin(str.c_str());

		if (fin.fail())
		{
			
			SimulationControl::wConsole("Failed: ",12);
			SimulationControl::wConsolex("\t", 12);
			SimulationControl::wConsole(str.c_str(), 15);
			
			SimulationControl::wConsoleRight(" file is not in the current working directory.", 9);
			//std::cout << std::setw(80) << std::right << " file is not in the current working directory." << std::endl;
			SimulationControl::exitOnError();
		}

		SimulationControl::wConsole("Success :", 10);
		SimulationControl::wConsolex("\t", 12);
		SimulationControl::wConsole(str.c_str(), 15);
		SimulationControl::wConsoleRight(" file is in the current working directory.", 9);
		it++;
	}
	SimulationControl::wConsoleSpacer();
}

void InputFiles::printInputFilesToLST()
{
	Writer * lstWriter = Writer::instance("LST");
	std::string logLine;
	logLine.append("\n\n\n\n\n           F I L E   U N I T   A S S I G N M E N T S\n\n");
	logLine.append("             INPUT UNITS : ");
	logLine.append("\n              INP FILE (MAIN INPUT)                     ASSIGNED TO ");
	logLine.append(filesToRead["INP"]);
	logLine.append("\n              ICS FILE (INITIAL CONDITIONS)             ASSIGNED TO ");
	logLine.append(filesToRead["ICS"]);
	logLine.append("\n              BCS FILE (TIME-VAR. BND. COND.)           ASSIGNED TO ");
	logLine.append(filesToRead["BCS"]);
	logLine.append("\n\n             OUTPUT UNITS : ");
	logLine.append("\n              SMY FILE (RUN SUMMARY)                    ASSIGNED TO ");
	logLine.append(filesToWrite["SMY"]);
	logLine.append("\n              LST FILE (GENERAL OUTPUT)                 ASSIGNED TO ");
	logLine.append(filesToWrite["LST"]);
	logLine.append("\n              RST FILE (RESTART DATA)                   ASSIGNED TO ");
	logLine.append(filesToWrite["RST"]);
	logLine.append("\n              NOD FILE (NODEWISE OUTPUT)                ASSIGNED TO ");
	logLine.append(filesToWrite["NOD"]);
	logLine.append("\n              ELE FILE (VELOCITY OUTPUT)                ASSIGNED TO ");
	logLine.append(filesToWrite["ELE"]);
	logLine.append("\n              OBS FILE (OBSERVATION OUTPUT)             ASSIGNED TO ");
	logLine.append(filesToWrite["OBS"]);
	logLine.append("\n\n              NAMES FOR OBS AND OBC FILES WILL BE GENERATED AUTOMATICALLY FROM THE BASE NAMES LISTED ABOVE AND SCHEDULE NAMES");
	logLine.append("\n              LISTED LATER IN THIS FILE.UNIT NUMBERS ASSIGNED TO THESE FILES WILL BE THE FIRST AVAILABLE NUMBERS GREATER THAN");
	logLine.append("\n              OR EQUAL TO THE VALUES LISTED ABOVE IN PARENTHESES.");
	lstWriter->add_line(logLine);

}

// Read props.inp file for layer information and simulation time step divide
void InputFiles::readPropsINP()
{
	std::string mline;
	const char * del = " "; // space delimiter
	std::string propsFile;
	propsFile.append(input_directory);
	propsFile.append("props.inp");
	std::ifstream propsFileStream;


	propsFileStream.open(propsFile.c_str());
	if (propsFileStream.fail()){ 
		std::cerr << "Error: " << strerror(errno) << "\n" << std::endl;
		std::cout << "Please make sure 'props.inp' is in your work directory \n" << std::endl;
		SimulationControl::exitOnError();
	}
	filesToRead["PROPS"] = "props.inp";
	std::getline(propsFileStream, mline);
	std::getline(propsFileStream, mline);
	std::vector<char> lineBuffer(mline.begin(), mline.end());
	lineBuffer.push_back('\0');
	Storage::instance()->set_pstar(std::stod(strtok(lineBuffer.data(), " ")));
	Storage::instance()->set_gconst(std::stod(strtok(NULL, " ")));
	Storage::instance()->set_temp(std::stod(strtok(NULL, " ")));
	Storage::instance()->set_smwh(std::stod(strtok(NULL, " ")));
	Storage::instance()->set_time_step_divide(std::stod(strtok(NULL, " ")));
	Storage::instance()->set_water_table(std::stod(strtok(NULL, " ")));
	std::getline(propsFileStream, mline);
	std::getline(propsFileStream, mline);
	lineBuffer.assign(mline.begin(), mline.end());
	lineBuffer.push_back('\0');
	int noLayers = std::stoi(strtok(lineBuffer.data(), " "));
	Storage::instance()->set_number_of_layers(noLayers);
	std::getline(propsFileStream, mline);
	for (int i = 0; i < noLayers; i++)
	{
		std::getline(propsFileStream, mline);
		lineBuffer.assign(mline.begin(), mline.end());
		lineBuffer.push_back('\0');
		Storage::instance()->add_layerData(lineBuffer);
	}

}


// Destructor Object
InputFiles::~InputFiles()
{
}



