#include "stdafx.h"
#include "Parser.h"
#include "SimulationControl.h"
#include "Timer.h"
#include "DataSet.h"
#include "Storage.h"
#include <cctype>

class Storage;

struct membuf : std::streambuf {
	membuf(char const* base, size_t size) {
		char* p(const_cast<char*>(base));
		this->setg(p, p, p + size);
	}
};
struct imemstream : virtual membuf, std::istream {
	imemstream(char const* base, size_t size)
		: membuf(base, size)
		, std::istream(static_cast<std::streambuf*>(this)) {
	}
};

Parser * Parser::m_PInstance = nullptr;


Parser* Parser::instance()
{
	if (m_PInstance == nullptr)
	{
		std::cout << "instance() must be called with input directory first." << std::endl;
		SimulationControl::exitOnError();
	}
	return m_PInstance;
}


Parser* Parser::instance(std::string inputDirectory,Storage * storage)
{
	if (m_PInstance == nullptr)
		m_PInstance = new Parser(inputDirectory,storage);
	return m_PInstance;
}

Parser::Parser(std::string inputDirectory,Storage * storage_):input_directory(inputDirectory),storage(storage_)
{
	DataSetStart.reserve(32);
	DataSetEnd.reserve(32);
	dataSets.reserve(25);
}

void Parser::mapFile(std::string fileName)
{
	std::string inpFile;
	inpFile.append(input_directory);
	inpFile.append(fileName);

	int ssize = inpFile.length();
	TCHAR *param = new TCHAR[ssize + 1];
	param[ssize] = 0;
	std::copy(inpFile.begin(), inpFile.end(), param);
	HANDLE hMapFile;
	HANDLE hFile;
	hFile = CreateFile(param,               // file to open
		GENERIC_READ | GENERIC_WRITE,          // open for reading
		FILE_SHARE_READ,       // share for reading
		NULL,                  // default security
		OPEN_EXISTING,         // existing file only
		FILE_ATTRIBUTE_NORMAL | FILE_FLAG_OVERLAPPED, // normal file
		NULL);                 // no attr. template


	if (hFile == INVALID_HANDLE_VALUE)
	{

		_tprintf(TEXT("\t Terminal failure: unable to open file \"%s\" for read.\n"), param);
		return;
	}
	hMapFile = CreateFileMapping(hFile, NULL, PAGE_READWRITE, 0, 0, NULL);
	if (hMapFile == NULL)
	{
		_tprintf(TEXT("\t Could not create file mapping object (%d).\n"),
			GetLastError());

	}
	mapViewOfFile = (char*)MapViewOfFile(hMapFile,   // handle to map object
		FILE_MAP_ALL_ACCESS, // read/write permission
		0,
		0,
		0);
	std::cout << "\t Successfully mapped " << fileName << " file." << std::endl;
}

void Parser::findDataSetPositionsInMap()
{
	int const size = strlen(mapViewOfFile);
	char* str_start;
	char* str_end;
	
	for (int i = 0; i <= size; ++i)
	{
		if (mapViewOfFile[i] == '#'){
			str_start = mapViewOfFile + i;
			str_end = mapViewOfFile + i + 6;
			std::string key(str_start, str_end);
			if ((key == "#Start") || (key == "# Data"))
				DataSetStart.push_back(i);
		}
	}

	for (int i : DataSetStart)
	{
		for (int j = 1; j < 512; j++)
		{
			if (mapViewOfFile[i + j] == '\n')
			{
				DataSetEnd.push_back(i + j);
				break;
			}
		}
			
	}

	if (DataSetStart.size() != DataSetEnd.size())
		std::cout << "Check Input File .. There are Errors " << std::endl;

	// pair data Set values
	std::vector<std::pair<int, int>> dataSetLocation;
	std::vector <std::pair<int, int>> dataSetOmit;
	for (int i = 0; i < DataSetStart.size(); i++){
		std::pair<int, int> p = { DataSetStart[i], DataSetEnd[i] };
		dataSetLocation.push_back(p);
	}
	dataSetOmit.push_back(dataSetLocation[2]);
	dataSetOmit.push_back(dataSetLocation[4]);
	dataSetOmit.push_back(dataSetLocation[12]);


	str_start = mapViewOfFile;
	int dslSize = dataSetLocation.size();
	std::pair<int, int> p = { 0, dataSetLocation[0].first };
	dataSets.push_back(p);
	bool bOmit;
	for (int i = 0; i < dslSize - 1; i++){
		bOmit = false;
		for (std::pair<int, int> omit : dataSetOmit){
			if ((dataSetLocation[i].second + 1 == omit.first) || (dataSetLocation[i].first == omit.first)){
				bOmit = true;
				break;
			}
		}
		if (!bOmit){
			if (dataSetLocation[i].second + 1 != dataSetLocation[i + 1].first){
				p = { dataSetLocation[i].second + 1, dataSetLocation[i + 1].first - 1 };
			}
			else{
				p = { dataSetLocation[i].second + 1, dataSetLocation[i + 1].first };
			}
			dataSets.push_back(p);
		}
	}

	p = { dataSetLocation[dslSize - 1].second + 1, size - 1 };
	dataSets.push_back(p);

}

void Parser::extractDataSets()
{
	
	parseDataSet_1();
	parseDataSet_2A();
	parseDataSet_2B();
	parseDataSet_3();
	parseDataSet_4();
	parseDataSet_5();
	parseDataSet_6();
	parseDataSet_7A();
	parseDataSet_7B();
	parseDataSet_7C();
	parseDataSet_8ABC();
	parseDataSet_8D();
	parseDataSet_8E_9_10_11();
	parseDataSet_12_13_14A();
	parseDataSet_14B();
	parseDataSet_15A();
	parseDataSet_15B();
	parseDataSet_17_18();
	parseDataSet_19_20();
	parseDataSet_22();
}


void Parser::parseDataSet_1()
{
	int size = dataSets[0].second - dataSets[0].first;
	char * str_start = mapViewOfFile;
	char * str_end;
	std::string line;
	std::vector<std::string> lines;
	for (int i = 0; i < size; i++)
	{
		if (mapViewOfFile[i] == '\r')
		{
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	for (std::string str : lines)
		storage->addTittle(str.c_str());
}

void Parser::parseDataSet_2A()
{
	int size = dataSets[1].second - dataSets[1].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[1].first;
	char * str_end;
	for (int i = dataSets[1].first; i < dataSets[1].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 2A should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	// Remove ' from string
	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin();i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);
				
		}
		*i = tmp;
	}

	std::vector<std::vector<char>> parsed_strings(5);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}
	// SUTRA VERSION 2.2 SOLUTE TRANSPORT

	storage->set_sutra_string(parsed_strings[0]);
	storage->set_version_string(parsed_strings[1]);
	storage->set_version_num_string(parsed_strings[2]);
	storage->set_simulation_type_string(parsed_strings[3]);
	storage->set_transport_string(parsed_strings[4]);
}

void Parser::parseDataSet_2B()
{
	int size = dataSets[2].second - dataSets[2].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[2].first;
	char * str_end;
	for (int i = dataSets[2].first; i < dataSets[2].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 2B should be one line." << std::endl;
		SimulationControl::exitOnError();
	}

	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
	}
	std::vector<std::vector<char>> parsed_strings(6);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}
	storage->set_mesh_dim_string(parsed_strings[0]);
	storage->set_mesh_type_string(parsed_strings[1]);
	if (parsed_strings.size() == 5)
	{
		storage->set_nn1_string(parsed_strings[3]);
		storage->set_nn2_string(parsed_strings[4]);
	} else if (parsed_strings.size() == 6)
	{
		storage->set_nn1_string(parsed_strings[3]);
		storage->set_nn2_string(parsed_strings[4]);
		storage->set_nn3_string(parsed_strings[5]);
	}
}

void Parser::parseDataSet_3()
{
	int size = dataSets[3].second - dataSets[3].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[3].first;
	char * str_end;
	for (int i = dataSets[3].first; i < dataSets[3].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 3 should be one line." << std::endl;
		SimulationControl::exitOnError();
	}

	std::vector<std::vector<char>> parsed_strings(7);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_nn(convertInt(parsed_strings[0]));
	storage->set_ne(convertInt(parsed_strings[1]));
	storage->set_npbc(convertInt(parsed_strings[2]));
	storage->set_nubc(convertInt(parsed_strings[3]));
	storage->set_nsop(convertInt(parsed_strings[4]));
	storage->set_nsou(convertInt(parsed_strings[5]));
	storage->set_nobs(convertInt(parsed_strings[6]));
}


void Parser::parseDataSet_4()
{
	int size = dataSets[4].second - dataSets[4].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[4].first;
	char * str_end;
	for (int i = dataSets[4].first; i < dataSets[4].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 4 should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
	}
	std::vector<std::vector<char>> parsed_strings(7);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_simulation_condition_string(parsed_strings[0]);
	storage->set_flow_type_string(parsed_strings[1]);
	storage->set_transport_type_string(parsed_strings[3]);
	storage->set_simulation_start_string(parsed_strings[5]);
	storage->set_solution_storage(convertInt(parsed_strings[6]));

}
void Parser::parseDataSet_5(){
	int size = dataSets[5].second - dataSets[5].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[5].first;
	char * str_end;
	for (int i = dataSets[5].first; i < dataSets[5].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 5 should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	std::vector<std::vector<char>> parsed_strings(3);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_up(stod(std::string(parsed_strings[0].begin(), parsed_strings[0].end())));
	storage->set_gnup(stod(std::string(parsed_strings[1].begin(), parsed_strings[1].end())));
	storage->set_gnuu(stod(std::string(parsed_strings[2].begin(), parsed_strings[2].end())));

}

void Parser::parseDataSet_6()
{
	int size = dataSets[6].second - dataSets[6].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[6].first;
	char * str_end;
	for (int i = dataSets[6].first; i < dataSets[6].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
		tmp.clear();
	}

	if (lines[lines.size() - 1].size() == 1 && lines[lines.size() - 1][0] != '-'){
		std::cout << "Data Set 6 must be ended with '-' after defined temporal conditions." << std::endl;
		SimulationControl::exitOnError();
	}


	std::vector<std::vector<char>> parsed_strings(3);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_nsch(convertInt(parsed_strings[0]));
	storage->set_npcyc(convertInt(parsed_strings[1]));
	storage->set_nucyc(convertInt(parsed_strings[2]));

	for (int i = 1; i < lines.size()-1; i++)
	{
		storage->add_temporal_control(lines[i]);
	}

}

void Parser::parseDataSet_7A()
{
	int size = dataSets[7].second - dataSets[7].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[7].first;
	char * str_end;
	for (int i = dataSets[7].first; i < dataSets[7].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 7A should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	std::vector<std::vector<char>> parsed_strings(3);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}
	storage->set_itrmax(convertInt(parsed_strings[0]));
	if (parsed_strings[1].size() >0)
	{
		storage->set_rpmax(convertInt(parsed_strings[1]));
		storage->set_rumax(convertInt(parsed_strings[2]));
	}

}

void Parser::parseDataSet_7B()
{
	int size = dataSets[8].second - dataSets[8].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[8].first;
	char * str_end;
	for (int i = dataSets[8].first; i < dataSets[8].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 7B should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
	}
	std::vector<std::vector<char>> parsed_strings(3);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_p_solver_string(parsed_strings[0]);
	if (parsed_strings[1].size() > 0)
	{
		storage->set_max_p_iterations(convertInt(parsed_strings[1]));
		storage->set_p_tolerance(stod(std::string(parsed_strings[2].begin(), parsed_strings[2].end())));
	}
}

void Parser::parseDataSet_7C()
{
	int size = dataSets[9].second - dataSets[9].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[9].first;
	char * str_end;
	for (int i = dataSets[9].first; i < dataSets[9].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 7C should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
	}
	std::vector<std::vector<char>> parsed_strings(3);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_u_solver_string(parsed_strings[0]);
	if (parsed_strings[1].size() > 0)
	{
		storage->set_max_u_iterations(convertInt(parsed_strings[1]));
		storage->set_u_tolerance(stod(std::string(parsed_strings[2].begin(), parsed_strings[2].end())));
	}
}

void Parser::parseDataSet_8ABC()
{
	int size = dataSets[10].second - dataSets[10].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[10].first;
	char * str_end;
	for (int i = dataSets[10].first; i < dataSets[10].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 3)
	{
		std::cout << "Data Set 8A,8B and 8C should be one line each." << std::endl;
		SimulationControl::exitOnError();
	}

	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
		tmp.clear();
	}

	std::vector<std::vector<char>> data_set_8a(13);
	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (data_set_8a[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_8a[it].push_back(c);
		if (c == ' ')
			it++;
	}


	std::vector<std::vector<char>> data_set_8b(13);
	it = 0;
	for (char c : lines[1])
	{
		if (it >0)
		{
			if (data_set_8b[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_8b[it].push_back(c);
		if (c == ' ')
			it++;
	}

	std::vector<std::vector<char>> data_set_8c(11);
	it = 0;
	for (char c : lines[2])
	{
		if (it >0)
		{
			if (data_set_8c[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_8c[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_simulation_output_every(convertInt(data_set_8a[0]));
	storage->set_node_output_every(convertInt(data_set_8b[0]));
	storage->set_element_output_every(convertInt(data_set_8c[0]));

	for (int i = 1; i < data_set_8a.size()-3; i++)
		storage->add_simulation_output_controls(data_set_8a[i]);

	for (int i = 1; i < data_set_8b.size()-3; i++)
		storage->add_node_output_headers(data_set_8b[i]);

	for (int i = 1; i < data_set_8c.size()-3; i++)
		storage->add_element_output_headers(data_set_8c[i]);



}

void Parser::parseDataSet_8D()
{
	int size = dataSets[11].second - dataSets[11].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[11].first;
	char * str_end;
	for (int i = dataSets[11].first; i < dataSets[11].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() != storage->get_nobs() + 2)
	{
		std::cout << "Error in Data Set 8D.. " << std::endl;
		SimulationControl::exitOnError();
	}

	if (lines[lines.size() - 1].size() != 1 || lines[lines.size() - 1][0] != '-')
	{
		std::cout << "Last line of Data Set 8D must be '-'.." << std::endl;
		SimulationControl::exitOnError();
	}
	storage->set_observation_first_line(convertInt(lines[0]));

	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
		tmp.clear();
	}
	for (int i = 1; i < lines.size() - 1; i++){
		std::vector<std::vector<char>> parsed_strings(13);
		int it = 0;
		for (char c : lines[i])
		{
			if (it > 0)
			{
				if (parsed_strings[it - 1].empty())
					it = it - 1;
			}
			if (c != ' ')
				parsed_strings[it].push_back(c);
			if (c == ' ')
				it++;
		}

		storage->add_obsData(parsed_strings);
	}
}

void Parser::parseDataSet_8E_9_10_11()
{
	int size = dataSets[12].second - dataSets[12].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[12].first;
	char * str_end;
	for (int i = dataSets[12].first; i < dataSets[12].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
		tmp.clear();
	}

	std::vector<std::vector<char>> data_set_8e(8);
	int it = 0;
	for (char c : lines[0])
	{
		if (it > 0)
		{
			if (data_set_8e[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_8e[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_nbcfpr(convertInt(data_set_8e[0]));
	storage->set_nbcspr(convertInt(data_set_8e[1]));
	storage->set_nbcppr(convertInt(data_set_8e[2]));
	storage->set_nbcupr(convertInt(data_set_8e[3]));
	storage->set_cinact(data_set_8e[4][0]);
	std::vector<std::vector<char>> data_set_9(10);
	it = 0;
	for (char c : lines[1])
	{
		if (it > 0)
		{
			if (data_set_9[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_9[it].push_back(c);
		if (c == ' ')
			it++;
	}
	storage->set_compfl(std::stod(std::string(data_set_9[0].begin(), data_set_9[0].end())));
	storage->set_cw(std::stod(std::string(data_set_9[1].begin(), data_set_9[1].end())));
	storage->set_sigmaw(std::stod(std::string(data_set_9[2].begin(), data_set_9[2].end())));
	storage->set_rhow0(std::stod(std::string(data_set_9[3].begin(), data_set_9[3].end())));
	storage->set_urhow0(std::stod(std::string(data_set_9[4].begin(), data_set_9[4].end())));
	storage->set_drwdu(std::stod(std::string(data_set_9[5].begin(), data_set_9[5].end())));
	storage->set_visc0(std::stod(std::string(data_set_9[6].begin(), data_set_9[6].end())));
	std::vector<std::vector<char>> data_set_10(7);

	it = 0;
	for (char c : lines[2])
	{
		if (it > 0)
		{
			if (data_set_10[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_10[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_compfl(stod(std::string(data_set_10[0].begin(), data_set_10[0].end())));
	storage->set_cs(stod(std::string(data_set_10[1].begin(), data_set_10[1].end())));
	storage->set_sigmas(stod(std::string(data_set_10[2].begin(), data_set_10[2].end())));
	storage->set_rhos(stod(std::string(data_set_10[3].begin(), data_set_10[3].end())));
	

	std::vector<std::vector<char>> data_set_11(8);

	it = 0;
	for (char c : lines[3])
	{
		if (it > 0)
		{
			if (data_set_11[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_11[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_adsorption_string(data_set_11[0]);
}

void Parser::parseDataSet_12_13_14A()
{
	int size = dataSets[13].second - dataSets[13].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[13].first;
	char * str_end;
	for (int i = dataSets[13].first; i < dataSets[13].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	std::vector<char> tmp;
	for (std::vector<std::vector<char> >::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		for (std::vector<char>::iterator j = i->begin(); j != i->end(); ++j)
		{
			if (*j != '\'')
				tmp.push_back(*j);

		}
		*i = tmp;
		tmp.clear();
	}


	std::vector<std::vector<char>> data_set_12(20);
	int it = 0;
	for (char c : lines[0])
	{
		if (it > 0)
		{
			if (data_set_12[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_12[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_prodf0(stod(std::string(data_set_12[0].begin(), data_set_12[0].end())));
	storage->set_prods0(stod(std::string(data_set_12[1].begin(), data_set_12[1].end())));
	storage->set_prodf1(stod(std::string(data_set_12[2].begin(), data_set_12[2].end())));
	storage->set_prods1(stod(std::string(data_set_12[3].begin(), data_set_12[3].end())));

	std::vector<std::vector<char>> data_set_13(20);
	it = 0;
	for (char c : lines[1])
	{
		if (it > 0)
		{
			if (data_set_13[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_13[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_gravx(stod(std::string(data_set_13[0].begin(), data_set_13[0].end())));
	storage->set_gravy(stod(std::string(data_set_13[1].begin(), data_set_13[1].end())));
	storage->set_gravz(stod(std::string(data_set_13[2].begin(), data_set_13[2].end())));
	std::vector<std::vector<char>> data_set_14A(20);
	it = 0;
	for (char c : lines[2])
	{
		if (it > 0)
		{
			if (data_set_14A[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			data_set_14A[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_scalx(stod(std::string(data_set_14A[1].begin(), data_set_14A[1].end())));
	storage->set_scaly(stod(std::string(data_set_14A[2].begin(), data_set_14A[2].end())));
	storage->set_scalz(stod(std::string(data_set_14A[3].begin(), data_set_14A[3].end())));
	storage->set_porfac(stod(std::string(data_set_14A[4].begin(), data_set_14A[4].end())));
}

void Parser::parseDataSet_14B()
{
	int size = dataSets[14].second - dataSets[14].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	//lines.reserve(storage->get_nn());
	storage->reserveNodeData();
	char * str_start = mapViewOfFile + dataSets[14].first;
	char * str_end;
	for (int i = dataSets[14].first; i < dataSets[14].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			//lines.push_back(line);
			storage->add_nodeData(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	storage->reserveNodeData();

	/*for (int i = 0; i < lines.size();i++){
		std::vector<std::vector<char>> data_set_14B(6);
		int it = 0;
		for (char c : lines[i])
		{
			if (it > 0)
			{
				if (data_set_14B[it - 1].empty())
					it = it - 1;
			}
			if (c != ' ')
				data_set_14B[it].push_back(c);
			if (c == ' ')
				it++;
		}
		storage->add_nodeData(data_set_14B);
	}*/
}


void Parser::parseDataSet_15A()
{
	
	int size = dataSets[15].second - dataSets[15].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	char * str_start = mapViewOfFile + dataSets[15].first;
	char * str_end;
	for (int i = dataSets[15].first; i < dataSets[15].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			lines.push_back(line);
			str_start = mapViewOfFile + i + 2;
		}
	}

	if (lines.size() > 1)
	{
		std::cout << "Data Set 15A should be one line." << std::endl;
		SimulationControl::exitOnError();
	}
	std::vector<std::vector<char>> parsed_strings(20);

	int it = 0;
	for (char c : lines[0])
	{
		if (it >0)
		{
			if (parsed_strings[it - 1].empty())
				it = it - 1;
		}
		if (c != ' ')
			parsed_strings[it].push_back(c);
		if (c == ' ')
			it++;
	}

	storage->set_element_props(parsed_strings);

}

void Parser::parseDataSet_15B()
{


	int size = dataSets[16].second - dataSets[16].first;
	std::vector<char> line;
	std::vector<std::vector<char>> lines;
	//lines.reserve(storage->get_ne());
	storage->reserveElementData();
	char * str_start = mapViewOfFile + dataSets[16].first;
	char * str_end;
	for (int i = dataSets[16].first; i < dataSets[16].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			//lines.push_back(line);
			storage->add_elementData(line);
			str_start = mapViewOfFile + i + 2;
		}
	}


	//storage->reserveElementData();

	//for (int i = 0; i < lines.size(); i++){
	//	std::vector<std::vector<char>> data_set_15B(14);
	//	int it = 0;
	//	for (char c : lines[i])
	//	{
	//		if (it > 0)
	//		{
	//			if (data_set_15B[it - 1].empty())
	//				it = it - 1;
	//		}
	//		if (c != ' ')
	//			data_set_15B[it].push_back(c);
	//		if (c == ' ')
	//			it++;
	//	}
	//	storage->add_elementData(data_set_15B);
	//}
}

void Parser::parseDataSet_17_18()
{
	if (storage->get_nsop() != 0)
	{
		int size = dataSets[17].second - dataSets[17].first;
		std::vector<char> line;
		std::vector<std::vector<char>> lines;
	//	lines.reserve(storage->get_nsop());
		storage->reserve_nsop_data();
		char * str_start = mapViewOfFile + dataSets[17].first;
		char * str_end;
		for (int i = dataSets[17].first; i < dataSets[17].second; i++){
			if (mapViewOfFile[i] == '\r'){
				str_end = mapViewOfFile + i;
				line.assign(str_start, str_end);
				//line.push_back('\0');
				//lines.push_back(line);
				storage->add_nsop_data(line);
				str_start = mapViewOfFile + i + 2;
			}
		}


	}
	if (storage->get_nsou() != 0)
	{
		int size = dataSets[18].second - dataSets[18].first;
		std::vector<char> line;
		std::vector<std::vector<char>> lines;
		//lines.reserve(storage->get_nsou());
		storage->reserve_nsou_data();
		char * str_start = mapViewOfFile + dataSets[18].first;
		char * str_end;
		for (int i = dataSets[18].first; i < dataSets[18].second; i++){
			if (mapViewOfFile[i] == '\r'){
				str_end = mapViewOfFile + i;
				line.assign(str_start, str_end);
				//line.push_back('\0');
				//lines.push_back(line);
				storage->add_nsou_data(line);
				str_start = mapViewOfFile + i + 2;
			}
		}
		
	}
	
}

void Parser::parseDataSet_19_20()
{
	if (storage->get_npbc() != 0)
	{
		int size = dataSets[19].second - dataSets[19].first;
		std::vector<char> line;
		std::vector<std::vector<char>> lines;
		//lines.reserve(storage->get_npbc());
		storage->reserve_npbc_data();
		char * str_start = mapViewOfFile + dataSets[19].first;
		char * str_end;
		for (int i = dataSets[19].first; i < dataSets[19].second; i++){
			if (mapViewOfFile[i] == '\r'){
				str_end = mapViewOfFile + i;
				line.assign(str_start, str_end);
				//line.push_back('\0');
				//lines.push_back(line);
				storage->add_npbc_data(line);
				str_start = mapViewOfFile + i + 2;
			}
		}

		
	}

	if (storage->get_nubc() != 0)
	{
		int size = dataSets[20].second - dataSets[20].first;
		std::vector<char> line;
		std::vector<std::vector<char>> lines;
		//lines.reserve(storage->get_nubc());
		storage->reserve_nubc_data();
		char * str_start = mapViewOfFile + dataSets[20].first;
		char * str_end;
		for (int i = dataSets[20].first; i < dataSets[20].second; i++){
			if (mapViewOfFile[i] == '\r'){
				str_end = mapViewOfFile + i;
				line.assign(str_start, str_end);
				//line.push_back('\0');
				//lines.push_back(line);
				storage->add_nubc_data(line);
				str_start = mapViewOfFile + i + 2;
			}
		}
		

	}
	
}
void Parser::parseDataSet_22()
{
	int size = dataSets[21].second - dataSets[21].first;
	std::vector<char> line;
	//std::vector<std::vector<char>> lines;
	//lines.reserve(storage->get_ne());
	storage->reserve_incidence_data();
	char * str_start = mapViewOfFile + dataSets[21].first;
	char * str_end;
	for (int i = dataSets[21].first; i < dataSets[21].second; i++){
		if (mapViewOfFile[i] == '\r'){
			str_end = mapViewOfFile + i;
			line.assign(str_start, str_end);
			//line.push_back('\0');
			//lines.push_back(line);
			storage->add_incidence_data(line);
			str_start = mapViewOfFile + i + 2;
		}
	}
}
void Parser::RemoveCharFromString(char * p, char c)
{
	if (p == NULL)
		return;
	char * pDest = p;
	while (*p)
	{
		if (*p != c)
			*pDest++ = *p;
		p++;
	}
	*pDest = '\n';
}


double Parser::convertDouble(std::vector<char> str)
{
	int cSize = str.size();
	int ival = 0;
	int ictr = 0;
	int dctr = 0;
	double dval = 0;
	double val = 0;
	int dotIndx = -1;
	int mult = 1;
	if (cSize > 0)
	{
		for (int i = 0; i < cSize; i++)
		{
			if (str[i] == '.'){
				dotIndx = i;
				break;
			}
		}

		if (dotIndx != -1){

			if (cSize > 1 && (str[0] != '-' && str[0] != '+')){
				ictr = dotIndx;
				ictr--;
				for (int j = 0; j < dotIndx; j++)
				{
					int b = str[j] - '0';
					ival = ival + b * pow(10, ictr);
					ictr--;
				}
				dctr = 1;
				for (int k = dotIndx + 1; k < cSize; k++)
				{
					int b = str[k] - '0';
					dval = dval + b / pow(10, dctr);
					dctr++;
				}
				val = ival + dval;
			}
			else
			{
				if (str[0] == '-')
					mult = -1;

				ictr = dotIndx-1;
				ictr--;
				for (int j = 1; j < dotIndx; j++)
				{
					int b = str[j] - '0';
					ival = ival + b * pow(10, ictr);
					ictr--;
				}
				dctr = 1;
				for (int k = dotIndx + 1; k < cSize; k++)
				{
					int b = str[k] - '0';
					dval = dval + b / pow(10, dctr);
					dctr++;
				}
				val = ival + dval;
				val = mult * val;
			}
		}
	}
	return val;
}

int Parser::convertInt(std::vector<char> cInt)
{
	int cSize = cInt.size();
	int val = 0;
	int mult = 1;
	if (cSize > 0)
	{
		if (cSize > 1 && (cInt[0] != '-' && cInt[0] != '+')){
			cSize--;
			for (char i : cInt)
			{
				int b = i - '0';
				val = val + b * pow(10, cSize);
				cSize--;
			}
		} else
		{
			if (cInt[0] == '-')
				mult = -1;
			cSize-=2;
			for (int i = 1; i < cInt.size();i++)
			{
				int b = cInt[i] - '0';
				val = val + b * pow(10, cSize);
				val = mult*val;
				cSize--;
			}
		}
	} 
	return val;		
}
Parser::~Parser()
{
}

