#ifndef PARSER_H
#define PARSER_H
#pragma once

class SimulationControl;
class InputFiles;
class Storage;

class Parser
{
public:
	static Parser* instance();
	static Parser* instance(std::string inputDirectory, Storage * storage);
	void mapFile(std::string fileName);
	void unmapFile(){ UnmapViewOfFile(mapViewOfFile); };
	void findDataSetPositionsInMap();
	void extractDataSets();
	void RemoveCharFromString(char * p, char c);
	int convertInt(std::vector<char> cInt);
	double convertDouble(std::vector<char> cDouble);
	void parseDataSet_1();
	void parseDataSet_2A();
	void parseDataSet_2B();
	void parseDataSet_3();
	void parseDataSet_4();
	void parseDataSet_5();
	void parseDataSet_6();
	void parseDataSet_7A();
	void parseDataSet_7B();
	void parseDataSet_7C();
	void parseDataSet_8ABC();
	void parseDataSet_8D();
	void parseDataSet_8E_9_10_11();
	void parseDataSet_12_13_14A();
	void parseDataSet_14B();
	void parseDataSet_15A();
	void parseDataSet_15B();
	void parseDataSet_17_18();
	void parseDataSet_19_20();
	void parseDataSet_22();
	void extractICS();
private:
	static Parser* m_PInstance;
	Storage * storage;
	std::string input_directory;
	char* mapViewOfFile;
	Parser(std::string inputDirectory,Storage * storage);
	std::vector<int>DataSetStart;
	std::vector<int>DataSetEnd;
	std::vector<std::pair<int, int>>dataSets;
	~Parser();
};
#endif

