#pragma once
class DataSet
{
public:
	void free();
	char * data;
	DataSet(std::string dataSetName);
	~DataSet();
};

