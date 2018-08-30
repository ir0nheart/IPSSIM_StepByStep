#include "stdafx.h"
#include "Writer.h"
#include <mutex>

std::vector<Writer *> Writer::instances = std::vector<Writer *>();

Writer * Writer::findWriter(std::string name)
{
	Writer * ptr = nullptr;
	for (Writer * wr : Writer::instances)
	{
		if (wr->name == name){
			ptr = wr;
			break;
		}
	}
	return ptr;
}
Writer * Writer::instance(std::string name)
{
	Writer * ptr = nullptr;
	if (instances.empty())
	{
		Writer * ins = new Writer(name);
		instances.push_back(ins);
		ptr = ins;
	} else
	{
		ptr = findWriter(name);
		if (ptr == nullptr)
		{
			Writer * ins = new Writer(name);
			instances.push_back(ins);
			ptr = ins;
		}
	}

	return ptr;
}

void Writer::write_to_file(std::string str)
{
	std::ofstream outfile;
	if (!newRun)
	{
		outfile.open(filename, std::ios::app);
		outfile << str << std::endl;
	} else
	{
		outfile.open(filename, std::ios::out);
		outfile << str << std::endl;
		newRun = false;
	}

}

void Writer::run()
{	
	std::string str = "";
	checkFILE:
	if (filename.empty()){
		std::this_thread::sleep_for(std::chrono::seconds(1));
		goto checkFILE;
	}

	checkContainer:
	if (!writeContainer.empty())
	{
		str = writeContainer[0];
		std::lock_guard<std::mutex> guard(mtx);
		writeContainer.pop_front();
		write_to_file(str);
		goto checkContainer;
	} else
	{
		std::this_thread::sleep_for(std::chrono::seconds(1));
		goto checkContainer;
	}

	std::cout << " I am Writer " << name << " and running.." << std::endl;
}
Writer::Writer(std::string name)
{
	this->name = name;
}


Writer::~Writer()
{
}
