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
		outfile.close();
	} else
	{
		outfile.open(filename, std::ios::out);
		outfile << str << std::endl;
		outfile.close();
		newRun = false;
	}

}
void Writer::write_to_file_a(std::string str)
{
	std::ofstream outfile;
	outfile.open(filename, std::ios::app);
	outfile << str << std::endl;
	


}

void Writer::run()
{	
	std::string str = "";
	checkFILE:
	if (filename.empty()){
		std::this_thread::sleep_for(std::chrono::seconds(1));
		goto checkFILE;
	}

	if (!done_writing){
	checkContainer:

		
			if (!writeContainer.empty())
			{
				str = writeContainer[0];
				std::lock_guard<std::mutex> guard(mtx);
				writeContainer.pop_front();
				write_to_file(str);
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
				goto checkContainer;
			}
			
			std::this_thread::sleep_for(std::chrono::seconds(1));
			
		goto checkContainer;
		
	}
}
Writer::Writer(std::string name)
{
	this->name = name;
	done_writing = false;
	newRun = true;
	isObs = false;
}


Writer::~Writer()
{
}
