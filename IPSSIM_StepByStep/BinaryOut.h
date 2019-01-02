#pragma once
#include <fstream>
#include <cstring>

class BinaryOut : public std::ofstream
{
public:
	BinaryOut(const char* fname) :
		std::ofstream(fname, std::ios::binary)
	{}

};

