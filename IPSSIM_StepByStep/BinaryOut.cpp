#include "stdafx.h"
#include "BinaryOut.h"

// use a generic func to use for everything which not handled in special way
template < typename DataType >
BinaryOut& operator << (BinaryOut& out, const DataType& data)
{
	out.write((char*)&data, sizeof(DataType));
	return out;
}

// if pointer type, write not pointer but values which pointer points to:
template < typename DataType >
BinaryOut& operator << (BinaryOut& out, const DataType*& data)
{
	out.write((char*)data, sizeof(DataType));
	return out;
}

// special case for char ptr ( old style c string ) which ends with '\0'
// use old style c here ( no std::string is involved )
BinaryOut& operator << (BinaryOut& out, const char* ptr)
{
	out.write(ptr, strlen(ptr));
	return out;
}


