#include "stdafx.h"
#include "obsPoint.h"


obsPoint::obsPoint(std::string _name, double _x, double _y, double _z, std::string _schedule, std::string _formt) :
name(_name), x(_x), y(_y), z(_z), schedule(_schedule), format(_formt)
{
}


obsPoint::~obsPoint()
{
}
