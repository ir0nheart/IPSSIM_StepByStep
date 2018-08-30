#include "stdafx.h"
#include "Schedule.h"


Schedule::Schedule(std::string name):schedule_name(name)
{
}

Schedule::Schedule(std::string name, std::string type):schedule_name(name),schedule_type(type)
{
}


Schedule::~Schedule()
{
}
