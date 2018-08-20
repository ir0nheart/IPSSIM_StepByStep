#ifndef _TIMER_H
#define _TIMER_H
#include <ctime>
#pragma once
class Timer
{
public:
	Timer();
	void reset();
	operator double() const;
	~Timer();
private:
	std::clock_t m_start;
};
#endif