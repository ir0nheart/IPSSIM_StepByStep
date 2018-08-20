#include "stdafx.h"
#include "Timer.h"


Timer::Timer() :m_start(std::clock()) {
}

Timer::operator double() const {

	return (std::clock() - m_start) / static_cast<double>(CLOCKS_PER_SEC);
}
void Timer::reset(){
	m_start = std::clock();
}
Timer::~Timer()
{
}
