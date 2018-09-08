#ifndef SCHEDULE_H
#define SCHEDULE_H
#pragma once
class Schedule
{
public:
	Schedule(std::string name);
	Schedule(std::string name, std::string type);
	bool get_elapsed() const{ return elapsed; }
	bool get_sbased() const{ return sbased; }
	void set_step_at(int ind, int step){ step_time[ind].first = step; }
	void set_time_at(int ind, double time){ step_time[ind].second = time; }
	std::vector<std::pair<int, double>> get_step_time() const{ return step_time; }
	void add_step_time(int step, double time)
	{
		step_time.push_back(std::pair<int, double>(step, time));
	}
	void set_elapsed(bool val){ elapsed = val; }
	void set_sbased(bool val){ sbased = val; }
	std::string get_type() const{ return schedule_type; }
	std::string get_name() const{ return schedule_name; }
	int get_step_list_size()const { return step_time.size(); }
	double get_max_time(){ return step_time[step_time.size() - 1].second; }
private:
	std::string schedule_name;
	std::string schedule_type;
	std::vector<std::pair<int, double>> step_time;
	bool sbased;
	bool elapsed;
	~Schedule();
};
#endif