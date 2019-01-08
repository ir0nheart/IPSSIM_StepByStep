#ifndef WRITER_H
#define WRITER_H
#pragma once
#include <mutex>

class Writer
{
public:
	static Writer * instance(std::string name);
	virtual void run();
	void set_filename(std::string fname){ filename = fname; }
	void add_line(std::string str){ writeContainer.push_back(str); }
	void set_done(){ done_writing = true; }
	void set_obs(bool b){ isObs = b; }
	void set_newRun(bool b){ newRun = b; }
private:
	std::mutex mtx;
	void write_to_file(std::string str);
	void write_to_file_a(std::string str);
	static std::vector<Writer *> instances;
	static Writer * findWriter(std::string name);
	std::string filename;
	std::deque<std::string> writeContainer;
	std::string name;
	bool done_writing;
	bool newRun;
	bool isObs;
	Writer(std::string name);
	~Writer();
};
#endif

