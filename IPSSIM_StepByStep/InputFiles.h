#ifndef INPUT_FILES_H
#define INPUT_FILES_H
#include "stdafx.h"
#pragma once
class InputFiles
{
public:
	static InputFiles * instance();
	static InputFiles * instance(std::string inputDirectory);

	void setFilesForReading();
	void setFilesForWriting();
	void printAllFileInformation();

	void getFileList();
	std::string getInputDirectory();
	void checkInputFiles();
	void printInputFilesToLST();
	void readPropsINP();
	std::unordered_map<std::string, std::string> getFilesForReading();
private:
	InputFiles();
	InputFiles(std::string inputDirectory);
	~InputFiles();
	static InputFiles * m_pInstance;
	std::string input_directory;
	std::unordered_map<std::string, std::string> inputFileMap;
	std::unordered_map<std::string, std::string> filesToRead;
	std::unordered_map<std::string, std::string> filesToWrite;
};
#endif

