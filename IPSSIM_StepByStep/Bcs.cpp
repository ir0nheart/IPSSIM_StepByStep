#include "stdafx.h"
#include "Bcs.h"


Bcs::Bcs()
{
}


Bcs::~Bcs()
{
}

void Bcs::setTimeStep(int ts){ this->timeStep = ts; }
void Bcs::setScheduleName(std::string scheduleName){ this->scheduleName = scheduleName; }
void Bcs::setBCSID(std::string BcsID){ this->BCSID = BcsID; }
void Bcs::setNumberOfQINC(int val){ this->numberOfQINC = val; }
void Bcs::setNumberOfUINC(int val){ this->numberOfQUINC = val; }
void Bcs::setNumberOfPBC(int val){ this->numberOfPBC = val; }
void Bcs::setNumberOfUBC(int val){ this->numberOfUBC = val; }
void Bcs::addNode(int val){ this->node.push_back(val); }
void Bcs::addQINC(double val){ this->QINC.push_back(val); }
void Bcs::addQUINC(double val){ this->QUINC.push_back(val); }
void Bcs::addUINC(double val){ this->UINC.push_back(val); }
void Bcs::addPBC(double val){ this->PBC.push_back(val); }
void Bcs::addUBC(double val){ this->UBC.push_back(val); }
void Bcs::addIsQINC(bool val){ this->isQINC.push_back(val); }
void Bcs::addIsQUINC(bool val){ this->isQUINC.push_back(val); }
void Bcs::addIsUINC(bool val){ this->isUINC.push_back(val); }
void Bcs::addIsPBC(bool val){ this->isPBC.push_back(val); }
void Bcs::addIsUBC(bool val){ this->isUBC.push_back(val); }

int Bcs::getTimeStep(){ return this->timeStep; }
std::string Bcs::getScheduleName(){ return this->scheduleName; }
std::string Bcs::getBCSID(){ return this->BCSID; }
int Bcs::getNumberOfQINC(){ return this->numberOfQINC; }
int Bcs::getNumberOfQUINC(){ return this->numberOfQUINC; }
int Bcs::getNumberOfPBC(){ return this->numberOfPBC; }
int Bcs::getNumberOfUBC(){ return this->numberOfUBC; }
std::vector<int>Bcs::getNodes(){ return this->node; }
std::vector<double>Bcs::getQINC(){ return this->QINC; }
std::vector<double>Bcs::getQUINC(){ return this->QUINC; }
std::vector<double>Bcs::getUINC(){ return this->UINC; }
std::vector<double>Bcs::getPBC(){ return this->PBC; }
std::vector<double>Bcs::getUBC(){ return this->UBC; }

std::vector<bool>Bcs::getIsQINC(){ return this->isQINC; }
std::vector<bool>Bcs::getIsQUINC(){ return this->isQUINC; }
std::vector<bool>Bcs::getIsUINC(){ return this->isUINC; }
std::vector<bool>Bcs::getIsPBC(){ return this->isPBC; }
std::vector<bool>Bcs::getIsUBC(){ return this->isUBC; }