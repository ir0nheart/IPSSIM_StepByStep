#include "stdafx.h"
#include "Node.h"
#include <algorithm>


Node::Node(int & num, int &reg, double &xx, double &yy, double &zz, double &ppor):node_num(num), nreg(reg), x(xx), y(yy), z(zz),por(ppor)
{
}


void Node::make_neighbors_unique()
{
	sort(neighbors.begin(), neighbors.end());
	neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

Node::~Node()
{
}
