#include "stdafx.h"
#include "Element.h"

double Element::GXLOC[8] = { -1, 1, 1, -1, -1, 1, 1, -1 };
double Element::GYLOC[8] = { -1, -1, 1, 1, -1, -1, 1, 1 };
double Element::GZLOC[8] = { -1, -1, -1, -1, 1, 1, 1, 1 };
Element::Element()
{
}



Element::~Element()
{
}
