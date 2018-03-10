/*
 * Copyright (C) 2004-2018 David Bernstein <david.h.bernstein@gmail.com>
 *
 * This file is part of Mesh.
 *
 * Mesh is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mesh is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mesh.  If not, see <http://www.gnu.org/licenses/>.
*/
 
#include "edge.h"
#include "triangle.h"

using namespace NAMESPACE;
using namespace std;

void Edge::Bisect(Element *pChild[2], Vertex *pVMid)
{
	BisectionData *pBD = GetBisectionData();
	BisectionData *pBD0 = pChild[0]->GetBisectionData();
	BisectionData *pBD1 = pChild[1]->GetBisectionData();
	
	pBD0->SetParent(this);
	pBD1->SetParent(this);
	pBD->SetChildren(pChild[0], pChild[1]);
			
	pChild[0]->SetVertices(PV(0), pVMid);
	pChild[1]->SetVertices(pVMid, PV(1));
	
	// set neighbors
	pChild[0]->SetNeighbor(0, pChild[1]);
	pChild[1]->SetNeighbor(1, pChild[0]);
    
	pChild[0]->SetNeighbor(1, PN(1));
	if (PN(1)->Type() != BLOB) 
        PN(1)->ReplaceNeighbor(this, pChild[0]);
	
    pChild[1]->SetNeighbor(0, PN(0));
	if (PN(0)->Type() != BLOB) 
        PN(0)->ReplaceNeighbor(this, pChild[1]);
    
	
	return;
}



Point Edge::MidPoint() const
{
	return (*PV(0) + *PV(1)) * 0.5;
}



bool Edge::ContainsPoint(const Point &p) const
{
    // this function assumes that this edge is located on the X axis
    // the Y and Z coordinates of p are ignored
    
   	return  (p.X() > PV(0)->X()) && (p.X() < PV(1)->X());
}



void Edge::OutwardAreaNormal(Array<Point> &unitNormal) const
{
	unitNormal.SetSize(2);
	
	unitNormal[0] = *PV(1) - *PV(0);
	unitNormal[1] = *PV(0) - *PV(1);
	
	unitNormal[0].Normalize();
	unitNormal[1].Normalize();
	
	return;
}



bool Edge::IntersectsEdge(Vertex *pV[2], Point &p) const
{
	// this function currently only works for straight edges
	Point d0 = *mpVertex[1] - *mpVertex[0];
	Point d1 = *(pV[0]) - *(pV[1]);
	Point d2 = *(pV[0]) - *mpVertex[0];
	
	double m[2][2], rhs[2], u[2];
	m[0][0] = d0.X();
	m[0][1] = d1.X();
	m[1][0] = d0.Y();
	m[1][1] = d1.Y();
	
	rhs[0] = d2.X();
	rhs[1] = d2.Y();
	
	double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
	if (det == 0.0)
		return false;
	
	double fac = 1.0 / det;
	u[0] = fac * (m[1][1] * rhs[0] - m[0][1] * rhs[1]);
	if ((u[0] >= 0.0) && (u[0] <= 1.0)) {
		u[1] = fac * (m[0][0] * rhs[1] - m[1][0] * rhs[0]);
		if ((u[1] >= 0.0) && (u[1] <= 1.0)) {
			p = d0 * u[0] + *mpVertex[0];
			return true;
		}
	}
	
	return false;
}



bool Edge::operator==(const Edge &e) const
{
	ThrowException("Edge::operator== called");
	return false;
}



bool Edge::operator<(const Edge &e) const
{
	ThrowException("Edge::operator< called");
	return false;
}



bool Edge::operator!=(const Edge &e) const
{
	ThrowException("Edge::operator!= called");
	return false;
}



bool Edge::operator>(const Edge &e) const
{
	ThrowException("Edge::operator> called");
	return false;
}
