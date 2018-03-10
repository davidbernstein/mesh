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
 
#include <iostream>
#include <fstream>
#include "element.h"
#include "edge.h"
#include "meshconstants.h"

using namespace std;
using namespace NAMESPACE;

// copy constructor
Element::Element(const Element &e) : DataPacket(e)
{
	mID = e.mID;
	mMaterialNumber = e.mMaterialNumber;
	mpVertex = e.mpVertex;
	mpNeighbor = e.mpNeighbor;
	mpListIterator = e.mpListIterator;
	
	return;
}

	
	
Element::~Element()
{
    //DeleteListIterator();
    return;
}



// operators
Element& Element::operator=(const Element &e)
{
	// check for assignment to self
	if (this == &e)
		return *this;

	// copy the Point part of v
	DataPacket::operator=(e);
	
	mID = e.mID;
	mMaterialNumber = e.mMaterialNumber;
	mpVertex = e.mpVertex;
	mpNeighbor = e.mpNeighbor;
	
	// don't copy list iterator
	
	return *this;
}



void Element::SetNeighbors()
{
	AllocateNeighborArray();
		
	Array<Vertex*> face;
	Array<Element*> common;
	for (short i = 0; i < NumNeighbors(); ++i) {
        GetFace(i, face);
		if (face.Size() > 2)
			face[0]->CommonNeighbors(common, face[1], face[2], Type());
		else
			face[0]->CommonNeighbors(common, face[1], NULL, Type());
		
		switch (common.Size()) {
		case 1:
			if (common[0] != this)
				ThrowException("Element::SetNeighbors : bad vertex neighbor data");
			break;
			
		case 2:
			if (common[0] == this)
				SetNeighbor(i, common[1]);
			else
				SetNeighbor(i, common[0]);
			break;
		
		default:
			ThrowException("Element::SetNeighbors : bad vertex neighbor data");
			break;
		}
	}
	
	
	return;
}



void Element::Ring(Array<Element*> &ring, short edgeIndex) const
{	
	Vertex *pV[2];
	GetEdge(edgeIndex, pV);
	
	set<Element*> ringSet, waiting;
	ringSet.insert((Element*) this);
	waiting.insert((Element*) this);

	while (!waiting.empty()) {
		Element *pETop = *waiting.begin();
		waiting.erase(waiting.begin());
		
		for (short i = 0; i < pETop->NumNeighbors(); ++i) {
			if (pETop->mpNeighbor[i]->Type() != BLOB) {
				if (pETop->mpNeighbor[i]->EdgeIndex(pV[0], pV[1]) != -1) {
					if (ringSet.insert(pETop->mpNeighbor[i]).second) {
						waiting.insert(pETop->mpNeighbor[i]);
					}
				}
			}
		}
	}
	
	ring = ringSet;
	
	
	return;
}



void Element::Neighborhood(set<Element*> &neighborhood, short neighborhoodSize, ElementType eType) const
{
	neighborhood.clear();
	
	if (neighborhoodSize < 0) 
		return;
	
	neighborhood.insert((Element*) this);
	
	if (neighborhoodSize == 0) 
		return;
		
	set<Vertex*> vertexAll, vertexRing;
	
	for (short i = 0; i < mpVertex.Size(); ++i) {
		vertexAll.insert(mpVertex[i]);
		vertexRing.insert(mpVertex[i]);
	}
    
	set<Element*>::iterator iNs;
	set<Vertex*>::iterator iVs;
	for (short j = 1; j <= neighborhoodSize; ++j) {
		for (iVs = vertexRing.begin(); iVs != vertexRing.end(); ++iVs) 
			(*iVs)->AddNeighborsToSet(neighborhood, eType);
	
		vertexRing.clear();
		for (iNs = neighborhood.begin(); iNs != neighborhood.end(); ++iNs) {
			for (short i = 0; i < (*iNs)->NumVertices(); ++i) {
				if (vertexAll.find((*iNs)->PV(i)) == vertexAll.end())
					vertexRing.insert((*iNs)->PV(i));
			}
		}
		
		for (iVs = vertexRing.begin(); iVs != vertexRing.end(); ++iVs) 
			vertexAll.insert(*iVs);
	}
	
	
	return;
}



bool Element::Conforming() const
{
	for (short i = 0; i < NumNeighbors(); ++i) {
        if (PN(i)->Type() != BLOB) {
            if (PN(i)->GetBisectionData()->IsALeaf() == false)
                return false;
        }
	}
	
	return true;
}



void Element::Centroid(Point &centroid) const
{
	centroid.SetPosition(0.0, 0.0, 0.0);
	
	if (Type() != BLOB) {
		for (short i = 0; i < mpVertex.Size(); ++i) 
			centroid += *(mpVertex[i]);

		centroid /= (double) mpVertex.Size();
	}

	return;
}



short Element::LongestEdge() const
{	
	// longest edge with well-defined tie-breaking rule
	Vertex *pV[2];
	Vertex *pVCurrentMax[2];
	
	double maxDist = -1.0;
	short index = -1;
	for (short i = 0; i < NumEdges(); ++i) {
		GetEdge(i, pV);
		double d2 = pV[0]->SqrDistanceToPoint(pV[1]);
		
		if (FloatComparison(d2, maxDist, FLOAT_EQUALITY, LONGEST_EDGE_TOL)) {
			// this is the tie breaker
			pair<Vertex*, Vertex*> p0, pCM;
			p0 = (pV[0] > pV[1]) ? make_pair(pV[0], pV[1]) : make_pair(pV[1], pV[0]);
			pCM = (pVCurrentMax[0] > pVCurrentMax[1]) ? make_pair(pVCurrentMax[0], pVCurrentMax[1]) : make_pair(pVCurrentMax[1], pVCurrentMax[0]);
			if (p0 > pCM) {
				index = i;
				maxDist = d2;
				pVCurrentMax[0] = pV[0];
				pVCurrentMax[1] = pV[1];
			}
		}
		
		if (FloatComparison(d2, maxDist, FLOAT_GREATER_THAN, LONGEST_EDGE_TOL)) {
			index = i;
			maxDist = d2;
			pVCurrentMax[0] = pV[0];
			pVCurrentMax[1] = pV[1];
		}
	}

	return index;
}
    


void Element::CorrectVertices(const Array<Vertex*> &pV)
{
	for (short i = 0; i < mpVertex.Size(); ++i) 	
		mpVertex[i] = pV[mpVertex[i]->ID()];
	
	return;
}



void Element::AddVerticesToSet(std::set<Vertex*> &vertexSet) const
{
	for (short i = 0; i < mpVertex.Size(); ++i)
		vertexSet.insert(mpVertex[i]);
	
	return;
}



void Element::ReplaceVertex(const Vertex *pVOld, const Vertex *pVNew)
{	
	for (short i = 0; i < mpVertex.Size(); ++i) {
		if (mpVertex[i] == pVOld) {
			mpVertex[i] = (Vertex*) pVNew;
			return;
		}
	}
	
	return;
}



void Element::UpdateVertexNeighbors()
{
	for (short i = 0; i < mpVertex.Size(); ++i) 
		mpVertex[i]->AddNeighbor(this);
	
	return;
}



void Element::RemoveFromVertexNeighborSets()
{
	for (short i = 0; i < mpVertex.Size(); ++i) {
		if (mpVertex[i] != NULL) 
			mpVertex[i]->RemoveNeighbor(this);
	}
	
	return;
}



void Element::Print(string prefix) const
{
    cout << prefix << " " << this << endl;
    PrintVertices();
    //PrintNeighbors();
    
    return;
}



void Element::PrintVertices() const
{
	for (short j = 0; j < mpVertex.Size(); ++j) {
		cout << mpVertex[j] << " ";
        cout << "(" << mpVertex[j]->X() << ", " << mpVertex[j]->Y() << ", " << mpVertex[j]->Z() << ") ";
    }
    
	cout << endl;
	
	return;
}



void Element::PrintNeighbors() const
{
	cout << "this = " << this << endl;
	cout << "my vertices ";
	PrintVertices();
	for (short j = 0; j < mpNeighbor.Size(); ++j) {
		cout << mpNeighbor[j];
        if (mpNeighbor[j] != NULL)
            mpNeighbor[j]->PrintVertices();
    }
    
	cout << endl;
	
	return;
}



void Element::WriteVertexIDsToFile(ofstream &file) const
{
	for (short j = 0; j < mpVertex.Size(); ++j)
		file << mpVertex[j]->ID() << " ";
			
	file << endl;
	
	return;
}



void Element::CheckVertices() const
{
	for (short i = 0; i < mpVertex.Size(); ++i) {
		if (mpVertex[i] == NULL)
			ThrowException("Element::CheckVertices : NULL vertex");
		
		if (mpVertex[i]->NumNeighbors() > 0) {	
			if (mpVertex[i]->IsANeighbor(this) == false) 
				ThrowException("Element::CheckVertices : not a vertex neighbor");
		}
	}
	
	return;
}



void Element::CheckNeighbors() const
{
    Array<Vertex*> pVFace;

	for (short i = 0; i < mpNeighbor.Size(); ++i) {
		if (mpNeighbor[i] == NULL)
			ThrowException("Element::CheckNeighbors : NULL neighbor");
		
		if (mpNeighbor[i]->Type () != BLOB) {
			if (mpNeighbor[i]->IsANeighbor(this) == false)
				ThrowException("Element::CheckNeighbors : neighbor not a neighbor");
            
            GetFace(i, pVFace);
            for (short j = 0; j < pVFace.Size(); ++j) {
                if (mpNeighbor[i]->IsAVertex(pVFace[j]) == false)  {
                	cout << i << " " << j << endl;
                	PrintVertices();
                	mpNeighbor[i]->PrintVertices();
                	for (short k = 0; k < pVFace.Size(); ++k) {
                		cout << k << " " << mpNeighbor[i]->WhichVertex(pVFace[k]) << " " << pVFace[k] << endl;
                	}
                    ThrowException("Element::CheckNeighbors : neighbor not in correct position");
            	}
			}
		}
	}
	
	return;
}

