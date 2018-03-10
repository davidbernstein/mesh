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
 
#include "vertex.h"
#include "element.h"
#include "mesh.h"

using namespace std;
using namespace NAMESPACE;


// copy constructor
Vertex::Vertex(const Vertex &v) : Point(v)
{
	mNeighborSet = v.mNeighborSet;
	
	// don't copy list iterator
	mpListIterator = NULL;

	return;
}



Vertex& Vertex::operator=(const Vertex &v)
{
	// check for assignment to self
	if (this == &v)
		return *this;

	// copy the Point part of v
	Point::operator=(v);

	mNeighborSet = v.mNeighborSet;
	
	mpListIterator = NULL;


	return *this;
}



Vertex& Vertex::operator=(const Point &p)
{
	// copy the Point part of v
	Point::operator=(p);

	EraseNeighbors();
	
	mpListIterator = NULL;

	return *this;
}



void Vertex::Copy(const Vertex &v)
{
	// check for assignment to self
	if (this == &v)
		return;

	// copy the Point part of v
	Point::operator=(v);

	mNeighborSet = v.mNeighborSet;
	
	return;
}



void Vertex::RemoveFromNeighbors()
{
	VertexNeighborSet::iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) 
		(*iNs)->ReplaceVertex(this, NULL);
	
	return;
}



void Vertex::CorrectNeighbors(Array<Element*> &pE)
{
	VertexNeighborSet oldSet = mNeighborSet;
	
	EraseNeighbors();
	
	VertexNeighborSet::iterator iNs;
	for (iNs = oldSet.begin(); iNs != oldSet.end(); ++iNs) {
		AddNeighbor(pE[(*iNs)->ID()]);
	}
	
	return;
}



void Vertex::CorrectNeighborVertex(const Vertex *pV)
{
	VertexNeighborSet::iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if ((*iNs)->IsAVertex(this) == false)
			(*iNs)->ReplaceVertex(pV, this);
	}

	return;
}



void Vertex::CommonNeighbors(Array<Element*> &common, const Vertex *pV0, const Vertex *pV1, ElementType eType) const
{
	VertexNeighborSet commonSet;
	CommonNeighbors(commonSet, pV0, pV1, eType);
	common = commonSet;
	
	return;
}
	


void Vertex::CommonNeighbors(VertexNeighborSet &common, const Vertex *pV0, const Vertex *pV1, ElementType eType) const
{	
	common = mNeighborSet;
	
	if (pV0 != NULL) 
		pV0->NeighborIntersection(common, eType);
	
	if (pV1 != NULL) 
		pV1->NeighborIntersection(common, eType);

	return;
}



void Vertex::NeighborIntersection(VertexNeighborSet &common, ElementType eType) const
{
	VertexNeighborSet tmp;	
	insert_iterator<VertexNeighborSet > tmpInsert(tmp, tmp.begin());
    
    // tmp is the intersection of mNeighborSet and common          
	set_intersection(mNeighborSet.begin(), mNeighborSet.end(), common.begin(), common.end(), tmpInsert);
		
	// remove all elements of the wrong type
	if (eType != ALL_ELEMENT_TYPES) {
		VertexNeighborSet::iterator iTmp = tmp.begin();
		while (iTmp != tmp.end()) {
			VertexNeighborSet::iterator iPrev = iTmp;
			++iTmp;
			if ((*iPrev)->Type() != eType) 
				tmp.erase(iPrev);
		}
	}
	
	common = tmp;
	
	
	return;
}



void Vertex::EraseNeighbors(ElementType eType)
{
	if (eType == ALL_ELEMENT_TYPES) {
		mNeighborSet.clear();
		return;
	}
	
	VertexNeighborSet::iterator iNs = mNeighborSet.begin();
	while (iNs != mNeighborSet.end()) {
		VertexNeighborSet::iterator iPrev = iNs;
		++iNs;
		
		if ((*iPrev)->Type() == eType) 
			mNeighborSet.erase(iPrev);
	}
	
	return;
}



Element* Vertex::CommonEdge(const Vertex *pV, bool throwException) const
{
	VertexNeighborSet::const_iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if ((*iNs)->Type() == EDGE) {
			Vertex *pVOpp = ((*iNs)->PV(0) == this) ? (*iNs)->PV(1) : (*iNs)->PV(0);
			if (pVOpp == pV)
				return *iNs;
		}
	}
	
    if (throwException)
        ThrowException("Vertex::CommonEdge : edge not found");	
        
	return NULL;
}



void Vertex::AddNeighborsToSet(set<Element*> &elementSet, ElementType eType) const
{
	VertexNeighborSet::const_iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if (((*iNs)->Type() == eType) || (eType == ALL_ELEMENT_TYPES))
			elementSet.insert(*iNs);
	}
	
	return;
}



void Vertex::GetNeighbors(Array<Element*> &pNeighbor, ElementType eType) const
{
    list<Element*> neighborList;
    VertexNeighborSet::const_iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if (((*iNs)->Type() == eType) || (eType == ALL_ELEMENT_TYPES))
            neighborList.push_back(*iNs);
	}
    
    pNeighbor = neighborList;
    
    return;
}



bool Vertex::IsANeighbor(const Element &element) const
{
    VertexNeighborSet::const_iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if ((*iNs)->Type() == element.Type()) {
            bool found = false;
            short i = 0;
            while (!found && (i < element.NumVertices())) {
                found = !(*iNs)->IsAVertex(element.PV(i));
                ++i;
            }
            if (found == false)
                return true;
        }
	}
    
    return false;
}



double Vertex::AverageEdgeLength() const
{
    double sum = 0.0;
    short count = 0;
    
    VertexNeighborSet::const_iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if ((*iNs)->Type() == EDGE) {
            Edge *pEdge = (Edge*) *iNs;
			sum += pEdge->Length();
            ++count;
        }
	}
    
    if (count == 0)
        ThrowException("Vertex::AverageEdgeLength : no edge neighbors");
    
    return sum / count;
}

		

void Vertex::CheckNeighbors() const
{
	VertexNeighborSet::const_iterator iNs;
	for (iNs = mNeighborSet.begin(); iNs != mNeighborSet.end(); ++iNs) {
		if ((*iNs)->Type() != BLOB) {
			if ((*iNs)->IsAVertex(this) == false)
				ThrowException("Vertex::CheckNeighbors : bad neighbor set");
		}
	}	

	return;
}



void Vertex::CheckListIterator() const
{
	if (mpListIterator == NULL)
		return;
			
	list<Vertex>::iterator iEtmp = *(list<Vertex>::iterator *) mpListIterator;
	Vertex *pTmp = &*iEtmp;
	if (pTmp != this)
		ThrowException("Vertex::CheckListIterator : bad vertex iterator");
		
	return;
}



bool Vertex::operator==(const Vertex &v) const
{
	ThrowException("Vertex::operator== called");
	return false;
}



bool Vertex::operator<(const Vertex &v) const
{
	ThrowException("Vertex::operator< called");
	return false;
}



bool Vertex::operator!=(const Vertex &v) const
{
	ThrowException("Vertex::operator!= called");
	return false;
}



bool Vertex::operator>(const Vertex &v) const
{
	ThrowException("Vertex::operator> called");
	return false;
}
