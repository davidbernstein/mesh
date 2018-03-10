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

#ifndef _vertex_h_
#define _vertex_h_

#include "point.h"
#include "array.h"
#include "meshenums.h"
#include "meshtypes.h"

namespace NAMESPACE {
	class Element;
}

namespace NAMESPACE {
	class Vertex : public Point {	
	public:
		// Constructor
		Vertex(double x = 0.0, double y = 0.0, double z = 0.0);
		Vertex(const Point &p);

		// Copy constructor
		Vertex(const Vertex &v);

		// Destructor
		~Vertex(void);

		// operators
		Vertex& operator=(const Vertex &v);
		Vertex& operator=(const Point &p);
		void Copy(const Vertex &v);
		
		// neighbor functions
		bool NumNeighbors(void) const;
        void GetNeighbors(Array<Element*> &pNeighbor, ElementType eType = ALL_ELEMENT_TYPES) const;
		void AddNeighbor(Element *pE);
		void RemoveNeighbor(const Element *pE);
		void ReplaceNeighbor(const Element *pEOld, const Element *pENew);
		void EraseNeighbors(ElementType eType = ALL_ELEMENT_TYPES);
		void AddNeighborsToSet(std::set<Element*> &elementSet, ElementType eType = ALL_ELEMENT_TYPES) const;
		bool IsANeighbor(const Element *pE) const;
        bool IsANeighbor(const Element &element) const;
		void CorrectNeighbors(Array<Element*> &pE);
		void RemoveFromNeighbors(void);
		void CommonNeighbors(Array<Element*> &common, const Vertex *pV0, const Vertex *pV1, ElementType eType) const;
		void CorrectNeighborVertex(const Vertex *pV);
		
		// list iterator
		void SetListIterator(void *pVertexListIterator);
		void* GetListIterator(void) const;
        
        // geometry
        double AverageEdgeLength(void) const;
		Element* CommonEdge(const Vertex *pV, bool throwException = true) const;

		// checking
		void CheckNeighbors(void) const;
		void CheckListIterator(void) const;

		// operators needed by STL
		bool operator==(const Vertex &v) const;
		bool operator<(const Vertex &v) const;
		bool operator!=(const Vertex &v) const;
		bool operator>(const Vertex &v) const;
		
		// private member functions
	private:
		void NeighborIntersection(VertexNeighborSet &common, ElementType type) const;
		void CommonNeighbors(VertexNeighborSet &common, const Vertex *pV0 = NULL, const Vertex *pV1 = NULL, 
							 ElementType eType = ALL_ELEMENT_TYPES) const;

	private:
		// list of element neighbors
		VertexNeighborSet mNeighborSet; 
		
		// this keeps track of the position of this vertex on a list
		void *mpListIterator;
	};



	inline Vertex::Vertex(double x, double y, double z)
	{
		mX = x;
		mY = y;
		mZ = z;

		mID = -1;
		mpListIterator = NULL;

		return;
	}



	inline Vertex::Vertex(const Point &p)
	{
		mX = p.X();
		mY = p.Y();
		mZ = p.Z();
		mID = p.ID();
		
		mpListIterator = NULL;
		
		return;
	}
	
	
	
	inline Vertex::~Vertex()
	{
            if (mpListIterator != NULL) {
                delete ((std::list<Vertex>::iterator *) mpListIterator);
                mpListIterator = NULL;
            }
		
            return;
	}
	
	
	
	inline bool Vertex::NumNeighbors() const
	{
		return mNeighborSet.size();
	}
	


	inline void Vertex::AddNeighbor(Element *pE)
	{
		mNeighborSet.insert(pE);
		return;
	}



	inline void Vertex::RemoveNeighbor(const Element *pE)
	{
        if (!mNeighborSet.empty())
            mNeighborSet.erase((Element*) pE);
            
		return;
	}



	inline void Vertex::ReplaceNeighbor(const Element *pEOld, const Element *pENew)
	{
		RemoveNeighbor(pEOld);
		AddNeighbor((Element *) pENew);
		return;
	}



	inline bool Vertex::IsANeighbor(const Element *pE) const
	{
		return mNeighborSet.find((Element*) pE) != mNeighborSet.end();
	}



	inline void* Vertex::GetListIterator() const
	{
        return mpListIterator;
	}
	
	
	
	inline void Vertex::SetListIterator(void *pListIterator)
	{
		if (mpListIterator != NULL)
			ThrowException("Vertex::SetListIterator : already on a list");
			
		mpListIterator = new std::list<Vertex>::iterator;
		*(std::list<Vertex>::iterator *) mpListIterator = *(std::list<Vertex>::iterator *) pListIterator;

		return;
	}
}


#endif // _vertex_h_

