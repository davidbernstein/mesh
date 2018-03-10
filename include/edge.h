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
 
#ifndef _edge_h_
#define _edge_h_

#include "element.h"

namespace NAMESPACE {
	class Edge : public Element {
	public:
		// Constructors
		Edge(void);

		// Destructor
        ~Edge(void) {DeleteListIterator();};
        void DeleteListIterator(void);
                
		// type
		ElementType Type(void) const {return EDGE;};
		
		// vertex functions
		short NumVertices(void) const {return 2;};
        Vertex* OppositeVertex(const Vertex *pV) const;
		
		// neighbors
		short NumNeighbors(void) const {return 2;};
		
		// edges
		short NumEdges(void) const {return 1;};
		void GetEdge(short index, Vertex *pV[2]) const;
	
		// geometry
		double Length(void) const;
		Point MidPoint(void) const;
		short Dimension(void) const {return 1;};
		double Volume(void) const;
        void FaceArea(Array<double> &faceArea) const;
        void GetFace(short index, Array<Vertex*> &pV) const;
        bool ContainsPoint(const Point &p) const;
        void OutwardAreaNormal(Array<Point> &unitNormal) const;
        bool IntersectsEdge(Vertex *pV[2], Point &p) const;

		// bisection
		void Bisect(Element *pChild[2], Vertex *pVNew);
		
		// boundary
		bool OnBoundary(void) const;
		
		// list iterator
		void MakeListIterator(void *pListIterator);
	
		// operators needed by STL
		bool operator==(const Edge &e) const;
		bool operator<(const Edge &e) const;
		bool operator!=(const Edge &e) const;
		bool operator>(const Edge &e) const;
	};
	
	
	
	inline Edge::Edge()
	{
		AllocateVertexArray();
		return;
	}
	


	inline double Edge::Length() const
	{
		return PV(0)->DistanceToPoint(PV(1));
	}
	
	
	
	inline double Edge::Volume() const
	{
		return Length();
	}
	
	
	
    inline void Edge::FaceArea(Array<double> &faceArea) const
    {
    	faceArea.SetSize(2);
    	faceArea[0] = 1.0;
    	faceArea[1] = 1.0;
    	
    	return;
    }
    
	
	
	inline void Edge::GetFace(short index, Array<Vertex*> &pV) const
	{
		// faces are stored in the following order:
		// face[0] is opposite to vertex 0
		// face[1] is opposite to vertex 1
		
		pV.SetSize(1);
		
        switch (index) {
        case 0:
            pV[0] = PV(1);
            break;
        
        case 1:
            pV[0] = PV(0);
            break;
       
        default:
            ThrowException("Edge::GetFace : bad face index");
            break;
        }
        
		return;
	}
	
	
	inline void Edge::GetEdge(short index, Vertex *pV[2]) const
	{
		if (index == 0) {
			pV[0] = PV(0);
			pV[1] = PV(1);
		}
		else {
			ThrowException("Edge::GetEdge : bad index");
        }
		
		return;
	}



    inline Vertex* Edge::OppositeVertex(const Vertex *pV) const
    {
        if (pV == PV(0))
            return PV(1);
        
        if (pV == PV(1))
            return PV(0);
        
        ThrowException("Edge::OppositeVertex : not a vertex");
        return NULL;
    }
        
	
	
	inline void Edge::MakeListIterator(void *pListIterator)
	{
        if (mpListIterator != NULL) 
            ThrowException("Edge::MakeListIterator : already on a list");
                
        std::list<Edge>::iterator *pLi = new std::list<Edge>::iterator;
        *pLi = *(std::list<Edge>::iterator *) pListIterator;
        mpListIterator = (std::list<Edge>::iterator *) pLi;

		return;
	}
        
        
        
    inline void Edge::DeleteListIterator()
    {
        delete ((std::list<Edge>::iterator *) mpListIterator);
        return;
    }
}


#endif // _edge_h_

