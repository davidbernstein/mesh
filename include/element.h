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
 
#ifndef _element_h_
#define _element_h_

#include "array.h"
#include "vertex.h"
#include "datapacket.h"

namespace NAMESPACE {
	class Element : public DataPacket {
	public:
		Element(void);
		virtual ~Element(void);
		Element(const Element &e);
		Element& operator=(const Element &e);

		// id
		void SetID(long id);
		long ID(void) const;
		
		// pure virtual functions
		virtual ElementType Type(void) const = 0;
		virtual short NumVertices(void) const = 0;
		virtual short NumNeighbors(void) const = 0;
		virtual short NumEdges(void) const = 0;
		virtual void MakeListIterator(void *pListIterator) = 0;
		virtual short Dimension(void) const = 0;
        				   
		// vertices
		void AllocateVertexArray(short n = 0);
		Vertex* PV(short index) const;
		void SetVertex(short index, const Vertex *pV);
		void SetVertices(const Vertex *pV0 = NULL, const Vertex *pV1 = NULL, const Vertex *pV2 = NULL, const Vertex *pV3 = NULL,
						 const Vertex *pV4 = NULL, const Vertex *pV5 = NULL, const Vertex *pV6 = NULL, const Vertex *pV7 = NULL); 
		void AddVerticesToSet(std::set<Vertex*> &vertexSet) const;
		void UpdateVertexNeighbors(void);
		void RemoveFromVertexNeighborSets(void);
		bool IsAVertex(const Vertex *pV) const;
		short WhichVertex(const Vertex *pV) const;
		void ReplaceVertex(const Vertex *pVOld, const Vertex *pVNew);
		void CorrectVertices(const Array<Vertex*> &pV);
        
		// geometry functions
        virtual void GetEdge(short index, Vertex *pV[2]) const = 0;
        virtual short EdgeIndex(const Vertex *pV0, const Vertex *pV1) const;
        virtual void GetFace(short index, Array<Vertex*> &pV) const;
        virtual short FaceIndex(const Vertex *pV0, const Vertex *pV1 = NULL, const Vertex *pV2 = NULL) const;
		void Centroid(Point &centroid) const;
		short LongestEdge(void) const;
		virtual double Volume(void) const;	
        virtual double Area(void) const;
        virtual void FaceArea(Array<double> &faceArea) const;
        virtual void OutwardAreaNormal(Array<Point> &unitNormal) const;
		virtual void Orient(void);
		virtual bool Oriented(void) const;
		virtual Vertex* EdgeMidPoint(short index) const;
        virtual bool ContainsPoint(const Point &p) const;
        
		// element neighbors
		void AllocateNeighborArray(void);
		Element* PN(short index) const;
		void SetNeighbors(void);
		void SetNeighbor(short index, const Element *pN);
		bool IsANeighbor(const Element *pE) const;
		short NeighborIndex(const Element *pE) const;
		void ReplaceNeighbor(const Element *pOld, const Element *pNew, bool checkForNeighbor = false);
		void Ring(Array<Element*> &ring, short edgeIndex) const;
		bool OnBoundary(void) const;
		void Neighborhood(std::set<Element*> &neighborhood, short neighborhoodSize, ElementType eType) const;
		
		// material number
		void SetMaterialNumber(short m);
		short MaterialNumber(void) const;
        
		// refinement
        virtual void Bisect(Element *pChild[2], Vertex *pVNew);
		bool Conforming(void) const;
        virtual Vertex* Merge(Element *pSibling, Element *pParent);
                
        // list iterator
        void* GetListIterator(void) const;
		
		// useful IO functions
		void WriteVertexIDsToFile(std::ofstream &file) const;
		void WriteToFile(const std::string &fileName) const;
		void WriteToFile(std::ofstream &file) const;
        void Print(std::string prefix = "") const;

		// checking
		void CheckVertices(void) const;
		void CheckNeighbors(void) const;
        
    private:
        // useful IO functions
        void PrintVertices(void) const;
        void PrintNeighbors(void) const;
        
	protected:
		// id number of element
		long mID;
		
		// material number
		short mMaterialNumber;
		
		// vertices 
		Array<Vertex*> mpVertex;
		
		// neighbors
		Array<Element*> mpNeighbor;
                
		// position on a list
		void *mpListIterator;
	};
	
	
	
	inline Element::Element()
	{
		mID = -1;
		mMaterialNumber = -1;
		mpListIterator = NULL;

		return;
	}



	inline void Element::AllocateVertexArray(short n)
	{
		short numVertices = NumVertices();
		if (numVertices > -1)
			mpVertex.SetSize(numVertices);
		else
			mpVertex.SetSize(n);
		
		for (short i = 0; i < mpVertex.Size(); ++i)
			mpVertex[i] = NULL;
			
		return;
	}
	


	inline void Element::AllocateNeighborArray()
	{
		mpNeighbor.SetSize(NumNeighbors());
		
		for (short i = 0; i < mpNeighbor.Size(); ++i)
			mpNeighbor[i] = NULL;
			
		return;
	}
	


	inline void Element::SetID(long id)
	{
		mID = id;
		return;
	}



	inline long Element::ID() const
	{
		return mID;
	}
	
	
	
	inline void Element::SetMaterialNumber(short m)
	{
		mMaterialNumber = m;
		return;
	}
	
	
	
	inline short Element::MaterialNumber() const
	{
		return mMaterialNumber;
	}
	
	
	
	inline Vertex* Element::PV(short index) const
	{
		return mpVertex[index];
	}
	
	
	
	inline void Element::SetVertex(short index, const Vertex *pV)
	{
		mpVertex[index] = (Vertex *) pV;
		return;
	}
	
	
	
	inline void Element::SetVertices(const Vertex *pV0, const Vertex *pV1, const Vertex *pV2, const Vertex *pV3,
						 			 const Vertex *pV4, const Vertex *pV5, const Vertex *pV6, const Vertex *pV7)
	{
		// covers all cases up to Brick
		
		if (pV0 != NULL)
			mpVertex[0] = (Vertex *) pV0;
			
		if (pV1 != NULL)
			mpVertex[1] = (Vertex *) pV1;
			
		if (pV2 != NULL)
			mpVertex[2] = (Vertex *) pV2;
			
		if (pV3 != NULL)
			mpVertex[3] = (Vertex *) pV3;
		
		if (pV4 != NULL)
			mpVertex[4] = (Vertex *) pV4;
			
		if (pV5 != NULL)
			mpVertex[5] = (Vertex *) pV5;
			
		if (pV6 != NULL)
			mpVertex[6] = (Vertex *) pV6;
			
		if (pV7 != NULL)
			mpVertex[7] = (Vertex *) pV7;
		
		return;
	}



	inline bool Element::IsAVertex(const Vertex *pV) const
	{
		return WhichVertex(pV) != -1;
	}



	inline short Element::WhichVertex(const Vertex *pV) const
	{
		for (short i = 0; i < mpVertex.Size(); ++i) {
			if (mpVertex[i] == pV)
				return i;
		}

		return -1;
	}

	
	
	inline Element* Element::PN(short index) const
	{
		return mpNeighbor[index];
	}
	
	
	
	inline void Element::SetNeighbor(short index, const Element *pN)
	{
		mpNeighbor[index] = (Element *) pN;
		return;
	}
	
	
	
	inline void Element::ReplaceNeighbor(const Element *pOld, const Element *pNew, bool checkForNeighbor)
	{
		for (short i = 0; i < mpNeighbor.Size(); ++i) {
			if (mpNeighbor[i] == pOld) {
				mpNeighbor[i] = (Element*) pNew;
				return;
			}
		}
        
        if (checkForNeighbor == true) {
        	//BisectionData *pBD = pOld->GetBisectionData();
        	//bool a = IsANeighbor(pBD->Parent());
        	//bool b = pBD->Parent()->IsANeighbor(this);
        	PrintNeighbors();
        	pOld->PrintNeighbors();
            ThrowException("Element::ReplaceNeighbor : neighbor not found");
       	}
		
		return;
	}



	inline bool Element::IsANeighbor(const Element *pE) const
	{
		for (short i = 0; i < mpNeighbor.Size(); ++i) {
			if (mpNeighbor[i] == pE)
				return true;
		}
	
		return false;
	}



	inline short Element::NeighborIndex(const Element *pE) const
	{
		for (short i = 0; i < mpNeighbor.Size(); ++i) {
			if (mpNeighbor[i] == pE)
				return i;
		}
	
		return -1;
	}
	
	
	
	inline void Element::Bisect(Element *pChild[2], Vertex *pVNew)
	{
		ThrowException("Element::Bisect : element type not supported");
		return;
	}
	
	
	
	inline short Element::FaceIndex(const Vertex *pV0, const Vertex *pV1, const Vertex *pV2) const
	{
		ThrowException("Element::FaceIndex : element type not supported");
		return -1;
	}
	
	
	
	inline short Element::EdgeIndex(const Vertex *pV0, const Vertex *pV1) const
	{
		ThrowException("Element::EdgeIndex : element type not supported");
		return -1;
	}
	
	
	
	inline void Element::Orient()
	{
		ThrowException("Element::Orient : element type not supported");
		return;
	}
	
	
	
	inline bool Element::Oriented() const
	{
		ThrowException("Element::Oriented : element type not supported");
		return false;
	}
	
	
	
	inline Vertex* Element::EdgeMidPoint(short index) const
	{
		ThrowException("Element::EdgeMidPoint : element type not supported");
		return NULL;
	}
    
    
    
    inline bool Element::ContainsPoint(const Point &p) const
    {
        ThrowException("Element::ContainsPoint : element type not supported");
		return false;

    }
    
    
    
    inline double Element::Volume() const
    {
        // elements which have volume should override this function
        ThrowException("Element::Volume : volume not defined for element type " + ConvertIntegerToString(Type()));
        return -1.0;
    }
    
    
    
    inline double Element::Area(void) const
    {
        // elements which have area should override this function
        ThrowException("Element::Area : not defined for element type " + ConvertIntegerToString(Type()));
        return -1.0;
    }
    
    
    
    inline void Element::FaceArea(Array<double> &faceArea) const
    {
        ThrowException("Element::FaceArea : not defined for element type " + ConvertIntegerToString(Type()));
        return;
    }
    
    
    
    inline void Element::OutwardAreaNormal(Array<Point> &unitNormal) const
    {
        ThrowException("Element::OutwardUnitNormal : not defined for element type " + ConvertIntegerToString(Type()));
        return;
    }
    
    
    
    inline Vertex* Element::Merge(Element *pSibling, Element *pParent)
    {
        ThrowException("Element::Merge : element type not supported");
        return NULL;
    }
    
        
          
    inline void* Element::GetListIterator() const
    {
        return mpListIterator;
    }
    
    
    
    inline bool Element::OnBoundary(void) const
    {
        for (short i = 0; i < mpNeighbor.Size(); ++i) {
            if (mpNeighbor[i]->Type() == BLOB)
                return true;
        }
        
        return false;
    }
    
    
    /*
    inline void Element::GetEdge(short index, Vertex *pV[2]) const
    {   
        ThrowException("Element::GetEdge : element type not supported");
        return;
    }
    */
    
    
    inline void Element::GetFace(short index, Array<Vertex*> &pV) const
    {   
        ThrowException("Element::GetFace : element type not supported");
        return;
    }
}


#endif // _element_h_
