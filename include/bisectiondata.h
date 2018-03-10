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
 
#ifndef _bisectiondata_h_
#define _bisectiondata_h_

#include "dataobject.h"
#include "utility.h"

#include <string>

namespace NAMESPACE {
	class Element;
	class Vertex;
}

namespace NAMESPACE {
	class BisectionData : public DataObject {
	public:
		BisectionData(void);
		~BisectionData(void);
        BisectionData(const BisectionData &bd);
    	BisectionData& operator=(const BisectionData &bd);
		
		DataType Type(void) const {return BISECTION_DATA;}
		long Size(void) const {return -1;};
		void SetSize(short arraySize) {return;};
		
		// edge markings
		short RefinementEdge(void) const;
		void SetRefinementEdge(short e);
		
		// level
        short Level(void) const;
        void IncrementLevel(void);
		
		// parents and children
		void SetParent(const Element *pE);
        Element* Parent(void) const;
		void SetChildren(const Element *pE0, const Element *pE1);
		bool IsALeaf(void) const;
		bool IsRoot(void) const;
		void GetChildren(Element *pChild[2]) const;
		void EraseChildren(void);
		
		// midpoint vertices
		void AllocateMidpointArray(short numEdges);
		void SetMidpoint(short edgeIndex, const Vertex* pV);
		Vertex* RefinementEdgeMidPoint(void) const;
		Vertex* Midpoint(short edgeIndex) const;
		
	private:
		// edge markings
		short mRefinementEdge;
		
		// level
        short mLevel;
		
		// parents and children
		Element *mpParent;
		Element **mpChild;
		
		// edge midpoint vertices
		Array<Vertex*> mpMidpoint;
	};
	
	
	
	inline BisectionData::BisectionData()
	{
		mRefinementEdge = -1;
		mLevel = 0;
		
		mpParent = NULL;
		mpChild = NULL;
		
			
		return;
	}
	
	
	
	inline BisectionData::~BisectionData()
	{	
		EraseChildren();
		return;
	}
	
    
    
    inline BisectionData::BisectionData(const BisectionData &bd)
    {
        *this = bd;
        return;
    }
    
    
    
    inline BisectionData& BisectionData::operator=(const BisectionData &bd)
    {
        mRefinementEdge = bd.mRefinementEdge;
		mLevel = bd.mLevel;
		
		mpParent = bd.mpParent;
        
        if (bd.mpChild != NULL) 
            SetChildren(bd.mpChild[0], bd.mpChild[1]);
        else
            mpChild = NULL;
        
        mpMidpoint = bd.mpMidpoint;
        
        return *this;
    }
    
	
	
	inline void BisectionData::EraseChildren()
	{
		if (mpChild != NULL) {
			delete [] mpChild;
			mpChild = NULL;
		}
			
		return;
	}
    
    
    
    inline short BisectionData::Level() const
    {
        return mLevel;
    }
	
    
    
    inline void BisectionData::IncrementLevel()
    {
        ++mLevel;
        return;
    }
    
	
	
	inline short BisectionData::RefinementEdge() const
	{
		return mRefinementEdge;
	}
	
	
	
	inline void BisectionData::SetRefinementEdge(short e)
	{
		mRefinementEdge = e;
		return;
	}
	
	
	
	inline void BisectionData::SetParent(const Element *pE)
	{
		mpParent = (Element*) pE;
		return;
	}
	
    
    
    inline Element* BisectionData::Parent() const
    {
        return mpParent;
    }
    
	
	
	inline void BisectionData::SetChildren(const Element *pE0, const Element *pE1)
	{
		if (mpChild == NULL)
			mpChild = new Element*[2];
		
		mpChild[0] = (Element*) pE0;
		mpChild[1] = (Element*) pE1;
		
		return;
	}
	
	
	
	inline void BisectionData::GetChildren(Element *pChild[2]) const
	{
		pChild[0] = mpChild[0];
		pChild[1] = mpChild[1];
		
		return;
	}
	
	
	
	inline bool BisectionData::IsALeaf() const
	{
		return mpChild == NULL;
	}
	
	
	
	inline bool BisectionData::IsRoot() const
	{
		return mpParent == NULL;
	}
	
	
	
	inline void BisectionData::AllocateMidpointArray(short numEdges)
	{
		mpMidpoint.SetSize(numEdges);
		
		for (short i = 0; i < numEdges; ++i)
			mpMidpoint[i] = NULL;
			
		return;
	}
	
	
	
	inline void BisectionData::SetMidpoint(short edgeIndex, const Vertex* pV)
	{
		if (mpMidpoint[edgeIndex] == NULL)
			mpMidpoint[edgeIndex] = (Vertex*) pV;
		else
			ThrowException("BisectionData::SetMidpoint : midpoint already set");
		
		return;
	}
	
	
	
	inline Vertex* BisectionData::Midpoint(short edgeIndex) const
	{
		return mpMidpoint[edgeIndex];
	}
	
	
	
	inline Vertex* BisectionData::RefinementEdgeMidPoint() const
	{
		return mpMidpoint[mRefinementEdge];
	}
}


#endif // _bisectiondata_h_

