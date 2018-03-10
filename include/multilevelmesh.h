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
 
#ifndef _multilevelmesh_h_
#define _multilevelmesh_h_

#include "mesh.h"

namespace NAMESPACE {
	class MultiLevelMesh : public Mesh {	
    public:
        MultiLevelMesh(void) { };
        virtual ~MultiLevelMesh(void) { };
        
		// refinement
		void RefineToConformity(RefinementQueue &refQ);
		void RefineToConformity(Element *pE); 
		void Bisect(Element *pE);
		Vertex* RefinementEdgeMidPoint(const Element *pE);
		void SetBisectionData(void);
        void SetRefinementListIDs(long value);
        
        // unrefinement
        void UnRefineToConformity(UnRefinementQueue &unRefQ);
		void UnRefineToConformity(Element *pE, UnRefinementQueue &unRefQ);
		bool UnBisect(Element *pE, UnRefinementQueue &unRefQ);
        
        // bisection data
        void UpdateBisectionData(Element *pE);
        
        // testing
        void TestRefinement(void);
        void TestUnRefinement(void);
        void RefineAll(short level = 1);
        void PrintRefinementQueue(short numTop = 0) const;
        void PrintUnRefinementQueue(short numTop = 0) const;
        
    protected:
		// refinement queue
		RefinementQueue mRefinementQueue;
        
        // unrefinement queue
        UnRefinementQueue mUnRefinementQueue;
	};
}


#endif // _multilevelmesh_h_


