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
 
#ifndef _blob_h_
#define _blob_h_

#include "element.h"

namespace NAMESPACE {
	class Blob : public Element {
	public:
		// constructors
		Blob(void) { };
        ~Blob(void) {DeleteListIterator();};

		// type
		ElementType Type(void) const {return BLOB;};
		short NumVertices(void) const {return 0;};
		
		// neighbors
		short NumNeighbors(void) const {return 0;};
		
		// edges
		short NumEdges(void) const {return 0;};
        void GetEdge(short index, Vertex *pV[2]) const {};
        
        // geometry
        short Dimension(void) const {return 0;};
        
		// list iterator
		void MakeListIterator(void *pListIterator);
		void DeleteListIterator(void);
                
		// operators needed by STL
		bool operator==(const Blob &b) const;
		bool operator<(const Blob &b) const;
		bool operator!=(const Blob &b) const;
		bool operator>(const Blob &b) const;
	};
	
	
    
    inline void Blob::MakeListIterator(void *pListIterator)
    {
        if (mpListIterator != NULL) 
            ThrowException("Blob::MakeListIterator : already on a list");
                
        std::list<Blob>::iterator *pLi = new std::list<Blob>::iterator;
        *pLi = *(std::list<Blob>::iterator *) pListIterator;
        mpListIterator = (std::list<Blob>::iterator *) pLi;

        return;
    }
    
    
        
    inline void Blob::DeleteListIterator()
    {
        delete ((std::list<Blob>::iterator *) mpListIterator);
        return;
    }
}


#endif // _blob_h_

