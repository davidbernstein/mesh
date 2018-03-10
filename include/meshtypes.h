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
 
#ifndef _meshtypes_h_
#define _meshtypes_h_

#include <set>
#include "namespace.h"

namespace NAMESPACE {
	class Element;
	class Vertex;
}

namespace NAMESPACE {
    class RefinementQueueRelation {	
    public:
        bool operator()(const Element *pE0, const Element *pE1) const;
    };

    class UnRefinementQueueRelation {	
    public:
        bool operator()(const Element *pE0, const Element *pE1) const;
    };
}

typedef std::pair<NAMESPACE::Vertex*, NAMESPACE::Vertex*> VertexPair;
typedef std::set<NAMESPACE::Element*> VertexNeighborSet;
typedef std::set<NAMESPACE::Element*, NAMESPACE::RefinementQueueRelation> RefinementQueue;
typedef std::set<NAMESPACE::Element*, NAMESPACE::UnRefinementQueueRelation> UnRefinementQueue;

#endif // _meshtypes_h_
