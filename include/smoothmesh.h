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
 
#ifndef _smoothmesh_h_
#define _smoothmesh_h_

#include "mesh.h"

namespace NAMESPACE {
	class SmoothMesh : public Mesh {	
    public:
        // mesh generation routines for simple domains
		void MeshRectangularSolid(double xLength, double yLength, double zLength, short level = 1);
        void MeshRectangle(double xLength, double yLength, short level = 1);
    	void MeshCircle(double radius, short numLevels);
        void ConstructByRotation(SmoothMesh &base);
        
    protected:
		// mesh generation routines for simple domains
        void MeshStandardCube(Vertex *pV[8]);
		void MeshUnitCircle(short numLevels);
		void AddExteriorBlob(ElementType eType);
        void MakeVerticesByRotation(SmoothMesh &base, Array<CircularArray<Vertex*> > &pV);
        void MakeEdgesByRotation(SmoothMesh &base, Array<CircularArray<Vertex*> > &pV);
        void MakeTrianglesFromEdgeList(void);
        void MakeTetrahedraFromEdgeList(void);
        void StitchConicStrip(CircularArray<Vertex*> &pVSmaller, CircularArray<Vertex*> &pVLarger);
	};
}


#endif // _smoothmesh_h_




