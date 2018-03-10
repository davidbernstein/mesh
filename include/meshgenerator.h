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
 
#ifndef _meshgenerator_h_
#define _meshgenerator_h_

#include "mesh.h"
#include "circulararray.h"

namespace NAMESPACE {
	class MeshGenerator : public Mesh {	
    public:
        MeshGenerator(void) { };
        ~MeshGenerator(void) { };
        
        // mesh generation routines for 3D domains
        void MeshRectangularSolidWithTetrahedra(double xLength, double yLength, double zLength, short level = 1);
        void MeshRectangularSolidWithBricks(double xLength, double yLength, double zLength, 
        									short xDivisions = 1, short yDivisions = -1, short zDivisions = -1);
        void MeshSphere(double radius, short level = 1);
        void MeshSphereWithBricks(double radius, double brickSize);
        void MeshCylinder(double radius, double length, short level = 1);
    	void MeshPill(double waistLength, double capRadius, short numLevels);
        
        // 2D domains
        void MeshRectangle(double xLength, double yLength, short level = 1);
    	void MeshCircle(double radius, short numLevels);
    	void MeshHalfCircle(double radius, short numLevels);
    	void MeshHalfOval(double waistLength, double capRadius, short numLevels);
    	
        // 1D domains
        void MeshLineEvenly(double x0, double x1, long numIntervals);
        void MeshLineRandomly(double x0, double x1, long numIntervals);
        
        // polygon files 
        void MakeMeshFromPolygonFile(std::string fileName);
        
    private:
        // mesh generation routines for simple domains
        void MeshStandardCube(Vertex *pV[8]);
		void ConstructByRotation(MeshGenerator &base);
        void MakeVerticesByRotation(MeshGenerator &base, Array<CircularArray<Vertex*> > &pV);
        void MakeEdgesByRotation(MeshGenerator &base, Array<CircularArray<Vertex*> > &pV);
        void MakeTrianglesFromEdgeList(void);
        void MakeTetrahedraFromEdgeList(void);
        void StitchConicStrip(CircularArray<Vertex*> &pVSmaller, CircularArray<Vertex*> &pVLarger);
	};
}


#endif // _meshgenerator_h_


