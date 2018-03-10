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
 
#include <iostream>
#include "meshgenerator.h"
#include "constants.h"
#include "matrix.h"
#include "randomnumbergenerator.h"

using namespace std;
using namespace NAMESPACE;


void MeshGenerator::MeshRectangularSolidWithTetrahedra(double xLength, double yLength, double zLength, short level)
{
	if ((xLength <= 0.0) || (yLength <= 0.0) || (zLength <= 0.0))
		ThrowException("MeshGenerator::MeshRectangularSolid : non-positive side length");
		
	if (level < 1)
		ThrowException("MeshGenerator::MeshRectangularSolid : level must be positive");
		
	double length[3] = {xLength, yLength, zLength};
	short minIndex = 0;
	double minLength = length[0];
	for (short i = 1; i < 3; ++i) {
		if (length[i] < minLength) {
			minLength = length[i];
			minIndex = i;
		}
	}
	
	double xDivision = xLength, yDivision, zDivision;
	long xPts, yPts, zPts;
	switch (minIndex) {
	case 0:
		xPts = level + 1;
		xDivision = xLength / level;
		
		yPts = NearestInteger(yLength / xDivision) + 1;
		yDivision = yLength / (yPts - 1.0);
		
		zPts = NearestInteger(zLength / xDivision) + 1;
		zDivision = zLength / (zPts - 1.0);
		break;
		
	case 1:
		yPts = level + 1;
		yDivision = yLength / level;
		
		xPts = NearestInteger(xLength / yDivision) + 1;
		xDivision = xLength / (yPts - 1.0);
		
		zPts = NearestInteger(zLength / yDivision) + 1;
		zDivision = zLength / (zPts - 1.0);
		break;
		
	case 2:
		zPts = level + 1;
		zDivision = zLength / level;
		
		xPts = NearestInteger(xLength / zDivision) + 1;
		xDivision = zLength / (xPts - 1.0);
		
		yPts = NearestInteger(yLength / zDivision) + 1;
		yDivision = yLength / (yPts - 1.0);
		break;
	}	
	
	// make vertices
	Array<Array<Array<Vertex*> > > pV(xPts);
	for (short i = 0; i < xPts; ++i) {
		pV[i].SetSize(yPts);
		for (short j = 0; j < yPts; ++j) {
			pV[i][j].SetSize(zPts);
		}
	}
	
	Vertex vTmp;
	Point p0(0.0, 0.0, 0.0);
	double x, y, z;
	for (short i = 0; i < xPts; ++i) {
		x = p0.X() + xDivision * i;
		for (short j = 0; j < yPts; ++j) {
			y = p0.Y() + yDivision * j;
			for (short k = 0; k < zPts; ++k) {
				z = p0.Z() + zDivision * k;
				vTmp.SetPosition(x, y, z);
				pV[i][j][k] = AddVertex(vTmp);
			}
		}
	}
	
	// make tetrahedra
    mDefaultElementType = TETRAHEDRON;
	Vertex *pVCube[8];
	for (short i = 0; i < xPts - 1; ++i) {
		for (short j = 0; j < yPts - 1; ++j) {
			for (short k = 0; k < zPts - 1; ++k) {
				pVCube[0] = pV[i  ][j  ][k  ];
				pVCube[1] = pV[i+1][j  ][k  ];
				pVCube[2] = pV[i+1][j+1][k  ];
				pVCube[3] = pV[i  ][j+1][k  ];
				pVCube[4] = pV[i  ][j  ][k+1];
				pVCube[5] = pV[i+1][j  ][k+1];
				pVCube[6] = pV[i+1][j+1][k+1];
				pVCube[7] = pV[i  ][j+1][k+1];
				
				MeshStandardCube(pVCube);
			}
		}
	}
	
	SetElementNeighbors();
	
	// add blob neighbor to represent exterior of mesh
	AddExteriorBlob();
	
	// set material number
	SetMaterialNumber(0);
	
	
	return;
}



void MeshGenerator::MeshStandardCube(Vertex *pV[8])
{
	// creates tetrahedral mesh of the rectangular solid where
	// pV[0] is the "lower left" corner and pV[6] is the "upper right" corner
	
	Tetrahedron tet;
	tet.SetVertices(pV[0], pV[1], pV[3], pV[5]);
	AddElement(&tet);
	
	tet.SetVertices(pV[0], pV[5], pV[3], pV[7]);
	AddElement(&tet);
	
	tet.SetVertices(pV[0], pV[5], pV[7], pV[4]);
	AddElement(&tet);
	
	tet.SetVertices(pV[1], pV[2], pV[3], pV[6]);
	AddElement(&tet);
	
	tet.SetVertices(pV[1], pV[3], pV[5], pV[6]);
	AddElement(&tet);
	
	tet.SetVertices(pV[3], pV[7], pV[5], pV[6]);
	AddElement(&tet);
	
	return;
}



void MeshGenerator::MeshRectangularSolidWithBricks(double xLength, double yLength, double zLength, 
        										   short xDivisions, short yDivisions, short zDivisions)
{
	if ((xLength <= 0.0) || (yLength <= 0.0) || (zLength <= 0.0))
		ThrowException("MeshGenerator::MeshRectangularSolid : non-positive side length");
		
	if (xDivisions < 1)
		ThrowException("MeshGenerator::MeshRectangularSolid : number of divisions must be positive");
	
	if (yDivisions <= 0)
		yDivisions = xDivisions;
	
	if (zDivisions <= 0)
		zDivisions = xDivisions;
	
	double dx = xLength / xDivisions;
	double dy = yLength / yDivisions;
	double dz = zLength / zDivisions;
	
	short xPts = xDivisions + 1;
	short yPts = yDivisions + 1;
	short zPts = zDivisions + 1;
	
	// make vertices
	Array<Array<Array<Vertex*> > > pV(xPts);
	for (short i = 0; i < xPts; ++i) {
		pV[i].SetSize(yPts);
		for (short j = 0; j < yPts; ++j) {
			pV[i][j].SetSize(zPts);
		}
	}
	
	Vertex vTmp;
	Point p0(0.0, 0.0, 0.0);
	double x, y, z;
	for (short i = 0; i < xPts; ++i) {
		x = p0.X() + dx * i;
		for (short j = 0; j < yPts; ++j) {
			y = p0.Y() + dy * j;
			for (short k = 0; k < zPts; ++k) {
				z = p0.Z() + dz * k;
				vTmp.SetPosition(x, y, z);
				pV[i][j][k] = AddVertex(vTmp);
			}
		}
	}
	
	// make bricks
    mDefaultElementType = BRICK;
    Brick brick;
	for (short i = 0; i < xPts - 1; ++i) {
		for (short j = 0; j < yPts - 1; ++j) {
			for (short k = 0; k < zPts - 1; ++k) {
				Vertex *pV0 = pV[i  ][j  ][k  ];
				Vertex *pV1 = pV[i+1][j  ][k  ];
				Vertex *pV2 = pV[i+1][j+1][k  ];
				Vertex *pV3 = pV[i  ][j+1][k  ];
				Vertex *pV4 = pV[i  ][j  ][k+1];
				Vertex *pV5 = pV[i+1][j  ][k+1];
				Vertex *pV6 = pV[i+1][j+1][k+1];
				Vertex *pV7 = pV[i  ][j+1][k+1];
				
				brick.SetVertices(pV0, pV1, pV2, pV3, pV4, pV5, pV6, pV7);
				AddElement(&brick);
			}
		}
	}
	
	SetElementNeighbors();
	
	// add blob neighbor to represent exterior of mesh
	AddExteriorBlob();
	
	// set material number
	SetMaterialNumber(0);
	
	
	return;
}



void MeshGenerator::MeshSphereWithBricks(double radius, double brickSize)
{
	double diameter = 2.0 * radius;
	short divisions = NearestInteger(diameter / brickSize);
	if (divisions < 1)
		ThrowException("MeshGenerator::MeshSphereWithBricks : brick size too large");

	double dx = diameter / divisions;
	short numPts = divisions + 1;

	// make vertices
	Array<Array<Array<Vertex*> > > pV(numPts);
	for (short i = 0; i < numPts; ++i) {
		pV[i].SetSize(numPts);
		for (short j = 0; j < numPts; ++j) {
			pV[i][j].SetSize(numPts);
		}
	}
	
	Vertex vTmp;
	Point p0(0.0, 0.0, 0.0);
	double x, y, z;
	for (short i = 0; i < numPts; ++i) {
		x = p0.X() + dx * i;
		for (short j = 0; j < numPts; ++j) {
			y = p0.Y() + dx * j;
			for (short k = 0; k < numPts; ++k) {
				z = p0.Z() + dx * k;
				vTmp.SetPosition(x, y, z);
				pV[i][j][k] = AddVertex(vTmp);
			}
		}
	}
	
	// make bricks
    mDefaultElementType = BRICK;
    Brick brick;
    Point centroid, center(radius, radius, radius);
	for (short i = 0; i < numPts - 1; ++i) {
		for (short j = 0; j < numPts - 1; ++j) {
			for (short k = 0; k < numPts - 1; ++k) {
				Vertex *pV0 = pV[i  ][j  ][k  ];
				Vertex *pV1 = pV[i+1][j  ][k  ];
				Vertex *pV2 = pV[i+1][j+1][k  ];
				Vertex *pV3 = pV[i  ][j+1][k  ];
				Vertex *pV4 = pV[i  ][j  ][k+1];
				Vertex *pV5 = pV[i+1][j  ][k+1];
				Vertex *pV6 = pV[i+1][j+1][k+1];
				Vertex *pV7 = pV[i  ][j+1][k+1];
				
				brick.SetVertices(pV0, pV1, pV2, pV3, pV4, pV5, pV6, pV7);
				
				brick.Centroid(centroid);
				if ((center - centroid).Norm() <= radius)
					AddElement(&brick);
			}
		}
	}
	
	SetElementNeighbors();
	
	// add blob neighbor to represent exterior of mesh
	AddExteriorBlob();
	
	// set material number
	SetMaterialNumber(0);

	return;
}



void MeshGenerator::MeshSphere(double radius, short level)
{    
    MeshGenerator base;
    base.MeshHalfCircle(radius, level);    
    ConstructByRotation(base);
    
	// set material number
	SetMaterialNumber(0);
	
	return;
}



void MeshGenerator::MeshCylinder(double radius, double length, short level)
{
	MeshGenerator base;
    base.MeshRectangle(radius, length, level);    
    ConstructByRotation(base);

	// set material number
	SetMaterialNumber(0);
	
	return;
}



void MeshGenerator::MeshPill(double waistLength, double capRadius, short numLevels)
{
	MeshGenerator base;
	base.MeshHalfOval(waistLength, capRadius, numLevels);

	ConstructByRotation(base);

	// set material number
	SetMaterialNumber(0);
	
	return;
}



void MeshGenerator::MeshCircle(double radius, short numLevels)
{
	// mesh unit circle
	
	// make vertices
	Vertex *pVCenter = AddVertex();
	pVCenter->SetPosition(0.0, 0.0, 0.0);
	Array<CircularArray<Vertex*> > pV(numLevels + 1);
	for (short i = 1; i < numLevels + 1; ++i) {
		pV[i].SetSize(6 * i);
		double radius = ((double) i) / ((double) numLevels);
		for (short j = 0; j < pV[i].Size(); ++j) {
			pV[i][j] = AddVertex();
			double theta = 2.0 * PI * j / pV[i].Size();
			pV[i][j]->SetPosition(radius * cos(theta), radius * sin(theta), 0.0);
		}
	}
	
	// make first level
    mDefaultElementType = TRIANGLE;
	Element *pTriangle;
	short i;
	for (i = 0; i < pV[1].Size(); ++i) {
		pTriangle = AddElement(TRIANGLE);
		pTriangle->SetVertices(pVCenter, pV[1][i], pV[1][i+1]);
	}
	
	// add levels
	for (short l = 2; l < numLevels + 1; ++l) {
		// this adds the "inner triangles"
		short k = 0, j;
		for (j = 0; j < pV[l-1].Size(); ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % (l-1) == 0)
				++k;
			pTriangle->SetVertices(pV[l-1][j], pV[l-1][j+1], pV[l][k]);
			++k;
		}
		
		// this adds the "outer triangles"
		k = 1;
		for (j = 0; j < pV[l].Size(); ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % l == 0)
				--k;
			pTriangle->SetVertices(pV[l][j], pV[l][j+1], pV[l-1][k]);
			++k;
		}
	}
	
	// set element neighbors
	SetElementNeighbors();
		
	// add blob
	AddExteriorBlob();
	
	// set material number
	SetMaterialNumber(0);
	
	// scale up to radius
	Scale(radius);
	
	return;
}



void MeshGenerator::MeshHalfCircle(double radius, short numLevels)
{
	// mesh unit circle
	// make vertices
	Vertex *pVCenter = AddVertex();
	pVCenter->SetPosition(0.0, 0.0, 0.0);
	Array<Array<Vertex*> > pV(numLevels + 1);
	for (short i = 1; i < numLevels + 1; ++i) {
		pV[i].SetSize(3 * i + 1);
		double radius = ((double) i) / ((double) numLevels);
		for (short j = 0; j < pV[i].Size(); ++j) {
			pV[i][j] = AddVertex();
			double theta = PI * j / (pV[i].Size() - 1.0) - 0.5 * PI;
			pV[i][j]->SetPosition(radius * cos(theta), radius * sin(theta), 0.0);
		}
	}
	
	// make first level
    mDefaultElementType = TRIANGLE;
	Element *pTriangle;
	short i;
	for (i = 0; i < pV[1].Size() - 1; ++i) {
		pTriangle = AddElement(TRIANGLE);
		pTriangle->SetVertices(pVCenter, pV[1][i], pV[1][i+1]);
	}
	
	// add levels
	for (short l = 2; l < numLevels + 1; ++l) {
		// this adds the "inner triangles"
		short k = 0, j;
		for (j = 0; j < pV[l-1].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % (l-1) == 0)
				++k;
			pTriangle->SetVertices(pV[l-1][j], pV[l-1][j+1], pV[l][k]);
			++k;
		}
		
		// this adds the "outer triangles"
		k = 1;
		for (j = 0; j < pV[l].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % l == 0)
				--k;
			pTriangle->SetVertices(pV[l][j], pV[l][j+1], pV[l-1][k]);
			++k;
		}
	}
	
	// set element neighbors
	SetElementNeighbors();
		
	// add blob
	AddExteriorBlob();
	
	// set material number
	SetMaterialNumber(0);
	
	// scale up to radius
	Scale(radius);
	
	return;
}



void MeshGenerator::MeshHalfOval(double waistLength, double capRadius, short numLevels)
{
	// mesh top quarter circle
	short nCap = (numLevels + 1) * (numLevels + 2) / 2;
	short nWaist = nCap - 2;
	Array<Array<Vertex*> > pV(nCap);
	Array<Vertex*> pVTop(numLevels + 1);
	for (short i = 0; i < numLevels + 1; ++i) {
		pV[i].SetSize(i + 1);
		double radius = capRadius * ((double) i) / numLevels;
		for (short j = 0; j < pV[i].Size(); ++j) {
			pV[i][j] = AddVertex();
			
			double theta;
			if (pV[i].Size() == 1)
				theta = 0.0;
			else
				theta = 0.5 * PI * j / (pV[i].Size() - 1.0);
			
			pV[i][j]->SetPosition(radius * cos(theta), radius * sin(theta), 0.0);
			
			if (j == 0)
				pVTop[i] = pV[i][j];
		}
	}
	
	// add levels
    mDefaultElementType = TRIANGLE;
	Element *pTriangle = NULL;
	short i = 0;
	for (short l = 0; l < numLevels; ++l) {
		// this adds the "inner triangles"
		short k = 0, j;
		for (j = 0; j < pV[l].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % (l) == 0)
				++k;
			pTriangle->SetVertices(pV[l][j], pV[l][j+1], pV[l+1][k]);
			++k;
		}
		
		// this adds the "outer triangles"
		k = 1;
		for (j = 0; j < pV[l+1].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % (l+1) == 0)
				--k;
			pTriangle->SetVertices(pV[l+1][j], pV[l+1][j+1], pV[l][k]);
			++k;
		}
	}
	
	Point trans(0.0, waistLength, 0.0);
	Translate(trans);
	
	// add bottom quarter circle
	Array<Vertex*> pVBottom(numLevels + 1);
	for (short i = 0; i < numLevels + 1; ++i) {
		pV[i].SetSize(i + 1);
		double radius = capRadius * ((double) i) / numLevels;
		for (short j = 0; j < pV[i].Size(); ++j) {
			pV[i][j] = AddVertex();
			
			double theta;
			if (pV[i].Size() == 1)
				theta = 0.0;
			else
				theta = 0.5 * PI * j / (pV[i].Size() - 1.0);
			
			pV[i][j]->SetPosition(radius * cos(theta), -radius * sin(theta), 0.0);
			
			if (j == 0)
				pVBottom[i] = pV[i][j];
		}
	}
	
	// add levels
	for (short l = 0; l < numLevels; ++l) {
		// this adds the "inner triangles"
		short k = 0, j;
		for (j = 0; j < pV[l].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % (l) == 0)
				++k;
			pTriangle->SetVertices(pV[l][j], pV[l][j+1], pV[l+1][k]);
			++k;
		}
		
		// this adds the "outer triangles"
		k = 1;
		for (j = 0; j < pV[l+1].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			if (j % (l+1) == 0)
				--k;
			pTriangle->SetVertices(pV[l+1][j], pV[l+1][j+1], pV[l][k]);
			++k;
		}
	}
	
	// add waist
	double xDivision = capRadius / (numLevels);
	
	if (waistLength / xDivision <= 0.5)
		ThrowException("MeshGenerator::MeshHalfOval : waist is too small, either increase waist or number of levels");
	
	short yLevels = NearestInteger(1.0 + waistLength / xDivision);
	double yDivision = waistLength / (yLevels - 1.0);
	
	pV.SetSize(yLevels);
	pV[0].SetSize(numLevels + 1);
	pV[yLevels - 1].SetSize(numLevels + 1);
	for (short j = 0; j < numLevels + 1; ++j) {
		pV[yLevels - 1][j] = pVTop[j];
		pV[0][j] = pVBottom[j];
	}
	for (short i = 1; i < yLevels - 1; ++i) {
		double y = i * yDivision;
		pV[i].SetSize(numLevels + 1);
		for (short j = 0; j < pV[i].Size(); ++j) {
			pV[i][j] = AddVertex();
			double x = j * xDivision;
			pV[i][j]->SetPosition(x, y, 0.0);
		}
	}
	
	for (short i = 0; i < pV.Size() - 1; ++i) {
		for (short j = 0; j < pV[i].Size() - 1; ++j) {
			pTriangle = AddElement(TRIANGLE);
			pTriangle->SetVertices(pV[i][j], pV[i+1][j], pV[i+1][j+1]);
			
			pTriangle = AddElement(TRIANGLE);
			pTriangle->SetVertices(pV[i][j], pV[i+1][j+1], pV[i][j+1]);
		}
	}
	
	// set neighbors
	SetElementNeighbors();
	
	// add blob
	AddExteriorBlob();

	// set material number
	SetMaterialNumber(0);
	
	
	return;
}



void MeshGenerator::MeshRectangle(double xLength, double yLength, short level)
{
	if ((xLength <= 0.0) || (yLength <= 0.0))
		ThrowException("MeshGenerator::MeshRectangle : non-positive side length");
		
	if (level < 1)
		ThrowException("MeshGenerator::MeshRectangle : level must be positive");
	
	//double minLength = (xLength < yLength) ? xLength : yLength;
	short minIndex = (xLength < yLength) ? 0 : 1;

	double xDivision = xLength, yDivision;
	long xPts, yPts;
	switch (minIndex) {
	case 0:
		xPts = level + 1;
		xDivision = xLength / level;
		
		yPts = NearestInteger(yLength / xDivision) + 1;
		yDivision = yLength / (yPts - 1.0);
		break;
		
	case 1:
		yPts = level + 1;
		yDivision = yLength / level;
		
		xPts = NearestInteger(xLength / yDivision) + 1;
		xDivision = xLength / (xPts - 1.0);
		break;
	}	
	
	// make vertices
	Matrix<Vertex*> pV(xPts, yPts);
	
	Vertex vTmp;
	Point p0(0.0, 0.0, 0.0);
	double x, y;
	for (short i = 0; i < xPts; ++i) {
		x = p0.X() + xDivision * i;
		for (short j = 0; j < yPts; ++j) {
			y = p0.Y() + yDivision * j;
			vTmp.SetPosition(x, y, 0.0);
			pV(i, j) = AddVertex(vTmp);
		}
	}
	
	// make triangles
    mDefaultElementType = TRIANGLE;
	Triangle triangle;
	for (short i = 0; i < xPts - 1; ++i) {
		for (short j = 0; j < yPts - 1; ++j) {
			triangle.SetVertices(pV(i, j), pV(i+1, j), pV(i+1, j+1));
			AddElement(&triangle);
			
			triangle.SetVertices(pV(i, j), pV(i+1, j+1), pV(i, j+1));
			AddElement(&triangle);
		}
	}
	
	SetElementNeighbors();
	
	// add blob
	AddExteriorBlob();

	// set material number
	SetMaterialNumber(0);
	
	return;
}



void MeshGenerator::ConstructByRotation(MeshGenerator &base)
{
	// assumptions coming into this function:
	// i)	base mesh is in XY plane
	// ii)	axis of rotation is Y axis
	
	// construct vertices
	// here the vertices are rotated about the Y axis to create vertices in 3D 
	// there is an implicit transformation between cylindrical and cartesian coordinates
	// in the cylindrical coordinates the base mesh X is the radius and the base mesh Y
	// is the "Z"
    
	// construct vertices
    base.MakeEdgeList();
	Array<CircularArray<Vertex*> > pV;
	MakeVerticesByRotation(base, pV);
	MakeEdgesByRotation(base, pV);
	MakeTetrahedraFromEdgeList();
    
    // set default element
    mDefaultElementType = TETRAHEDRON;
    
    // erase edges
    base.EraseEdgeList();
    EraseEdgeList();
    
	SetElementNeighbors();
	
	// add blob neighbor to represent exterior of mesh
	AddExteriorBlob();
    
	return;
}



void MeshGenerator::MakeVerticesByRotation(MeshGenerator &base, Array<CircularArray<Vertex*> > &pV)
{
	base.SetVertexIDs();
	
	pV.SetSize(base.NumVertices());
    
	Vertex vTmp;
	long id = 0;
	
	Array<Vertex*> pBaseVertex;
	base.GetVertices(pBaseVertex);
	
	for (long i = 0; i < pBaseVertex.Size(); ++i) {
		double distMeasure = pBaseVertex[i]->AverageEdgeLength();
		long baseID = pBaseVertex[i]->ID();
		double circumference = 2.0 * PI * pBaseVertex[i]->X();
		long number3DVertices = 1;
		if (circumference / distMeasure > 1.0) 
			number3DVertices = NearestInteger(circumference / distMeasure);
        
		//if ((*iVl).NeighborOnYZPlane() && n > 7) 
		//	n = 7;
        
		//if ((*iVl).NeighborOnYZPlane() && n < 5) 
		//	n = 5;
		
		pV[baseID].SetSize(number3DVertices);
		
		if (number3DVertices == 1) {
			// in this case the new vertex is on the Z axis
			vTmp.SetPosition(0.0, 0.0, pBaseVertex[i]->Y());
			vTmp.SetID(id);
			//vTmp.SetOnBoundaryType((*iVl).GetOnBoundaryType());
			pV[baseID][0] = AddVertex(vTmp);
		}
		else {
			for (short j = 0; j < number3DVertices; ++j) {
				// this is an angle around the 3D Z axis
				double angle = j * 2.0 * PI / number3DVertices;
				double x = pBaseVertex[i]->X() * cos(angle);
				double y = pBaseVertex[i]->X() * sin(angle);
				double z = pBaseVertex[i]->Y();
				vTmp.SetPosition(x, y, z);
				vTmp.SetID(id);
				//vTmp.SetOnBoundaryType((*iVl).GetOnBoundaryType());
				pV[baseID][j] = AddVertex(vTmp);
			}
		}
		
		++id;
	}
    
    
	return;
}



void MeshGenerator::MakeEdgesByRotation(MeshGenerator &base, Array<CircularArray<Vertex*> > &pV)
{	
	vector<Vertex*> pBoundary;
	list<Edge>::iterator iEl;
	
	Array<Element*> pBaseEdge;
	base.GetElements(pBaseEdge, EDGE);
	
	for (long i = 0; i < pBaseEdge.Size(); ++i) {
		long baseID[2];
		baseID[0] = pBaseEdge[i]->PV(0)->ID();
		baseID[1] = pBaseEdge[i]->PV(1)->ID();
        
        if (pV[baseID[0]].Size() < pV[baseID[1]].Size())
            StitchConicStrip(pV[baseID[0]], pV[baseID[1]]);
        else
            StitchConicStrip(pV[baseID[1]], pV[baseID[0]]);
    }
	
    SetVertexElementNeighbors(EDGE);
	
	return;
}



void MeshGenerator::StitchConicStrip(CircularArray<Vertex*> &pVSmaller, CircularArray<Vertex*> &pVLarger) 
{
    Edge edge;
    
    // case where both vertices are on axis
    if (pVLarger.Size() == 1) {
        edge.SetVertices(pVLarger[0], pVSmaller[0]);
        AddElement(&edge);
        return;
    }
    
    // case where only one vertex is on axis
    if (pVSmaller.Size() == 1) {
        for (long i = 0; i < pVLarger.Size(); ++i) {
            edge.SetVertices(pVLarger[i], pVSmaller[0]);
            AddElement(&edge);
            
            edge.SetVertices(pVLarger[i], pVLarger[i+1]);
            AddElement(&edge);
        }     
        
        return;   
    }
    
    // general case
    for (long i = 0; i < pVLarger.Size(); ++i) {
        edge.SetVertices(pVLarger[i], pVLarger[i+1]);
        AddElement(&edge);
    }
    
    for (long i = 0; i < pVSmaller.Size(); ++i) {
        edge.SetVertices(pVSmaller[i], pVSmaller[i+1]);
        AddElement(&edge);
    }
    
    long currentSmallerIndex;
    long prevSmallerIndex;
    for (long i = 0; i <= pVLarger.Size(); ++i) {
        if (i == 0) {
            currentSmallerIndex = 0;
            prevSmallerIndex = 0;
        }
        else {
            double d2Current = pVLarger[i]->SqrDistanceToPoint(pVSmaller[currentSmallerIndex]);
            double d2Next = pVLarger[i]->SqrDistanceToPoint(pVSmaller[currentSmallerIndex+1]);
            if (d2Next <= d2Current) {
                prevSmallerIndex = currentSmallerIndex;
                ++currentSmallerIndex;
            }
            else {
                prevSmallerIndex = currentSmallerIndex;
            }
        }
        
        
        edge.SetVertices(pVLarger[i], pVSmaller[currentSmallerIndex]);
        AddElement(&edge);
        
        // split quad
        if (currentSmallerIndex > prevSmallerIndex) {
            double d2Back = pVLarger[i]->SqrDistanceToPoint(pVSmaller[prevSmallerIndex]);
            double d2Up = pVLarger[i-1]->SqrDistanceToPoint(pVSmaller[currentSmallerIndex]);
            if (d2Back <= d2Up) 
                edge.SetVertices(pVLarger[i], pVSmaller[prevSmallerIndex]);
            else 
                edge.SetVertices(pVLarger[i-1], pVSmaller[currentSmallerIndex]);
            
            AddElement(&edge);
        }
    }
    
    
    return;
}



void MeshGenerator::MakeTetrahedraFromEdgeList()
{
    Tetrahedron tetrahedron;
    Array<Element*> pEdge;
    list<Vertex>::iterator iVl;
    for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
        (*iVl).GetNeighbors(pEdge, EDGE);
        for (short i = 0; i < pEdge.Size(); ++i) {
            Vertex *pVi = ((Edge*) pEdge[i])->OppositeVertex(&*iVl);
            for (short j = 0; j < i; ++j) {
                Vertex *pVj = ((Edge*) pEdge[j])->OppositeVertex(&*iVl);
                for (short k = 0; k < j; ++k) {
                    Vertex *pVk = ((Edge*) pEdge[k])->OppositeVertex(&*iVl);
                    if ((pVi->CommonEdge(pVj, false) != NULL) && (pVj->CommonEdge(pVk, false) != NULL) && (pVi->CommonEdge(pVk, false) != NULL)) {
                        tetrahedron.SetVertices(&*iVl, pVi, pVj, pVk);
                        if ((*iVl).IsANeighbor(tetrahedron) == false) {
                            Element *pNewTet = AddElement(&tetrahedron);
                            (*iVl).AddNeighbor(pNewTet);
                            pVi->AddNeighbor(pNewTet);
                            pVj->AddNeighbor(pNewTet);
                            pVk->AddNeighbor(pNewTet);
                        }
                    }
                }
            }
        }
    }
    
    // erase vertex element neighbors
    EraseVertexElementNeighbors();
    
    
	return;
}



void MeshGenerator::MakeTrianglesFromEdgeList()
{
    Triangle triangle;
    Array<Element*> pEdge;
    list<Vertex>::iterator iVl;
    for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
        (*iVl).GetNeighbors(pEdge, EDGE);
        for (short i = 0; i < pEdge.Size(); ++i) {
            Vertex *pVi = ((Edge*) pEdge[i])->OppositeVertex(&*iVl);
            for (short j = 0; j < i; ++j) {
                Vertex *pVj = ((Edge*) pEdge[j])->OppositeVertex(&*iVl);
                if (pVi->CommonEdge(pVj, false) != NULL) {
                    triangle.SetVertices(&*iVl, pVi, pVj);
                    if ((*iVl).IsANeighbor(triangle) == false) {
                        Element *pNewTriangle = AddElement(&triangle);
                        (*iVl).AddNeighbor(pNewTriangle);
                        pVi->AddNeighbor(pNewTriangle);
                        pVj->AddNeighbor(pNewTriangle);
                    }
                }
            }
        }
    }
    
	return;
}



void MeshGenerator::MeshLineEvenly(double x0, double x1, long numIntervals)
{
	// make vertices
	long xPts = numIntervals + 1;
	Array<Vertex*> pV(xPts);
	
	Vertex vTmp;
	double xDivision = 1.0 / numIntervals;
	for (long i = 0; i < xPts - 1; ++i) {
		vTmp.SetPosition(xDivision * i, 0.0, 0.0);
		pV[i] = AddVertex(vTmp);
	}
	vTmp.SetPosition(1.0, 0.0, 0.0);
	pV[xPts - 1] = AddVertex(vTmp);
	
	// scale
	for (long i = 0; i < xPts; ++i) 
		pV[i]->SetPosition(x0 + (x1 - x0) * pV[i]->X(), 0.0, 0.0);
	
	// make edges
    mDefaultElementType = EDGE;
	Edge edge;
	for (long i = 0; i < numIntervals; ++i) {
		edge.SetVertices(pV[i], pV[i+1]);
		AddElement(&edge);
	}
	
	// set edge neighbors
	Array<Element*> pE;
	GetElements(pE);
	for (long i = 0; i < pE.Size(); ++i) 
		pE[i]->AllocateNeighborArray();

	for (long i = 0; i < pE.Size() - 1; ++i) 
		pE[i]->SetNeighbor(0, pE[i+1]);
		
	for (long i = 1; i < pE.Size(); ++i) 
		pE[i]->SetNeighbor(1, pE[i-1]);
		
	/*
	Array<Element*> pE;
	GetElements(pE);
	for (long i = 1; i < pE.Size() - 1; ++i) {
		pE[i]->AllocateNeighborArray();
		pE[i]->SetNeighbor(1, pE[i-1]);
		pE[i]->SetNeighbor(0, pE[i+1]);
	}
	pE[0]->AllocateNeighborArray();
	pE[0]->SetNeighbor(0, pE[1]);
	pE[numIntervals - 1]->AllocateNeighborArray();
	pE[numIntervals - 1]->SetNeighbor(1, pE[numIntervals - 2]);
	*/
	
	// add blob
	AddExteriorBlob();

	// set material number
	SetMaterialNumber(0);
	
	return;
}



void MeshGenerator::MeshLineRandomly(double x0, double x1, long numIntervals)
{
	RandomNumberGenerator rng;
	
	// make vertices
	long xPts = numIntervals + 1;
	Array<Vertex*> pV(xPts);
	
	Vertex vTmp(0.0, 0.0, 0.0);
	pV[0] = AddVertex(vTmp);
	for (long i = 1; i < xPts; ++i) {
		double x = vTmp.X() + rng.RandomNumber(0.1, 1.0);
		vTmp.SetPosition(x, 0.0, 0.0);
		pV[i] = AddVertex(vTmp);
	}
	
	double total = pV[xPts - 1]->X();
	for (long i = 0; i < xPts - 1; ++i) 
		*pV[i] /= total;
	pV[xPts - 1]->SetPosition(1.0, 0.0, 0.0);
		
	// scale
	for (long i = 0; i < xPts; ++i) 
		pV[i]->SetPosition(x0 + (x1 - x0) * pV[i]->X(), 0.0, 0.0);
		
	// make edges
    mDefaultElementType = EDGE;
	Edge edge;
	for (long i = 0; i < numIntervals; ++i) {
		edge.SetVertices(pV[i], pV[i+1]);
		AddElement(&edge);
	}
	
	// set edge neighbors
	Array<Element*> pE;
	GetElements(pE);
	for (long i = 1; i < pE.Size() - 1; ++i) {
		pE[i]->AllocateNeighborArray();
		pE[i]->SetNeighbor(1, pE[i-1]);
		pE[i]->SetNeighbor(0, pE[i+1]);
	}
	pE[0]->AllocateNeighborArray();
	pE[0]->SetNeighbor(0, pE[1]);
	pE[numIntervals - 1]->AllocateNeighborArray();
	pE[numIntervals - 1]->SetNeighbor(1, pE[numIntervals - 2]);
	
	// add blob
	AddExteriorBlob();

	// set material number
	SetMaterialNumber(0);
	

	return;
}



void MeshGenerator::MakeMeshFromPolygonFile(string fileName)
{
	ifstream file;
    OpenFile(fileName, file);
			
	// read number of polygons
	long numPolygons;
	file >> numPolygons;

	Array<Array<Point> > polygon(numPolygons);
	short numVertices;
	double x, y;
	Point p;
	for (long i = 0; i < numPolygons; ++i) {
		file >> numVertices;
		polygon[i].SetSize(numVertices);
		for (short j = 0; j < numVertices; ++j) {
			file >> x;
			file >> y;
			p.SetPosition(x, y, 0.0);
			polygon[i][j] = p;
		}
		
		Element *pEtmp = AddElement(POLYGON);
		pEtmp->AllocateVertexArray(numVertices);
		pEtmp->SetMaterialNumber(0);
	}
	file.close();
	
	mDefaultElementType = POLYGON;
	
	Array<Element*> pE;
	GetElements(pE);
	
	for (long i = 0; i < pE.Size(); ++i) {	
		for (short j = 0; j < pE[i]->NumVertices(); ++j) {
			Vertex *pV = IsAVertex(polygon[i][j]);
			if (pV == NULL)
				pV = AddVertex(polygon[i][j]);
			
			pE[i]->SetVertex(j, pV);
		}
	}
	
	SetElementNeighbors();
	AddExteriorBlob();
	
	CheckAll();

	return;
}
