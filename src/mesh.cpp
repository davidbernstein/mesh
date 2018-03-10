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
#include <fstream>
#include "mesh.h"
#include "meshconstants.h"

using namespace std;
using namespace NAMESPACE;


// copy constructor
Mesh::Mesh(const Mesh &mesh)
{
	*this = mesh;
	return;
}



Mesh& Mesh::operator=(const Mesh &mesh)
{
	ThrowException("Mesh operator= currently broken, doesn't transfer bisection data correctly");
	
    // pMesh is non-const
	Mesh *pMesh = (Mesh*) &mesh;

	// erase old this
	Erase();

	// if mesh is empty then we're done
	if (mesh.NumVertices() == 0) 
		return *this;

	// save IDs of mesh Vertices and Elements
	Array<long> vertexID;
	pMesh->SaveVertexIDs(vertexID);
	
	Array<long> elementID;
	pMesh->SaveElementIDs(elementID);
	
	// set ids 
	pMesh->SetVertexIDs();
	pMesh->SetElementIDs();
	
	// copy lists (needs to be done manually because of list iterators)
	list<Vertex>::iterator iVl;
	for (iVl = pMesh->mVertexList.begin(); iVl != pMesh->mVertexList.end(); ++iVl)
		AddVertex(*iVl);
		
	Array<Element*> pE;
	pMesh->GetElements(pE);
	for (long i = 0; i < pE.Size(); ++i)
		AddElement(pE[i]);
		
	// correct vertex and element pointers
	Array<Vertex*> pV;
	GetVertices(pV);
	GetElements(pE);
	CorrectVertexNeighbors(pE);
	CorrectElementVertices(pV);

	// restore original IDs
	RestoreVertexIDs(vertexID);
	RestoreElementIDs(elementID);
	pMesh->RestoreVertexIDs(vertexID);
	pMesh->RestoreElementIDs(elementID);

	// primary element type
	mDefaultElementType = mesh.mDefaultElementType;


    return *this;
}



void Mesh::UpdateVertexNeighbors(Array<Element*> &pE)
{
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->UpdateVertexNeighbors();

	return;
}



void Mesh::SetVertexElementNeighbors(ElementType eType)
{
	// erase currents neighbors
	EraseVertexElementNeighbors();
	
    if (eType == DEFAULT_ELEMENT_TYPE)
        eType = mDefaultElementType;
        
	Array<Element*> pE;
	GetElements(pE, eType);
	for (long i = 0; i < pE.Size(); ++i) {
		pE[i]->UpdateVertexNeighbors();
    }
	
		
	return;
}



void Mesh::AddExteriorBlob(ElementType eType, short materialNumber)
{
    if (eType == DEFAULT_ELEMENT_TYPE)
        eType = mDefaultElementType;
    
	Array<Element*> pE;
	GetElements(pE, eType);
	
	Element *pOutside = AddElement(BLOB);
	pOutside->SetMaterialNumber(materialNumber);
	
	for (long i = 0; i < pE.Size(); ++i) {
		for (short j = 0; j < pE[i]->NumNeighbors(); ++j) {
			if (pE[i]->PN(j) == NULL)
		 		pE[i]->SetNeighbor(j, pOutside);
		}
	}
	
	return;
}



void Mesh::EraseVertexElementNeighbors()
{
	// erase currents neighbors
	list<Vertex>::iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl).EraseNeighbors();
	
	return;
}



void Mesh::SetElementNeighbors()
{
	SetVertexElementNeighbors();
	
	Array<Element*> pE;
	GetElements(pE);
	for (long i = 0; i < pE.Size(); ++i) {
		pE[i]->SetNeighbors();
    }
		
	EraseVertexElementNeighbors();
	
	
	return;
}



void Mesh::GetElements(Array<Element*> &pE, ElementType type, bool leavesOnly) const
{
	vector<Element*> v;
    v.reserve(NumElementsTotal());
	
	BisectionData *pBD;
    
    if (type == DEFAULT_ELEMENT_TYPE)
        type = mDefaultElementType;
	
	if ((type == EDGE) || (type == ALL_ELEMENT_TYPES)) {
		list<Edge>::const_iterator iEl;
		for (iEl = mEdgeList.begin(); iEl != mEdgeList.end(); ++iEl) {
			pBD = (*iEl).GetBisectionData(false);
			if (pBD != NULL) {
				if (pBD->IsALeaf())
					v.push_back((Element*) &*iEl);
			}
			else {
				v.push_back((Element*) &*iEl);
			}
		}
	}
	
	if ((type == TRIANGLE) || (type == ALL_ELEMENT_TYPES)) {
		list<Triangle>::const_iterator iTr;
		for (iTr = mTriangleList.begin(); iTr != mTriangleList.end(); ++iTr) {
			pBD = (*iTr).GetBisectionData(false);
			if (pBD != NULL) {
				if (pBD->IsALeaf())
					v.push_back((Element*) &*iTr);
			}
			else {
				v.push_back((Element*) &*iTr);
			}
		}
	}
	
	if ((type == POLYGON) || (type == ALL_ELEMENT_TYPES)) {
		list<Polygon>::const_iterator iPo;
		for (iPo = mPolygonList.begin(); iPo != mPolygonList.end(); ++iPo) {
			pBD = (*iPo).GetBisectionData(false);
			if (pBD != NULL) {
				if (pBD->IsALeaf())
					v.push_back((Element*) &*iPo);
			}
			else {
				v.push_back((Element*) &*iPo);
			}
		}
	}
	
	if ((type == TETRAHEDRON) || (type == ALL_ELEMENT_TYPES)) {
		list<Tetrahedron>::const_iterator iTe;
		for (iTe = mTetrahedronList.begin(); iTe != mTetrahedronList.end(); ++iTe) {
			pBD = (*iTe).GetBisectionData(false);
			if (pBD != NULL) {
				if (pBD->IsALeaf())
					v.push_back((Element*) &*iTe);
			}
			else {
				v.push_back((Element*) &*iTe);
			}
		}
	}
	
	if ((type == BRICK) || (type == ALL_ELEMENT_TYPES)) {
		list<Brick>::const_iterator iBe;
		for (iBe = mBrickList.begin(); iBe != mBrickList.end(); ++iBe) {
			pBD = (*iBe).GetBisectionData(false);
			if (pBD != NULL) {
				if (pBD->IsALeaf())
					v.push_back((Element*) &*iBe);
			}
			else {
				v.push_back((Element*) &*iBe);
			}
		}
	}
	
	if ((type == BLOB) || (type == ALL_ELEMENT_TYPES)) {
		list<Blob>::const_iterator iBl;
		for (iBl = mBlobList.begin(); iBl != mBlobList.end(); ++iBl) {
			pBD = (*iBl).GetBisectionData(false);
			if (pBD != NULL) {
				if (pBD->IsALeaf())
					v.push_back((Element*) &*iBl);
			}
			else {
				v.push_back((Element*) &*iBl);
			}
		}
	}	
	
	pE = v;
	
	
	return;
}



void Mesh::MakeEdgeList()
{
	// check if edge list is already active
	if (!mEdgeList.empty())
		return;
	
	Array<Element*> pE;
	GetElements(pE);
		
	set<VertexPair> edgeSet;
	for (long i = 0; i < pE.Size(); ++i) {
		Vertex *pV[2];
		for (short j = 0; j < pE[i]->NumEdges(); ++j) {
			pE[i]->GetEdge(j, pV);
			Vertex *pV0 = min(pV[0], pV[1]);
			Vertex *pV1 = max(pV[0], pV[1]);
			edgeSet.insert(make_pair(pV0, pV1));
		}
	}
	
	Edge eTmp;
	set<VertexPair>::iterator iEs;
	for (iEs = edgeSet.begin(); iEs != edgeSet.end(); ++iEs) {
		eTmp.SetVertices((*iEs).first, (*iEs).second);
		AddElement(&eTmp);
	}
	
	// add edges to vertex neighbor lists
	list<Edge>::iterator iEl;
	for (iEl = mEdgeList.begin(); iEl != mEdgeList.end(); ++iEl) 
		(*iEl).UpdateVertexNeighbors();
	

	return;
}



void Mesh::EraseEdgeList()
{
	// remove edges from vertex neighbor lists
	list<Vertex>::iterator iVl = mVertexList.begin();
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl).EraseNeighbors(EDGE);	
	
	// erase list
	mEdgeList.erase(mEdgeList.begin(), mEdgeList.end());
	
	return;
}



void Mesh::BoundingBox(Point &min, Point &max) const
{
	if (mVertexList.empty()) {
		ThrowException("Mesh::BoundingBox : no vertices");
		return;
	}

	list<Vertex>::const_iterator iVl = mVertexList.begin();
	min.SetPosition(*iVl);
	max.SetPosition(*iVl);
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		(*iVl).Min(min);
		(*iVl).Max(max);
	}


	return;
}



void Mesh::SetVertexIDs()
{
	list<Vertex>::iterator iVl;
	long i = 0;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		(*iVl).SetID(i);
		++i;
	}

	return;
}



void Mesh::SetVertexIDs(short value)
{
	list<Vertex>::iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl).SetID(value);

	return;
}



void Mesh::SetElementIDs(ElementType eType)
{
	Array<Element*> pE;
	GetElements(pE, eType);
	for (long i = 0; i < pE.Size(); ++i) 
		pE[i]->SetID(i);

	return;
}



void Mesh::SetElementIDs(short value, ElementType eType)
{
	Array<Element*> pE;
	GetElements(pE, eType);
	for (long i = 0; i < pE.Size(); ++i) 
		pE[i]->SetID(value);

	return;
}



void Mesh::Erase()
{
	mVertexList.erase(mVertexList.begin(), mVertexList.end());
	mEdgeList.erase(mEdgeList.begin(), mEdgeList.end());
	mTriangleList.erase(mTriangleList.begin(), mTriangleList.end());
	mPolygonList.erase(mPolygonList.begin(), mPolygonList.end());
	mTetrahedronList.erase(mTetrahedronList.begin(), mTetrahedronList.end());
	mBrickList.erase(mBrickList.begin(), mBrickList.end());
	mBlobList.erase(mBlobList.begin(), mBlobList.end());
	
	return;
}



Element* Mesh::AddElement(ElementType type)
{
	Edge edge;
	Triangle triangle;
	Tetrahedron tetrahedron;
	Brick brick;
	Blob blob;
	Polygon polygon;
	
	list<Edge>::iterator iEl;
	list<Triangle>::iterator iTr;
	list<Tetrahedron>::iterator iTe;
	list<Blob>::iterator iBl;
	list<Polygon>::iterator iPl;
	list<Brick>::iterator iBr;
	
	switch (type) {
	case EDGE:
		mEdgeList.push_back(edge);
		iEl = --mEdgeList.end();
		(*iEl).MakeListIterator(&iEl);
		return &*iEl;
		break;
		
	case TRIANGLE:
		mTriangleList.push_back(triangle);
		iTr = --mTriangleList.end();
		(*iTr).MakeListIterator(&iTr);
		return &*iTr;
		break;
		
	case POLYGON:
		mPolygonList.push_back(polygon);
		iPl = --mPolygonList.end();
		(*iPl).MakeListIterator(&iPl);
		return &*iPl;
		break;
		
	case TETRAHEDRON:
		mTetrahedronList.push_back(tetrahedron);
		iTe = --mTetrahedronList.end();
		(*iTe).MakeListIterator(&iTe);
		return &*iTe;
		break;
		
			
	case BRICK:
		mBrickList.push_back(brick);
		iBr = --mBrickList.end();
		(*iBr).MakeListIterator(&iBr);
		return &*iBr;
		break;
		
	case BLOB:
		mBlobList.push_back(blob);
		iBl = --mBlobList.end();
		(*iBl).MakeListIterator(&iBl);
		return &*iBl;
		break;
		
	default:
		ThrowException("Mesh::AddElement : bad element type");
		break;
	}
	
	return NULL;
}



Element* Mesh::AddElement(const Element *pE)
{
	Element *pENew = AddElement(pE->Type());
	
	pENew->SetID(pE->ID());
	
	for (short i = 0; i < pE->NumVertices(); ++i)
		pENew->SetVertex(i, pE->PV(i));
	
	return pENew;
}



void Mesh::EraseElement(Element *pE)
{
	// remove element from its vertex element neighbor lists
	pE->RemoveFromVertexNeighborSets();

	// erase element
	void *pListIterator = pE->GetListIterator();

	if (pListIterator == NULL) 
		ThrowException("Mesh::EraseElement : no element list iterator");
	
	switch (pE->Type()) {
	case EDGE:
		mEdgeList.erase(*(list<Edge>::iterator *) pListIterator);
		break;
		
	case TRIANGLE:
		mTriangleList.erase(*(list<Triangle>::iterator *) pListIterator);
		break;
			
	case POLYGON:
		mPolygonList.erase(*(list<Polygon>::iterator *) pListIterator);
		break;
		
	case TETRAHEDRON:
		mTetrahedronList.erase(*(list<Tetrahedron>::iterator *) pListIterator);
		break;
			
	case BRICK:
		mBrickList.erase(*(list<Brick>::iterator *) pListIterator);
		break;
		
	case BLOB:
		mBlobList.erase(*(list<Blob>::iterator *) pListIterator);
		break;
		
	default:
		ThrowException("Mesh::AddElement : bad element type");
		break;
	}
	
	return;
}



Vertex* Mesh::AddVertex(const Vertex &v)
{
	mVertexList.push_back(v);
	
	list<Vertex>::iterator iNew = --mVertexList.end();
	Vertex *pVNew = &*iNew;
	pVNew->SetListIterator(&iNew);

	return pVNew;
}



Vertex* Mesh::AddVertex()
{
	Vertex vTmp;
	return AddVertex(vTmp);
}



void Mesh::EraseVertex(Vertex *pV)
{
	if (pV == NULL)
		ThrowException("Mesh::EraseVertex : vertex is NULL");
		
	// remove this vertex from its element neighbors
	pV->RemoveFromNeighbors();
	
	// remove vertex
	void *pListIterator = pV->GetListIterator();

	if (pListIterator == NULL) 
            ThrowException("Mesh::EraseVertex : vertex not found");

	mVertexList.erase(*(list<Vertex>::iterator *) pListIterator);
	
	return;
}



Vertex* Mesh::IsAVertex(const Point &p) const
{
	list<Vertex>::const_iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		if ((*iVl).SameLocation(p))
			return (Vertex*) &*iVl;
	}
	
	return NULL;
}



void Mesh::Translate(const Point &trans)
{
	list<Vertex>::iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl) += trans;

	return;
}



void Mesh::Scale(double scaleFactor)
{
	list<Vertex>::iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl) *= scaleFactor;

	return;
}

		

void Mesh::GetVertices(Array<Vertex*> &pV)
{
	pV.SetSize(mVertexList.size());
	list<Vertex>::iterator iVl;
	long i = 0;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		pV[i] = &*iVl;
		++i;
	}
	
	return;
}



void Mesh::OrientElements(ElementType elementType)
{
	Array<Element*> pE;
	GetElements(pE);
	for (long i = 0; i < pE.Size(); ++i) {
		if (pE[i]->Type() == elementType)
			pE[i]->Orient();
	}
		
	return;
}


	
void Mesh::CorrectVertexNeighbors(Array<Element*> &pE) 
{
	list<Vertex>::iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl).CorrectNeighbors(pE);
	
	return;
}



void Mesh::CorrectElementVertices(const Array<Vertex*> &pV) 
{
	Array<Element*> pE;
	GetElements(pE);
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->CorrectVertices(pV);

	return;
}



void Mesh::SaveVertexIDs(Array<long> &id) const
{
	id.SetSize(mVertexList.size());
	list<Vertex>::const_iterator iVl;
	long i = 0;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		id[i] = (*iVl).ID();
		++i;
	}

	return;
}
		
		
		
void Mesh::RestoreVertexIDs(Array<long> &id)
{
	list<Vertex>::iterator iVl;
	long i = 0;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		(*iVl).SetID(id[i] );
		++i;
	}

	return;
}



void Mesh::SaveElementIDs(Array<long> &id) const
{
	Array<Element*> pE;
	GetElements(pE);
		
	id.SetSize(pE.Size());
	
	for (long i = 0; i < pE.Size(); ++i)
		id[i] = pE[i]->ID();
	
	return;
}
		
		
		
void Mesh::RestoreElementIDs(Array<long> &id)
{
	Array<Element*> pE;
	GetElements(pE);
	
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->SetID(id[i]);

	return;
}



long Mesh::NumElements(ElementType eType, bool leavesOnly) const
{	
	eType = (eType == DEFAULT_ELEMENT_TYPE) ? mDefaultElementType : eType;
	
	Array<Element*> pE;
	GetElements(pE, eType, leavesOnly);
	
	return pE.Size();
}



long Mesh::NumElementsTotal(ElementType eType) const
{	
	eType = (eType == DEFAULT_ELEMENT_TYPE) ? mDefaultElementType : eType;
	
	switch (eType) {
	case NO_ELEMENT_TYPE:
		return 0;
		break;
		
	case EDGE:
		return mEdgeList.size();
		break;
			
	case POLYGON:
		return mPolygonList.size();
		break;
			
	case TRIANGLE:
		return mTriangleList.size();
		break;
		
	case TETRAHEDRON:
		return mTetrahedronList.size();
		break;
		
	case BRICK:
		return mBrickList.size();
		break;
		
	case BLOB:
		return mBlobList.size();
		
	case ALL_ELEMENT_TYPES:
		return mEdgeList.size() + mTriangleList.size() + mTetrahedronList.size() + mBlobList.size();
		break;
	
	default:
		ThrowException("Mesh::NumElements : invalid element type");
		return -1;
	}
}



void Mesh::SetElementDataSize(short numData, ElementType elementType)
{
	Array<Element*> pE;
    GetElements(pE, elementType);
		
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->SetDataSize(numData);
		
	return;
}



void Mesh::AddElementData(DataType dataType, short arraySize, ElementType elementType)
{
	Array<Element*> pE;
    GetElements(pE, elementType);
		
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->AddDataObject(dataType, arraySize);

	return;
}



void Mesh::SetMaterialNumber(short m, ElementType elementType)
{
	Array<Element*> pE;
    GetElements(pE, elementType);
		
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->SetMaterialNumber(m);
		
	return;
}



double Mesh::Volume(ElementType elementType) const
{
	Array<Element*> pE;
    GetElements(pE, elementType);
		
	double sum = 0.0;
	for (long i = 0; i < pE.Size(); ++i)
		sum += pE[i]->Volume();
	
	return sum;
}



void Mesh::VolumeByMaterialNumber(Array<double> &volume, bool normalizeByTotal) const
{
	Array<Element*> pE;
    GetElements(pE);
	
	short numMaterials = CheckMaterialNumbers();
	
	volume.SetSize(numMaterials);
	for (short i = 0; i < numMaterials; ++i) 
		volume[i] = 0.0;
	
	double totalVolume = 0.0;
	for (long i = 0; i < pE.Size(); ++i) {
		double v = pE[i]->Volume();;
		volume[pE[i]->MaterialNumber()] += v;
		totalVolume += v;
	}
		
	for (short i = 0; i < numMaterials; ++i) 
		volume[i] /= totalVolume;
		
		
	return;
}



double Mesh::SurfaceArea() const
{
	if (Dimension() != 3)
		return 0.0;
	
	Array<Element*> pE;
    GetElements(pE);
	
	Array<Vertex*> pV;
	Triangle triangle;
	
	double area = 0.0;
	for (long i = 0; i < pE.Size(); ++i) {
		for (short j = 0; j < pE[i]->NumNeighbors(); ++j) {
			if (pE[i]->PN(j)->Type() == BLOB) {
				pE[i]->GetFace(j, pV);
				if (pV.Size() != 3)
					ThrowException("Mesh::SurfaceArea : face type not supported");
				
				triangle.SetVertices(pV[0], pV[1], pV[2]);
				area += triangle.Volume();
			}
		}
	}
	
	return area;
}



double Mesh::LongestEdgeLength() const
{
	Array<Element*> pE;
    GetElements(pE);
	
	Vertex *pV[2];
	
	double lMax = -1.0;
	for (long i = 0; i < pE.Size(); ++i) {
		for (short j = 0; j < pE[i]->NumEdges(); ++j) {
			pE[i]->GetEdge(j, pV);
			double length = pV[0]->DistanceToPoint(pV[1]);
			lMax = max(length, lMax);
		}
	}
	
	return lMax;
}



void Mesh::MaterialNumberTest()
{
	Array<Element*> pE;
    GetElements(pE);
		
	Point centroid;
	/*
	for (long i = 0; i < pE.Size(); ++i) {
		pE[i]->Centroid(centroid);
		
		
		//if (centroid.X() > 0.5)
		//	pE[i]->SetMaterialNumber(1);
		
		
		if ((centroid.X() > 2.0) && (centroid.X() < 3.0))
			pE[i]->SetMaterialNumber(1);
	}		
	*/
	
	Point center(5, 5, 5);
	bool found = false;
	long i = 0;
	while (!found) {
		pE[i]->Centroid(centroid);
		
		if ((centroid - center).Norm() < 1.5) {
			pE[i]->SetMaterialNumber(1);
			found = true;
		}
		
		++i;
	}	
	
	return;
}



void Mesh::WriteToFile(string fileName) const
{
	ofstream file(fileName.c_str());
	if (!file.good())
		ThrowException("Mesh::WriteToFile : was not able to create file " + fileName);
	
	Array<long> vertexID;
	SaveVertexIDs(vertexID);
	((Mesh*) this)->SetVertexIDs();
	
	Array<long> elementID;
	SaveElementIDs(elementID);
	((Mesh*) this)->SetElementIDs();
	
	// write vertices
	file << NumVertices() << endl;
	list<Vertex>::const_iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		file << (*iVl).ID() << " ";
		file << (*iVl).X() << " " << (*iVl).Y() << " " << (*iVl).Z() << endl;
	}
	
	// write elements
	Array<Element*> pE;
	GetElements(pE, ALL_ELEMENT_TYPES);	
	file << pE.Size() << endl;
	
	for (long i = 0; i < pE.Size(); ++i) {
		file << pE[i]->ID() << " ";
		file << pE[i]->Type() << " ";
		file << pE[i]->MaterialNumber() << " ";
		
		for (short j = 0; j < pE[i]->NumVertices(); ++j)
			file << pE[i]->PV(j)->ID() << " ";
			
		file << endl;	
	}
	
	// write element neighbors
	for (long i = 0; i < pE.Size(); ++i) {	
		for (short j = 0; j < pE[i]->NumNeighbors(); ++j)
			file << pE[i]->PN(j)->ID() << " ";
			
		file << endl;	
	}
	
	((Mesh*) this)->RestoreVertexIDs(vertexID);	
	((Mesh*) this)->RestoreElementIDs(elementID);	
		
	file.close();
	
	
	return;
}



void Mesh::ReadFromFile(string fileName)
{
	ifstream file;
    OpenFile(fileName, file);
		
	Erase();
	
	// read vertices
	long numVertices;
	file >> numVertices;
	Array<Vertex*> pV(numVertices);
	Vertex vTmp;
	double x, y, z;
	long id;
	for (long i = 0; i < numVertices; ++i) {
		file >> id;
		file >> x;
		file >> y;
		file >> z;
		
		vTmp.SetID(id);
		vTmp.SetPosition(x, y, z);
		pV[i] = AddVertex(vTmp);
	}
	
	// read elements
	long numElements;
	file >> numElements;
	short eType;
	long eID;
	short materialNumber;
	mDefaultElementType = ALL_ELEMENT_TYPES;
    Array<Element*> pE(numElements);
	for (long i = 0; i < numElements; ++i) {
		file >> eID;
		file >> eType;
		file >> materialNumber;
		
		ElementType type = (ElementType) eType;
		pE[i] = AddElement(type);
		
		pE[i]->SetID(eID);
		pE[i]->SetMaterialNumber(materialNumber);
		
		// read element vertices
		if (type == POLYGON) {
			char pLine[POLYGON_VERTEX_LINE_LENGTH];
			file.getline(pLine, POLYGON_VERTEX_LINE_LENGTH);
			string s = pLine;
			Array<long> a;
			ConvertStringToIntegerArray(s, a);
			
			pE[i]->AllocateVertexArray(a.Size());
			for (short j = 0; j < a.Size(); ++j) 
				pE[i]->SetVertex(j, pV[a[j]]);
		}
		else {
			short numVertices = pE[i]->NumVertices();
			for (short j = 0; j < numVertices; ++j) {
				file >> id;
				pE[i]->SetVertex(j, pV[id]);
			}
		}
		
		if (type != BLOB) {
			if (mDefaultElementType == ALL_ELEMENT_TYPES)
				mDefaultElementType = type;
			
			if (mDefaultElementType != type)
				mDefaultElementType = NO_ELEMENT_TYPE;
		}
	}

	// read element neighbors (if in file)
	if (file >> id) {	
		for (long i = 0; i < numElements; ++i) {
			pE[i]->AllocateNeighborArray();
			for (short j = 0; j < pE[i]->NumNeighbors(); ++j) {
				if ((i > 0) || (j > 0))
					file >> id;
				
				pE[i]->SetNeighbor(j, pE[id]);
			}
		}
	}
	else {
		SetElementNeighbors();
		AddExteriorBlob();
	}
	
	file.close();

	// check mesh for consistency
	CheckAll();

	return;
}



void Mesh::ReadElementPointSetData(const std::string fileName)
{
	ifstream file;
    OpenFile(fileName, file);
    
    Array<Element*> pE;
    GetElements(pE);
    
    PointSet *pPS;
    Point p;
    double x, y, z;
    for (long i = 0; i < pE.Size(); ++i) {
    	pPS = pE[i]->GetPointSet();
    	pPS->SetSize(1);
    	
    	file >> x;
    	file >> y;
    	file >> z;
    	
    	p.SetPosition(x, y, z);
    	pPS->SetPoint(0, p); 
    }
    
    file.close();
    
	return;
}



void Mesh::CheckAll() const
{
	CheckVertexElementNeighbors();
	CheckElementVertices();
	CheckListIterators();
	//CheckElementOrientation();
	CheckElementNeighbors();
	//CheckMaterialNumbers();
	
	return;
}



void Mesh::CheckVertexElementNeighbors() const
{
	list<Vertex>::const_iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) 
		(*iVl).CheckNeighbors();

	return;
}



void Mesh::CheckElementVertices() const
{
	Array<Element*> pE;
	GetElements(pE);
			
	for (long i = 0; i < pE.Size(); ++i)
		pE[i]->CheckVertices();

	return;
}



void Mesh::CheckElementOrientation() const
{
	Array<Element*> pE;
	GetElements(pE);
			
	for (long i = 0; i < pE.Size(); ++i) {
		if (pE[i]->Oriented() == false) {
			ThrowException("Mesh::CheckElementOrientation : inside out element");
		}
	}
		
	return;
}



void Mesh::CheckElementNeighbors() const
{
	Array<Element*> pE;
	GetElements(pE);
			
	for (long i = 0; i < pE.Size(); ++i) 
		pE[i]->CheckNeighbors();
	
	
	return;
}



void Mesh::CheckListIterators() const
{
	list<Vertex>::const_iterator iVl;
	for (iVl = mVertexList.begin(); iVl != mVertexList.end(); ++iVl) {
		void *pLi = (*iVl).GetListIterator();
		if (pLi == NULL)
			ThrowException("Mesh::CheckListIterators : NULL iterator");
			
		list<Vertex>::iterator iEtmp = *(list<Vertex>::iterator *) pLi;
		Vertex *pTmp = &*iEtmp;
		Vertex *pTmp1 = (Vertex *) &*iVl;
		if (pTmp != pTmp1)
			ThrowException("Mesh::CheckListIterators : bad vertex iterator");
	}
	
	list<Edge>::const_iterator iEl;
	for (iEl = mEdgeList.begin(); iEl != mEdgeList.end(); ++iEl) {
		void *pLi = (*iEl).GetListIterator();
		if (pLi == NULL)
			ThrowException("Mesh::CheckListIterators : NULL iterator");
			
		list<Edge>::iterator iEtmp = *(list<Edge>::iterator *) pLi;
		Edge *pTmp = &*iEtmp;
		Edge *pTmp1 = (Edge *) &*iEl;
		if (pTmp != pTmp1)
			ThrowException("Mesh::CheckListIterators : bad edge iterator");
	}
	
	list<Triangle>::const_iterator iTr;
	for (iTr = mTriangleList.begin(); iTr != mTriangleList.end(); ++iTr) {
		void *pLi = (*iTr).GetListIterator();
		if (pLi == NULL)
			ThrowException("Mesh::CheckListIterators : NULL iterator");
			
		list<Triangle>::iterator iEtmp = *(list<Triangle>::iterator *) pLi;
		Triangle *pTmp = &*iEtmp;
		Triangle *pTmp1 = (Triangle *) &*iTr;
		if (pTmp != pTmp1)
			ThrowException("Mesh::CheckListIterators : bad triangle iterator");
	}
	
	list<Polygon>::const_iterator iPl;
	for (iPl = mPolygonList.begin(); iPl != mPolygonList.end(); ++iPl) {
		void *pLi = (*iPl).GetListIterator();
		if (pLi == NULL)
			ThrowException("Mesh::CheckListIterators : NULL iterator");
			
		list<Polygon>::iterator iEtmp = *(list<Polygon>::iterator *) pLi;
		Polygon *pTmp = &*iEtmp;
		Polygon *pTmp1 = (Polygon *) &*iPl;
		if (pTmp != pTmp1)
			ThrowException("Mesh::CheckListIterators : bad polygon iterator");
	}
	
	list<Tetrahedron>::const_iterator iTe;
	for (iTe = mTetrahedronList.begin(); iTe != mTetrahedronList.end(); ++iTe) {
		void *pLi = (*iTe).GetListIterator();
		if (pLi == NULL)
			ThrowException("Mesh::CheckListIterators : NULL iterator");
			
		list<Tetrahedron>::iterator iEtmp = *(list<Tetrahedron>::iterator *) pLi;
		Tetrahedron *pTmp = &*iEtmp;
		Tetrahedron *pTmp1 = (Tetrahedron *) &*iTe;
		if (pTmp != pTmp1)
			ThrowException("Mesh::CheckListIterators : bad tetrahedron iterator");
	}
		
	return;
}



short Mesh::CheckMaterialNumbers() const
{
	Array<Element*> pE;
	GetElements(pE);
	
	set<short> materialNumberSet;
	
	for (long i = 0; i < pE.Size(); ++i) 
		materialNumberSet.insert(pE[i]->MaterialNumber());
		
	// material numbers must start from 0 and increment by 1
	set<short>::iterator iMs = materialNumberSet.begin();
	short m = *iMs;
	if (m != 0)
		ThrowException("Mesh::CheckMaterialNumbers : first material number is not 0, it is" + ConvertIntegerToString(m));
	
	++iMs;
	while (iMs != materialNumberSet.end()) {
		if (m != *iMs - 1)
			ThrowException("Mesh::CheckMaterialNumbers : material numbers not consecutive");
		
		m = *iMs;
		++iMs;
	}

		
	return materialNumberSet.size();
}



void Mesh::PrintVolumeByMaterialNumber() const
{
	Array<double> volume;
	VolumeByMaterialNumber(volume);
	
	double sum = 0.0;
	for (short i = 0; i < volume.Size(); ++i) {
		cout << "Volume of material " << i << " = " << volume[i] << endl;
		sum += volume[i];
	}
	
	cout << "Total volume = " << sum << endl;

	return;
}



void Mesh::PrintSurfaceArea() const
{
	if (Dimension() != 3)
		return;
	
	cout << "Surface area = " << SurfaceArea() << endl;
	
	return;
}



void Mesh::PrintLongestEdgeLength() const
{
	cout << "Longest edge length = " << LongestEdgeLength() << endl;
	return;
}



void Mesh::WriteMatlabEdgeFile()
{
	ofstream edgeFile("/Users/dave/Projects/ReDi/edges.dat");
	
	MakeEdgeList();
	
	list<Edge>::iterator iEl;
	for (iEl = mEdgeList.begin(); iEl != mEdgeList.end(); ++iEl) {
		edgeFile << (*iEl).PV(0)->X() << " " << (*iEl).PV(0)->Y() << " ";
		edgeFile << (*iEl).PV(1)->X() << " " << (*iEl).PV(1)->Y() << endl;
	}
	
	EraseEdgeList();
	edgeFile.close();

	return;
}

