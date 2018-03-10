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
 
#ifndef _mesh_h_
#define _mesh_h_

#include "vertex.h"
#include "edge.h"
#include "triangle.h"
#include "polygon.h"
#include "tetrahedron.h"
#include "blob.h"
#include "brick.h"
#include "matrix.h"

namespace NAMESPACE {
	class Mesh
	{
    public:
        Mesh(void);
		virtual ~Mesh(void) { };
        
        // vertices
        long NumVertices(void) const;
        void SaveVertexIDs(Array<long> &id) const;
        void SetVertexIDs(void);
        void SetVertexIDs(short value);
        void GetVertices(Array<Vertex*> &pV);
  		void RestoreVertexIDs(Array<long> &id);
  		void SetVertexElementNeighbors(ElementType eType = DEFAULT_ELEMENT_TYPE);
		void EraseVertexElementNeighbors(void);
		Vertex* IsAVertex(const Point &p) const;

        // elements
        long NumElements(ElementType eType = DEFAULT_ELEMENT_TYPE, bool leavesOnly = true) const;
        long NumElementsTotal(ElementType eType = DEFAULT_ELEMENT_TYPE) const;
        void GetElements(Array<Element*> &pE, ElementType eType = DEFAULT_ELEMENT_TYPE, bool leavesOnly = true) const;
        ElementType DefaultElementType(void) const;

        // geometry
        double Volume(ElementType eType = DEFAULT_ELEMENT_TYPE) const;
        void VolumeByMaterialNumber(Array<double> &volume, bool normalizeByTotal = false) const;
        double SurfaceArea(void) const;
		short Dimension(void) const;
		double LongestEdgeLength(void) const;
        
        // internal file format reader and writer
		void WriteToFile(std::string fileName) const;
		void ReadFromFile(std::string fileName);
        
        // consistency checking routines
		void CheckAll(void) const;
 		short CheckMaterialNumbers(void) const;
       		
		// element data
		void SetElementDataSize(short numData, ElementType elementType = DEFAULT_ELEMENT_TYPE);
		void AddElementData(DataType dataType, short arraySize = -1, ElementType elementType = DEFAULT_ELEMENT_TYPE);
		
		// material number
		void MaterialNumberTest(void);
		
		// other IO
		void WriteMatlabEdgeFile(void);
		void PrintVolumeByMaterialNumber(void) const;
		void PrintSurfaceArea(void) const;
		void PrintLongestEdgeLength(void) const;
		void ReadElementPointSetData(const std::string fileName);
		
	protected:
		void Erase(void);

		// copy constructor
		Mesh(const Mesh &mesh);

    	// operators
    	Mesh& operator=(const Mesh &mesh);

		// vertices
    	Vertex* AddVertex(const Vertex &v);
    	Vertex* AddVertex(void);
		void EraseVertex(Vertex *pV);
		void CorrectVertexNeighbors(Array<Element*> &pE);
        
		// elements
		Element* AddElement(ElementType type);
		Element* AddElement(const Element *pE);
		void EraseElement(Element *pE);
		void SetElementIDs(ElementType eType = ALL_ELEMENT_TYPES);
		void UpdateVertexNeighbors(Array<Element*> &pE);
		void CorrectElementVertices(const Array<Vertex*> &pV);
		void OrientElements(ElementType elementType);
		Element* FindElementByRegion(long regionNumber, ElementType eType) const;
		void SetElementNeighbors(void);
		void RestoreElementIDs(Array<long> &id);
		void SaveElementIDs(Array<long> &id) const;
		void SetElementIDs(short value, ElementType eType = ALL_ELEMENT_TYPES);
		
		// elementn neighbors
		void AddExteriorBlob(ElementType eType = DEFAULT_ELEMENT_TYPE, short materialNumber = 0);

		// edges
		void MakeEdgeList(void);
		void EraseEdgeList(void);	
		
		// geometrical functions 
		void Translate(const Point &trans);
		void Scale(double scaleFactor);
		void BoundingBox(Point &min, Point &max) const;
		
		// material number
		void SetMaterialNumber(short m, ElementType elementType = DEFAULT_ELEMENT_TYPE);
		
		// IO 		
        // MC
		void WriteMCHeader(std::ofstream &file) const;
		void WriteMCVertices(std::ofstream &file) const;
		void WriteMCElements(std::ofstream &file, Array<short> &elementGroupNumber, 
							 Array<Array<short> > &faceType) const;
		void WriteMCFooter(std::ofstream &file) const;
		void WriteMCFile(const std::string &fileName, Array<short> &elementGroupNumber, Array<Array<short> > &faceType) const;
		void WriteMCFile(std::ofstream &file, Array<short> &elementGroupNumber, Array<Array<short> > &faceType) const;
	
	protected:
		void CheckVertexElementNeighbors(void) const;
		void CheckElementVertices(void) const;
		void CheckListIterators(void) const;
		void CheckElementOrientation(void) const;
		void CheckElementNeighbors(void) const;
        
	protected:
		// list of vertices 
		std::list<Vertex> mVertexList;
		
		// list of edges
		std::list<Edge> mEdgeList;
		
		// list of triangles
		std::list<Triangle> mTriangleList;
		
		// list of polygones
		std::list<Polygon> mPolygonList;
		
		// list of tetrahedra
		std::list<Tetrahedron> mTetrahedronList;
		
		// list of bricks	
		std::list<Brick> mBrickList;
		
		// list of blobs
		std::list<Blob> mBlobList;
            
		// primary element type
		ElementType mDefaultElementType;
    };
	
	
	
	inline Mesh::Mesh() 
	{
		mDefaultElementType = NO_ELEMENT_TYPE;

		return;
	}
		
	

	inline long Mesh::NumVertices() const
	{
		return mVertexList.size();
	}
	
	

	inline short Mesh::Dimension() const
	{
		switch (mDefaultElementType) {
		case EDGE:
			return 1;
			break;
	
		case TRIANGLE:
			return 2;
			break;
			
		case POLYGON:
			return 2;
			break;
	
		case TETRAHEDRON:
			return 3;
			break;
		
		case BRICK:
			return 3;
			break;
			
		default:
			ThrowException("Mesh::Dimension : no default element type");
			return -1;
			break;
		}
	}
	
	
	
	inline ElementType Mesh::DefaultElementType() const
	{
		return mDefaultElementType;
	}
}


#endif // _mesh_h_	
