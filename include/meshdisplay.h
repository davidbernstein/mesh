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
 
#ifndef _meshdisplay_h_
#define _meshdisplay_h_

#include "multilevelmesh.h"
#include "xydata.h"

#include "vtkUnstructuredGrid.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkHexahedron.h"
#include "vtkProperty.h"
#include "vtkDataSetMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"

#include <string>

namespace NAMESPACE {
	class MeshDisplay {
    public:
		// Constructors
		MeshDisplay(MeshDisplayType displayType = DISPLAY_DEFAULT, short component = 0);
        
		// Destructor
        ~MeshDisplay(void);
        
        // functions for displaying meshes
        void ShowMesh(const Mesh &mesh);
        void SetDisplayType(MeshDisplayType displayType);
        void SetComponent(short component);
        
        // IO
        void SaveCurrentImageAsJPEG(const std::string fileName) const;
    	void SaveCurrentImageAsPNG(const std::string fileName) const;
    	
    private:
        void LoadMesh(const Mesh *pMesh);
        void SetElementColors(const Mesh &mesh);
        float* GetDefaultColors(const Mesh &mesh) const;
        float* GetConcentrationColors(const Mesh &mesh) const;
        float* GetOccupancyColors(const Mesh &mesh) const;
        float* GetMaterialColors(const Mesh &mesh) const;
        
    private:
    	// data type
    	MeshDisplayType mDisplayType;
    	short mComponent;
    	
    	// vtk data structures
        vtkUnstructuredGrid *mpVtkGrid;
        vtkPoints *mpVtkPoints;  
        vtkTriangle *mpVtkTriangle;
        vtkTetra *mpVtkTetra;
        vtkHexahedron *mpVtkHexahedron;
        
        vtkDataSetMapper *mpVtkGridMapper;
        vtkDataSetMapper *mpVtkWireMapper;

        vtkActor *mpVtkGridActor;
        vtkActor *mpVtkWireActor;
        
        vtkRenderWindow *mpVtkRenderWindow;
        vtkRenderer *mpVtkRenderer;
        vtkRenderWindowInteractor *mpVtkInteractor;
	};
	
	
	
	inline void MeshDisplay::SetDisplayType(MeshDisplayType displayType)
	{
		mDisplayType = displayType;
		return;
	}
	
	
	
    inline void MeshDisplay::SetComponent(short component)
    {
    	mComponent = component;
    	return;
    }
}


#endif // _meshdisplay_h_



