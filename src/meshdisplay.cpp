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
 
#include "meshdisplay.h"

#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPolyDataMapper.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"
#include "vtkWindowToImageFilter.h"

using namespace NAMESPACE;
using namespace std;

MeshDisplay::MeshDisplay(MeshDisplayType displayType, short component)
{
	mDisplayType = displayType;
    mComponent = component;
    	
    mpVtkGrid = vtkUnstructuredGrid::New();
    mpVtkPoints = vtkPoints::New();
    mpVtkTriangle = vtkTriangle::New();
    mpVtkTetra = vtkTetra::New();
    mpVtkHexahedron = vtkHexahedron::New();
    
    mpVtkGridMapper = vtkDataSetMapper::New();
    mpVtkWireMapper = vtkDataSetMapper::New();
    
    mpVtkGridActor = vtkActor::New();
    mpVtkWireActor = vtkActor::New();
    
    mpVtkRenderWindow = vtkRenderWindow::New();
    mpVtkRenderer = vtkRenderer::New();
    mpVtkInteractor = vtkRenderWindowInteractor::New();
    
    // set up
    mpVtkGridMapper->SetInput(mpVtkGrid);
    mpVtkWireMapper->SetInput(mpVtkGrid);
    mpVtkWireMapper->ScalarVisibilityOff();
    
    mpVtkGridActor->SetMapper(mpVtkGridMapper);
    mpVtkGridActor->GetProperty()->SetRepresentationToSurface();
    //mpVtkGridActor->GetProperty()->SetOpacity(0.7);
    
    mpVtkWireActor->SetMapper(mpVtkWireMapper);
    mpVtkWireActor->GetProperty()->SetRepresentationToWireframe();
    mpVtkWireActor->GetProperty()->SetColor(0, 0, 0);
        
    mpVtkGridMapper->SetResolveCoincidentTopologyPolygonOffsetParameters(0,1);
	mpVtkGridMapper->SetResolveCoincidentTopologyToPolygonOffset();
	mpVtkWireMapper->SetResolveCoincidentTopologyPolygonOffsetParameters(1,1);
	mpVtkWireMapper->SetResolveCoincidentTopologyToPolygonOffset();

    mpVtkRenderer->AddActor(mpVtkGridActor);
    
    // comment out this line to turn off edges
    mpVtkRenderer->AddActor(mpVtkWireActor);
    
    //mpVtkRenderer->SetBackground(0.95, 0.95, 0.95);
    mpVtkRenderer->SetBackground(1.0, 1.0, 1.0);
    
    mpVtkRenderWindow->AddRenderer(mpVtkRenderer);
    mpVtkRenderWindow->SetSize(500, 500);
    mpVtkRenderWindow->SetPosition(00, 50);
    
    // add interactor
    mpVtkInteractor->SetRenderWindow(mpVtkRenderWindow);
    
    
    return;
}



MeshDisplay::~MeshDisplay()
{
    mpVtkGrid->Delete();
    mpVtkPoints->Delete();
    mpVtkTriangle->Delete();
    mpVtkTetra->Delete();
    mpVtkHexahedron->Delete();
    
    mpVtkGridMapper->Delete();
    mpVtkWireMapper->Delete();
    
    mpVtkGridActor->Delete();
    mpVtkWireActor->Delete();
    
    mpVtkRenderer->Delete();
    mpVtkInteractor->Delete();
    mpVtkRenderWindow->Delete();
    
    return;
}



void MeshDisplay::ShowMesh(const Mesh &mesh)
{
    LoadMesh(&mesh);
    SetElementColors(mesh);
    
    vtkCamera *pVtkCamera = mpVtkRenderer->GetActiveCamera();
    pVtkCamera->Azimuth(30);  
    pVtkCamera->Elevation(30);
    pVtkCamera->Zoom(1.05);
    
    //SaveCurrentImageAsJPEG("/Users/dave/Projects/ReDi/mesh.jpg");
    SaveCurrentImageAsPNG("/Users/dave/Projects/ReDi/mesh.png");
    
    // render window and start interactor
    mpVtkRenderWindow->Render();
    mpVtkInteractor->Start();
    
    
	return;
}



void MeshDisplay::LoadMesh(const Mesh *pMesh)
{
    mpVtkGrid->Initialize();
    
    Mesh *pNonConstMesh = (Mesh*) pMesh;
    
    Array<long> vID;
    pNonConstMesh->SaveVertexIDs(vID);
    
    mpVtkPoints->SetNumberOfPoints(pMesh->NumVertices());
    list<Vertex>::iterator iVl;
    long i = 0;
    pNonConstMesh->SetVertexIDs();
    Array<Vertex*> pV;
    pNonConstMesh->GetVertices(pV);
    for (long i = 0; i < pV.Size(); ++i) 
        mpVtkPoints->InsertPoint(i, pV[i]->X(), pV[i]->Y(), pV[i]->Z());
    
    mpVtkGrid->SetPoints(mpVtkPoints);
    
    Array<Element*> pE;
    pNonConstMesh->GetElements(pE);
    mpVtkGrid->Allocate(pE.Size(), 1);
    for (i = 0; i < pE.Size(); ++i) {
        switch (pE[i]->Type()) {
            case TRIANGLE:
                mpVtkTriangle->GetPointIds()->SetId(0, pE[i]->PV(0)->ID());
                mpVtkTriangle->GetPointIds()->SetId(1, pE[i]->PV(1)->ID());
                mpVtkTriangle->GetPointIds()->SetId(2, pE[i]->PV(2)->ID());
                mpVtkGrid->InsertNextCell(mpVtkTriangle->GetCellType(), mpVtkTriangle->GetPointIds());
                break;
                
            case TETRAHEDRON:
                mpVtkTetra->GetPointIds()->SetId(0, pE[i]->PV(0)->ID());
                mpVtkTetra->GetPointIds()->SetId(1, pE[i]->PV(1)->ID());
                mpVtkTetra->GetPointIds()->SetId(2, pE[i]->PV(2)->ID());
                mpVtkTetra->GetPointIds()->SetId(3, pE[i]->PV(3)->ID());
                mpVtkGrid->InsertNextCell(mpVtkTetra->GetCellType(), mpVtkTetra->GetPointIds());
                break;
                   
            case BRICK:
                mpVtkHexahedron->GetPointIds()->SetId(0, pE[i]->PV(0)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(1, pE[i]->PV(1)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(2, pE[i]->PV(2)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(3, pE[i]->PV(3)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(4, pE[i]->PV(4)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(5, pE[i]->PV(5)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(6, pE[i]->PV(6)->ID());
                mpVtkHexahedron->GetPointIds()->SetId(7, pE[i]->PV(7)->ID());
                mpVtkGrid->InsertNextCell(mpVtkHexahedron->GetCellType(), mpVtkHexahedron->GetPointIds());
                break;
                
            default:
                ThrowException("Mesh::Display : element type not supported");
                break;
        }
    }
        
    pNonConstMesh->RestoreVertexIDs(vID);
    
    
    return;
}



void MeshDisplay::SetElementColors(const Mesh &mesh)
{        
	long numElements = mesh.NumElements();
    
    vtkFloatArray *pVtkFloatArray = vtkFloatArray::New();
    float *pScalar = NULL;
    
    switch (mDisplayType) {
    case DISPLAY_DEFAULT:
    	pScalar = GetDefaultColors(mesh);
    	break;
    	
    case DISPLAY_CONCENTRATION:
    	pScalar = GetConcentrationColors(mesh);
    	break;
    	
    case DISPLAY_OCCUPANCY:
    	pScalar = GetOccupancyColors(mesh);
    	break;
    	
    case DISPLAY_MATERIAL:
    	pScalar = GetMaterialColors(mesh);
    	break;
    	
    default:
    	ThrowException("Display::SetElementColors : unsupported display type");
    	break;
    }
    
    pVtkFloatArray->SetArray(pScalar, numElements, 0);
    mpVtkGrid->GetCellData()->SetScalars(pVtkFloatArray);
    
    pVtkFloatArray->Delete();
    //delete [] pScalar;
    
    return;
}



float* MeshDisplay::GetDefaultColors(const Mesh &mesh) const
{
	long numElements = mesh.NumElements();
	
    float *pScalar = new float[numElements];
    
    for (long i = 0; i < numElements; ++i) 
        pScalar[i] = 0.8;
    
	return pScalar;
}



float* MeshDisplay::GetConcentrationColors(const Mesh &mesh) const
{
	Array<Element*> pE;
    mesh.GetElements(pE);
    
    float *pScalar = new float[pE.Size()];
    
    float maxConcentration = -1.0;
    for (long i = 0; i < pE.Size(); ++i) {
    	Occupancy *pOcc = pE[i]->GetOccupancyData();
    	pScalar[i] = pOcc->NumOccupants(mComponent) / pE[i]->Volume();
    	
    	maxConcentration = max(pScalar[i], maxConcentration);
    }
    
    for (long i = 0; i < pE.Size(); ++i) {
    	pScalar[i] /= maxConcentration;
    }
    
	return pScalar;
}



float* MeshDisplay::GetOccupancyColors(const Mesh &mesh) const
{
	Array<Element*> pE;
    mesh.GetElements(pE);
    
    float *pScalar = new float[pE.Size()];
    
    float maxCount = -1.0;
    for (long i = 0; i < pE.Size(); ++i) {
    	Occupancy *pOcc = pE[i]->GetOccupancyData();
    	pScalar[i] = (float) pOcc->NumOccupants(mComponent);
    	
    	maxCount = max(pScalar[i], maxCount);
    }
    
    for (long i = 0; i < pE.Size(); ++i) {
    	pScalar[i] /= maxCount;
    }
    
	return pScalar;
}



float* MeshDisplay::GetMaterialColors(const Mesh &mesh) const
{
	Array<Element*> pE;
    mesh.GetElements(pE);
    
    float *pScalar = new float[pE.Size()];
    
    float maxCount = -1.0;
    for (long i = 0; i < pE.Size(); ++i) {
    	if (pE[i]->MaterialNumber() == 0)
    		pScalar[i] = 0.8;
    	
    	if (pE[i]->MaterialNumber() == 1)
    		pScalar[i] = 0.2;
    }
   
    
	return pScalar;
}



void MeshDisplay::SaveCurrentImageAsJPEG(const string fileName) const
{
    vtkJPEGWriter *pVtkJPEGWriter = vtkJPEGWriter::New();
    vtkWindowToImageFilter *pVtkWindowToImageFilter = vtkWindowToImageFilter::New();
    
    pVtkWindowToImageFilter->SetInput(mpVtkRenderWindow);
    
    pVtkJPEGWriter->SetInput(pVtkWindowToImageFilter->GetOutput());
    pVtkJPEGWriter->SetFileName(fileName.c_str());
    pVtkJPEGWriter->Write();
    
    pVtkJPEGWriter->Delete();
    pVtkWindowToImageFilter->Delete();
    
    return;
}



void MeshDisplay::SaveCurrentImageAsPNG(const string fileName) const
{
    vtkPNGWriter *pVtkPNGWriter = vtkPNGWriter::New();
    vtkWindowToImageFilter *pVtkWindowToImageFilter = vtkWindowToImageFilter::New();
    
    pVtkWindowToImageFilter->SetInput(mpVtkRenderWindow);
    
    pVtkPNGWriter->SetInput(pVtkWindowToImageFilter->GetOutput());
    pVtkPNGWriter->SetFileName(fileName.c_str());
    pVtkPNGWriter->Write();
    
    pVtkPNGWriter->Delete();
    pVtkWindowToImageFilter->Delete();
    
    return;
}
