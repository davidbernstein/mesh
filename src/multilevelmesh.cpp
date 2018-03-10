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
 
#include "multilevelmesh.h"
#include <iostream>

using namespace NAMESPACE;
using namespace std;

void MultiLevelMesh::RefineToConformity(RefinementQueue &refQ)
{
	while (!refQ.empty()) {
        RefinementQueue::iterator iTop = refQ.begin();
        refQ.erase(iTop);
        
        if ((*iTop)->GetBisectionData()->IsALeaf())
            RefineToConformity(*iTop);
	}
	
	return;
}



void MultiLevelMesh::RefineToConformity(Element *pE)
{
    mRefinementQueue.clear();
    mRefinementQueue.insert(pE);
        
	while (!mRefinementQueue.empty()) {
		Element *pTop = *mRefinementQueue.begin();
        mRefinementQueue.erase(mRefinementQueue.begin());
        
		if (pTop->GetBisectionData()->IsALeaf()) {
			Bisect(pTop);
		}
	}
	
	return;
}



void MultiLevelMesh::Bisect(Element *pE)
{	
	// make children
	Element *pChild[2];
	pChild[0] = AddElement(pE->Type());
	pChild[1] = AddElement(pE->Type());
	
	// this copies data associated with parent
	*pChild[0] = *pE;
	*pChild[1] = *pE;
	
	pE->InspectDataObjects();
	pChild[0]->InspectDataObjects();
	
	// find/add mid point vertex
	Vertex *pVMid = RefinementEdgeMidPoint(pE);
	
	pE->Bisect(pChild, pVMid);
	
	// put nonconforming neighbors on refinement list
	for (short i = 0; i < pE->NumNeighbors(); ++i) {
        if (pE->PN(i)->Conforming() == false) 
            mRefinementQueue.insert(pE->PN(i));
    }
    
    // put nonconforming children on refinement list
	if (pChild[0]->Conforming() == false)
		mRefinementQueue.insert(pChild[0]);
	
	if (pChild[1]->Conforming() == false)
		mRefinementQueue.insert(pChild[1]);
	
	// set bisection data of children
	UpdateBisectionData(pChild[0]);
	UpdateBisectionData(pChild[1]);
	
	
	return;
}



void MultiLevelMesh::UpdateBisectionData(Element *pChild)
{
	BisectionData *pBD = pChild->GetBisectionData();
	
	pBD->SetRefinementEdge(pChild->LongestEdge());
		 
	pBD->IncrementLevel();
	
	pBD->AllocateMidpointArray(pChild->NumEdges());
	
	Element *pParent = pBD->Parent();
	BisectionData *pBDParent = pParent->GetBisectionData();
	
	Vertex *pV[2];
	for (short i = 0; i < pChild->NumEdges(); ++i) {
		pChild->GetEdge(i, pV);
		short edgeIndex = pParent->EdgeIndex(pV[0], pV[1]);
		if (edgeIndex != -1) 
			pBD->SetMidpoint(i, pBDParent->Midpoint(edgeIndex));
	}
		
	
	return;
}



void MultiLevelMesh::UnRefineToConformity(UnRefinementQueue &unRefQ)
{
    //Display display;
    //display.ShowMesh(*this);
	while (!unRefQ.empty()) {
        Element *pTop = *unRefQ.begin();
        unRefQ.erase(unRefQ.begin());
        
        if (pTop->GetBisectionData()->IsALeaf())
            UnRefineToConformity(pTop, unRefQ);
        
        CheckAll();

        //display.ShowMesh(*this);
	}
	
	CheckAll();
    
    return;
}



void MultiLevelMesh::UnRefineToConformity(Element *pE, UnRefinementQueue &unRefQ)
{
    mUnRefinementQueue.clear();
    mUnRefinementQueue.insert(pE);
    
	while (!mUnRefinementQueue.empty()) {
		Element *pTop = *mUnRefinementQueue.begin();
        pTop->Print("top");
        PrintUnRefinementQueue();
        if (UnBisect(pTop, unRefQ))
            mUnRefinementQueue.erase(mUnRefinementQueue.begin());
    }
    
    return;
}



bool MultiLevelMesh::UnBisect(Element *pE, UnRefinementQueue &unRefQ)
{
    // is pE a leaf?
    BisectionData *pBD = pE->GetBisectionData();
    Element *pChild[2];
    if (pBD->IsALeaf() == false) {
        pBD->GetChildren(pChild);
        //BisectionData *pTmp = pChild[0]->FindBisectionData();
        mUnRefinementQueue.insert(pChild[0]);
        mUnRefinementQueue.insert(pChild[1]);
        return false;
    }
    
    // find sibling of pE
    Element *pParent = pBD->Parent();
    pParent->GetBisectionData()->GetChildren(pChild);
    Element *pSibling = (pChild[0] == pE) ? pChild[1] : pChild[0];
    
    // if sibling is not a leaf then put its children on mUnRefinementQueue
    if (pSibling->GetBisectionData()->IsALeaf() == false) {
        pSibling->GetBisectionData()->GetChildren(pChild);
        mUnRefinementQueue.insert(pChild[0]);
        mUnRefinementQueue.insert(pChild[1]);
        return false;
    }
    
    // merge pE and pSibling and determine if vertex needs to be erased
    Vertex *pVIsolated = pE->Merge(pSibling, pParent);
    if (pVIsolated != NULL)
        EraseVertex(pVIsolated);
    
    // add children of nonconforming neighbors to mUnRefinementQueue
    for (short i = 0; i < pParent->NumNeighbors(); ++i) {
        pParent->PN(i)->Print("neighbor");
        if (pParent->PN(i)->Conforming() == false) { 
            pParent->PN(i)->GetBisectionData()->GetChildren(pChild);
            pChild[0]->Print("child 0");
            pChild[1]->Print("child 1");
            mUnRefinementQueue.insert(pChild[0]);
            mUnRefinementQueue.insert(pChild[1]);
        }
    }    
    
    // remove sibling from unRefQ
    unRefQ.erase(pE);
    unRefQ.erase(pSibling);
    mUnRefinementQueue.erase(pSibling);
    
    // erase elements
    EraseElement(pE);
    EraseElement(pSibling);

        
    return true;
}



Vertex* MultiLevelMesh::RefinementEdgeMidPoint(const Element *pE)
{
	if (pE->Type() == EDGE) {
		Vertex *pVMid = AddVertex();
		pVMid->SetPosition(((Edge*) pE)->MidPoint());
		return pVMid;
	}
	
	// check if vertex is already there
	BisectionData *pBD = pE->GetBisectionData();
	Vertex *pVMid = pBD->RefinementEdgeMidPoint();
	if (pVMid != NULL)
		return pVMid;
	
	// otherwise make new vertex
	Vertex *pV[2];
	pE->GetEdge(pBD->RefinementEdge(), pV);
		
	pVMid = AddVertex();
	pVMid->SetPosition(((*pV[0]) + (*pV[1])) * 0.5);
	
    Array<Element*> ring;
    pE->Ring(ring, pBD->RefinementEdge());
    
	for (short i = 0; i < ring.Size(); ++i) 
		ring[i]->GetBisectionData()->SetMidpoint(ring[i]->EdgeIndex(pV[0], pV[1]), pVMid);
	
	return pVMid;
}



void MultiLevelMesh::SetBisectionData()
{
	Array<Element*> pE;
	GetElements(pE);
	
	for (long i = 0; i < pE.Size(); ++i) {
		BisectionData *pBD = pE[i]->GetBisectionData();
		pBD->SetRefinementEdge(pE[i]->LongestEdge());
		pBD->AllocateMidpointArray(pE[i]->NumEdges());
	}
	

	return;
}
   


void MultiLevelMesh::SetRefinementListIDs(long value)
{
    RefinementQueue::iterator iQ;
    for (iQ = mRefinementQueue.begin(); iQ != mRefinementQueue.end(); ++iQ)
        (*iQ)->SetID(value);
        
    return;
}



void MultiLevelMesh::RefineAll(short level)
{
    for (short i = 0; i < level; ++i) {
        Array<Element*> pE;
        GetElements(pE);
        
        RefinementQueue refQ;
        for (long i = 0; i < pE.Size(); ++i) {
            refQ.insert(pE[i]);
        }
    
        RefineToConformity(refQ);
    }
    
    //CheckAll();
    
    //MeshDisplay display;
    //display.ShowMesh(*this);
    
    return;
}



void MultiLevelMesh::TestRefinement()
{
	for (short j = 0; j < 2; ++j) {
		Array<Element*> pE;
		GetElements(pE);
	
		RefinementQueue refQ;
		Point c;
		for (long i = 0; i < pE.Size(); ++i) {
			pE[i]->Centroid(c);
			if (c.X() < 0.5)
				refQ.insert(pE[i]);
		}
        
		RefineToConformity(refQ);
	}
	
	
	return;
}



void MultiLevelMesh::TestUnRefinement()
{
    for (short j = 0; j < 8; ++j) {
		Array<Element*> pE;
		GetElements(pE);
        
		UnRefinementQueue unRefQ;
		Point c;
		for (long i = 0; i < pE.Size(); ++i) {
			pE[i]->Centroid(c);
			if ((c.X() > 2.0) && (c.X() < 3.0))
				unRefQ.insert(pE[i]);
		}
        
		UnRefineToConformity(unRefQ);
	}
    
    return;
}



void MultiLevelMesh::PrintRefinementQueue(short numTop) const
{
	long numPrinted = (numTop < 1) ? mRefinementQueue.size() : numTop;
	
	long i = 0;
	RefinementQueue::const_iterator iQ = mRefinementQueue.begin();
	while ((iQ != mRefinementQueue.end()) && (i < numPrinted)) {
		BisectionData *pBD = (*iQ)->GetBisectionData();
		cout << pBD->Level() << " " 
            << endl;
		
		++iQ;
		++i;
	}
	
	return;
}



void MultiLevelMesh::PrintUnRefinementQueue(short numTop) const
{
	long numPrinted = (numTop < 1) ? mUnRefinementQueue.size() : numTop;
	
	long i = 0;
	UnRefinementQueue::const_iterator iQ = mUnRefinementQueue.begin();
	while ((iQ != mUnRefinementQueue.end()) && (i < numPrinted)) {
        (*iQ)->Print();
		//BisectionData *pBD = (*iQ)->FindBisectionData();
		//cout << pBD->Level() << " "  << endl;
		
		++iQ;
		++i;
	}
	
	return;
}
