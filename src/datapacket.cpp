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
 
#include "datapacket.h"

using namespace NAMESPACE;
using namespace std;


// Destructor
DataPacket::~DataPacket() 
{ 
	Erase();
	return;
}


// copy constructor
DataPacket::DataPacket(const DataPacket &dp)
{
	*this = dp;
	return;
}



DataPacket& DataPacket::operator=(const DataPacket &dp)
{
	Erase();
	
	SetDataSize(dp.DataSize());
	
    BisectionData *pBD;
    Occupancy *pOccupancy;
    DiffusionFactor *pDF;
    PointSet *pPS;
    
	for (short i = 0; i < dp.DataSize(); ++i) {
		if (dp.mpDataObject[i] != NULL) {
			switch (dp.mpDataObject[i]->Type()) {
			case OCCUPANCY:
                pOccupancy = new Occupancy;
				*pOccupancy = *((Occupancy*) dp.mpDataObject[i]);
                mpDataObject[i] = pOccupancy;
				break;
				
			case BISECTION_DATA:
                pBD = new BisectionData;
				*pBD = *((BisectionData*) dp.mpDataObject[i]);
                mpDataObject[i] = pBD;
				break;
					
			case DIFFUSION_FACTOR:
                pDF = new DiffusionFactor;
				*pDF = *((DiffusionFactor*) dp.mpDataObject[i]);
                mpDataObject[i] = pDF;
				break;
					
			case POINT_SET:
                pPS = new PointSet;
				*pPS = *((PointSet*) dp.mpDataObject[i]);
                mpDataObject[i] = pPS;
				break;
				
			default:
				ThrowException("DataPacket::~DataPacket : bad data type");
				break;
			}
		}
	}
	
	return *this;
}



void DataPacket::Erase()
{
	for (short i = 0; i < mpDataObject.Size(); ++i) {
		if (mpDataObject[i] != NULL) {
			switch (mpDataObject[i]->Type()) {
			case OCCUPANCY:
				delete ((Occupancy*) mpDataObject[i]);
				break;
			
			case BISECTION_DATA:
				delete ((BisectionData*) mpDataObject[i]);
				break;
		
			case DIFFUSION_FACTOR:
				delete ((DiffusionFactor*) mpDataObject[i]);
				break;
				
			case POINT_SET:
				delete ((PointSet*) mpDataObject[i]);
				break;
				
			default:
				ThrowException("DataPacket::~DataPacket : bad data type");
				break;
			}
		}
	}
	
	mpDataObject.Erase();
	
	return;
}
	


void DataPacket::AddDataObject(DataType type, short arraySize)
{
	short i = 0;
	short index = -1;
	while ((i < mpDataObject.Size()) && (index == -1)) {
		if (mpDataObject[i] == NULL) 
			index = i;
		++i;
	}
		
	if (index == -1)
		ThrowException("DataPacket::AddDataObject : no space in array");
		
	switch (type) {
	case OCCUPANCY:
		mpDataObject[index] = new Occupancy;
		mpDataObject[index]->SetSize(arraySize);
		break;
	
	case BISECTION_DATA:
		mpDataObject[index] = new BisectionData;
		break;

	case DIFFUSION_FACTOR:
		mpDataObject[index] = new DiffusionFactor;
		break;
		
	case POINT_SET:
		mpDataObject[index] = new PointSet;
		break;
		
	default:
		ThrowException("DataPacket::AddDataObject : bad data type");
		break;
	}
	
	
	return;
}



void DataPacket::InspectDataObjects() const
{
	for (short i = 0; i < mpDataObject.Size(); ++i) {
		//DataObject *pDO = mpDataObject[i];
		//DataType type = mpDataObject[i]->Type();
	}
	
	return;
}