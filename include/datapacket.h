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
 
#ifndef _datapacket_h_
#define _datapacket_h_

#include "dataobject.h"
#include "occupancy.h"
#include "bisectiondata.h"
#include "diffusionfactor.h"
#include "pointset.h"

#include "array.h"

namespace NAMESPACE {
	class DataPacket {
	public:
		// Constructor
		DataPacket(void) { };
		DataPacket(const DataPacket &dp);
		DataPacket& operator=(const DataPacket &dp);

		// Destructor
		~DataPacket(void);
		void Erase(void);
		
		// size
		void SetDataSize(long numData);
		long DataSize(void) const;
		
		// data types
		void AddDataObject(DataType type, short arraySize);
		DataObject* GetDataObject(DataType type, bool throwException = true) const;
		
		// find data 
		Occupancy* GetOccupancyData(bool throwException = true) const;
		BisectionData* GetBisectionData(bool throwException = true) const;
		DiffusionFactor* GetDiffusionFactor(bool throwException = true) const;
		PointSet* GetPointSet(bool throwException = true) const;

		// checking
		void InspectDataObjects(void) const;
		
	private:
		Array<DataObject*> mpDataObject;
	};
	
	
	
	inline void DataPacket::SetDataSize(long numData)
	{
        Erase();
		mpDataObject.SetSize(numData);
		
		for (short i = 0; i < numData; ++i)
			mpDataObject[i] = NULL;
		
		return;
	}
	
	
	
	inline long DataPacket::DataSize() const
	{
		return mpDataObject.Size();
	}
	


	inline DataObject* DataPacket::GetDataObject(DataType type, bool throwException) const
	{
		for (short i = 0; i < mpDataObject.Size(); ++i) {
			if (mpDataObject[i] != NULL) {
				if (mpDataObject[i]->Type() == type)
					return mpDataObject[i];
			}
		}
	
		if (throwException)	{
			InspectDataObjects();
			ThrowException("DataPacket::FindDataObject : no data");
		}
	
		return NULL;
	}



	inline Occupancy* DataPacket::GetOccupancyData(bool throwException) const
	{
		return (Occupancy*) GetDataObject(OCCUPANCY, throwException);
	}
	
	
	
	inline BisectionData* DataPacket::GetBisectionData(bool throwException) const
	{
		return (BisectionData*) GetDataObject(BISECTION_DATA, throwException);
	}
	
	
	
	inline DiffusionFactor* DataPacket::GetDiffusionFactor(bool throwException) const
	{
		return (DiffusionFactor*) GetDataObject(DIFFUSION_FACTOR, throwException);
	}
	
	
	
	inline PointSet* DataPacket::GetPointSet(bool throwException) const
	{
		return (PointSet*) GetDataObject(POINT_SET, throwException);
	}
}


#endif // _datapacket_h_	
