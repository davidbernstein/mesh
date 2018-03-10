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
 
#ifndef _occupancy_h_
#define _occupancy_h_

#include "dataobject.h"
#include "array.h"
#include "utility.h"

namespace NAMESPACE {
	class Occupancy : public DataObject {
	public:
		Occupancy(void) { };
		
		long Size(void) const;
		DataType Type(void) const {return OCCUPANCY;}
		
		void SetSize(short arraySize);
		bool Occupied(void) const;
		short NumOccupants(short index = -1) const;
		void IncrementNumOccupants(short index, short change = 1);
		void SetNumOccupants(short index, short numOccupants);

	private:
		Array<short> mData;
	};
	
	
	
	inline void Occupancy::SetSize(short arraySize)
	{
		mData.SetSize(arraySize);
		return;
	}
	
	
	
	inline long Occupancy::Size() const
	{
		return mData.Size();
	}
	
	

	inline short Occupancy::NumOccupants(short index) const
	{
		return mData[index];
	}
	
	
	
	inline void Occupancy::SetNumOccupants(short index, short numOccupants)
	{
		if (numOccupants < 0)
			ThrowException("Occupancy::SetNumOccupants : numOccupants is negative");
			
		mData[index] = numOccupants;
		return;
	}
	


	inline void Occupancy::IncrementNumOccupants(short index, short change)
	{
		mData[index] += change;

		if (mData[index] < 0)
			ThrowException("Occupancy::IncrementNumOccupants : negative occupancy");
			
		return;
	}
}


#endif // _occupancy_h_

