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
 
#ifndef _boundarydata_h_
#define _boundarydata_h_

#include "dataobject.h"
#include "meshenums.h"

namespace NAMESPACE {
	class BoundaryData : public DataObject {
    public:
		BoundaryData(void);
		~BoundaryData(void) { };
		
		DataType Type(void) const {return BOUNDARY_DATA;};
        
    private:
        // type
        BoundaryType mType;
        
        // data
        Array<double> mpData;
	};
	
	
	
	inline BoundaryData::BoundaryData()
	{
		mType = NO_BOUNDARY_TYPE;
		return;
	}
}


#endif // _boundarydata_h_



