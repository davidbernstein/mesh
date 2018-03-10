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
 
#ifndef _dataobject_h_
#define _dataobject_h_

#include "meshenums.h"

namespace NAMESPACE {
	class DataObject {
	public:
		// constructors
		DataObject(void) { };

		// destructor
		virtual ~DataObject(void) { };
		
		// interface
		virtual DataType Type(void) const = 0;
		virtual void SetSize(short arraySize) = 0;
		virtual long Size(void) const = 0;
	};
}


#endif // _dataobject_h_
