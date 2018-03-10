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
 
#ifndef _meshenums_h_
#define _meshenums_h_

#include "namespace.h"

namespace NAMESPACE {
	// new element types should be added on to the end of this enum 
	// otherwise existing mesh files will be invalid
	enum ElementType{NO_ELEMENT_TYPE, DEFAULT_ELEMENT_TYPE, ALL_ELEMENT_TYPES, 
                     EDGE, TRIANGLE, TETRAHEDRON, BLOB, POLYGON, BRICK};
	
	enum DataType{OCCUPANCY, BISECTION_DATA, DIFFUSION_FACTOR, POINT_SET};
        
    enum BoundaryType{NO_BOUNDARY_TYPE, SPHERE};
        
    enum MeshDisplayType{DISPLAY_DEFAULT, DISPLAY_CONCENTRATION, DISPLAY_OCCUPANCY, DISPLAY_MATERIAL};
}

#endif // _meshenums_h_
