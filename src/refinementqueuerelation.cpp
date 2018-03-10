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
 
#include "meshtypes.h"
#include "element.h"

using namespace NAMESPACE;
using namespace std;

bool RefinementQueueRelation::operator()(const Element *pE0, const Element *pE1) const
{
	short level0 = pE0->GetBisectionData()->Level();
	short level1 = pE1->GetBisectionData()->Level();
	
	if (level0 != level1) 
		return level0 < level1;
	else 
		return pE0 < pE1;
}



bool UnRefinementQueueRelation::operator()(const Element *pE0, const Element *pE1) const
{
	short level0 = pE0->GetBisectionData()->Level();
	short level1 = pE1->GetBisectionData()->Level();
	
	if (level0 != level1) 
		return level0 > level1;
	else 
		return pE0 > pE1;
}
