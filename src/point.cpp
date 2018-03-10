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
 
#include <iostream>
#include <math.h>

#include "constants.h"
#include "point.h"
#include "utility.h"

using namespace NAMESPACE;
using namespace std;


double Point::Angle(const Point *pPt) const
{
	// returns the angle between this and pPt
	// 0 <= angle <= PI

	double normThis = Norm();
	double normPt = pPt->Norm();

	if ((normThis == 0.0) || (normPt == 0.0)) 
		ThrowException("Point::Angle : zero norm");

	return acos(this->Dot(*pPt)) / (normThis * normPt);
}



double Point::XYAngle() const
{
	// angle this point makes with origin of coordinate system
	// measured counterclockwise from the x axis
	
	if ((mX == 0.0) && (mY == 0.0))
		return 0.0;

	double hfpi = 0.5 * PI;
	double angle;

	if (mX != 0.0) {
		angle = atan(mY / mX);
	}
	else {
		return PI + (mY < 0.0 ? 1.0 : -1.0) * hfpi;
	}

	if (mX < 0.0)
		angle += PI;

	if (angle < 0.0)
		angle += 2.0 * PI;

	return angle;
}



bool Point::SameLocation(const Point &p) const
{
	if (FloatComparison(mX, p.mX, FLOAT_EQUALITY, FLOAT_COMPARE_TOL) == false)
		return false;
			
	if (FloatComparison(mY, p.mY, FLOAT_EQUALITY, FLOAT_COMPARE_TOL) == false)
		return false;
			
	return FloatComparison(mZ, p.mZ, FLOAT_EQUALITY, FLOAT_COMPARE_TOL);
}



void Point::Print() const
{
	cout << "this = " << this << endl;
	cout << "ID = " << mID << endl;
	cout << "location = (" << mX << ", " << mY << ", " << mZ << ")" << endl;

	return;
}



bool Point::operator==(const Point &p) const
{
	ThrowException("Point::operator== called");
	return false;
}



bool Point::operator!=(const Point &p) const
{
	ThrowException("Point::operator!= called");
	return false;
}



bool Point::operator<(const Point &p) const
{
	ThrowException("Point::operator< called");
	return false;
}



bool Point::operator>(const Point &p) const
{
	ThrowException("Point::operator> called");
	return false;
}
