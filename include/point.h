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
 
#ifndef _point_h_
#define _point_h_

#include "utility.h"
#include "constants.h"
#include <math.h>

namespace NAMESPACE {
	class Point
	{
	public:
		// Constructor
		Point(double x = 0.0, double y = 0.0, double z = 0.0);

		// Destructor
		virtual ~Point(void) { };

		// Sets and Gets
		void SetPosition(double x = 0.0, double y = 0.0, double z = 0.0);
		void SetPosition(const Point &p);
		double X(void) const;
		double Y(void) const;
		double Z(void) const;
		void SetID(long id);
		long ID(void) const;

		// operators
		Point& operator=(const Point &p);
		Point operator+(const Point &p) const;
		Point operator-(const Point &p) const;
		Point operator-(void);
		Point& operator+=(const Point &p);
		Point& operator-=(const Point &p);
		Point& operator*=(double a);
		Point operator*(double a) const;
		Point operator/(double a) const;
		Point& operator/=(double a);
	
		// geometry
		double SqrNorm(void) const;
		double Norm(void) const;
		void Normalize(void);
		double Dot(const Point &p) const;
		Point Cross(const Point &p) const;
		bool SameLocation(const Point &p) const;
		double SqrDistanceToPoint(const Point *pPt) const;
		double DistanceToPoint(const Point *pPt) const;
		double DistanceToPoint(const Point &p) const;
		double Angle(const Point *pPt) const;
		double XYAngle(void) const;
		double XYCross(const Point *pPt) const;
		
		// misc
		void Min(Point &pMin) const;
		void Max(Point &pMax) const;
		void Print(void) const;

		// operators needed by STL
		bool operator==(const Point &p) const;
		bool operator!=(const Point &p) const;
		bool operator<(const Point &p) const;
		bool operator>(const Point &p) const;

	protected:
		// coordinates
		double mX;
		double mY;
		double mZ;
		
		// id number
		long mID;
	};
	
	
	
	inline Point::Point(double x, double y, double z)
	{
		mX = x;
		mY = y;
		mZ = z;
		mID = -1;

		return;
	}



	inline void Point::SetPosition(double x, double y, double z)
	{
		mX = x;
		mY = y;
		mZ = z;
		
		return;
	}



	inline void Point::SetPosition(const Point &p)
	{
		mX = p.mX;
		mY = p.mY;
		mZ = p.mZ;
			
		return;
	}



	inline Point& Point::operator=(const Point &p)
	{
		// check for assignment to self
		if (this == &p)
			return *this;

		SetPosition(p);
		mID = p.mID;

		return *this;
	}
	


	inline double Point::X() const
	{
		return mX;
	}



	inline double Point::Y() const
	{
		return mY;
	}



	inline double Point::Z() const
	{
		return mZ;
	}
	


	inline void Point::SetID(long id)
	{
		mID = id;
		return;
	}



	inline long Point::ID() const
	{
		return mID;
	}



	inline double Point::SqrNorm() const
	{
		return mX * mX + mY * mY + mZ * mZ;
	}



	inline double Point::Norm() const
	{
		return sqrt(SqrNorm());
	}



	inline Point Point::operator+(const Point &p) const
	{
		return Point(mX + p.mX, mY + p.mY, mZ + p.mZ);
	}



	inline Point Point::operator-(const Point &p) const
	{
		return Point(mX - p.mX, mY - p.mY, mZ - p.mZ);
	}



	inline Point Point::operator*(double a) const
	{
		return Point(a * mX, a * mY, a * mZ);
	}



	inline Point Point::operator-()
	{
		return Point(-mX, -mY, -mZ);
	}



	inline Point& Point::operator+=(const Point &p)
	{
		mX += p.mX;
		mY += p.mY;
		mZ += p.mZ;

		return *this;
	}



	inline Point& Point::operator-=(const Point &p)
	{
		mX -= p.mX;
		mY -= p.mY;
		mZ -= p.mZ;

		return *this;
	}



	inline Point& Point::operator*=(double a)
	{
		mX *= a;
		mY *= a;
		mZ *= a;

		return *this;
	}
	


	inline Point Point::operator/(double a) const
	{
		if (a == 0.0) 
			ThrowException("Point::operator/ : division by zero");
	
		double b = 1.0 / a;
		
		return Point (b * mX, b * mY, b * mZ);
	}



	inline Point& Point::operator/=(double a)
	{
		if (a == 0.0) 
			ThrowException("Point::operator/= : division by zero");

		double b = 1.0 / a;
	
		mX *= b;
		mY *= b;
		mZ *= b;

		return *this;
	}



	inline double Point::Dot(const Point &p) const
	{
		return mX * p.mX + mY * p.mY + mZ * p.mZ;
	}



	inline void Point::Min(Point &pMin) const
	{
		pMin.mX = std::min(pMin.mX, mX);
		pMin.mY = std::min(pMin.mY, mY);
		pMin.mZ = std::min(pMin.mZ, mZ);
	
		return;
	}
	


	inline void Point::Max(Point &pMax) const
	{
		pMax.mX = std::max(pMax.mX, mX);
		pMax.mY = std::max(pMax.mY, mY);
		pMax.mZ = std::max(pMax.mZ, mZ);
	
		return;
	}



	inline double Point::DistanceToPoint(const Point &p) const
	{
		return sqrt(SqrDistanceToPoint(&p));
	}



	inline double Point::DistanceToPoint(const Point *pPt) const
	{
		return sqrt(SqrDistanceToPoint(pPt));
	}
	


	inline double Point::SqrDistanceToPoint(const Point *pPt) const
	{
		double diff = mX - pPt->mX; 
		double sum = diff * diff;
		
		diff = mY - pPt->mY; 
		sum += diff * diff;
		
		diff = mZ - pPt->mZ; 
		sum += diff * diff;
		
		return sum;
	}



	inline void Point::Normalize()
	{
		*this /= Norm();
		return;
	}
	

	
	inline Point Point::Cross(const Point &p) const
	{
		return Point(mY * p.mZ - mZ * p.mY,
					 mZ * p.mX - mX * p.mZ,
					 mX * p.mY - mY * p.mX);
	}
	


	inline double Point::XYCross(const Point *pPt) const
	{
		// cross product assuming the Z component of each point is zero
		// return value is Z component of cross product vector

		return mX * pPt->mY - mY * pPt->mX;
	}
}


#endif // _point_h_	
