#ifndef __TMV_HPP__
#define __TMV_HPP__

class Vec3 {

public:
	
	inline Vec3() {}
	inline Vec3( const Vec3& a ) : m_x(a.m_x), m_y(a.m_y), m_z(a.m_z) {}
	inline Vec3( double x, double y, double z ) : m_x(x), m_y(y), m_z(z) {}
	
	inline double& x() { return m_x; }
	inline double x() const { return m_x; }
	inline double& y() { return m_y; }
	inline double y() const { return m_y; }
	inline double& z() { return m_z; }
	inline double z() const { return m_z; }
	
	double& operator [](int i) { return i==0 ? m_x : (i==1 ? m_y : m_z); }
	double operator [](int i) const { return i==0 ? m_x : (i==1 ? m_y : m_z); }
	

    inline Vec3& operator=(const Vec3 &rhs) {

       // Only do assignment if RHS is a different object from this.
       if (this != &rhs) {
		   m_x = rhs.m_x; m_y = rhs.m_y; m_z = rhs.m_z;
       }

       return *this;
     }

    inline Vec3& operator+=(const Vec3& rhs)
    {
        m_x += rhs.m_x; m_y += rhs.m_y; m_z += rhs.m_z;
        return *this;
    }
    inline Vec3& operator-=(const Vec3& rhs)
    {
        m_x -= rhs.m_x; m_y -= rhs.m_y; m_z -= rhs.m_z;
        return *this;
    }
	inline double normsqr() const {
	    double norm = sqr(m_x) + sqr(m_y) + sqr(m_z);
	    return norm;
	}
	inline double norm() const {
	    return sqrt(normsqr());
	}
	
private:
	double m_x, m_y, m_z;

};

inline const Vec3 operator+(Vec3 lhs, const Vec3& rhs) {
    return lhs += rhs;
}
inline const Vec3 operator-(Vec3 lhs, const Vec3& rhs) {
    return lhs -= rhs;
}
inline const double dot(const Vec3& lhs, const Vec3& rhs) {
	return lhs.x()*rhs.x() + lhs.y()*rhs.y() + lhs.z()*rhs.z();
}
inline const Vec3 cross(const Vec3& lhs, const Vec3& rhs) {
	double c1 = lhs.y()*rhs.z() - lhs.z()*rhs.y();
	double c2 = lhs.z()*rhs.x() - lhs.x()*rhs.z();
	double c3 = lhs.x()*rhs.y() - lhs.y()*rhs.x();
	return Vec3(c1, c2, c3);
}

inline const Vec3 operator*(const double lambda, const Vec3 &v) {
    Vec3 result;
    result.x() = lambda * v.x();
    result.y() = lambda * v.y();
    result.z() = lambda * v.z();
    
    return result;
}


class Mat33 {
public:
    inline Mat33() {}
    inline Mat33(const Mat33& M) :  m_11 (M.m_11), m_12(M.m_12), m_13(M.m_13), m_21(M.m_21), m_22(M.m_22), m_23(M.m_23), m_31(M.m_31), m_32(M.m_32), m_33(M.m_33) {}
    
    inline Mat33(Vec3 v1, Vec3 v2, Vec3 v3) : m_11(v1.x()), m_21(v1.y()), m_31(v1.z()), m_12(v2.x()), m_22(v2.y()), m_32(v2.z()), m_13(v3.x()), m_23(v3.y()), m_33(v3.z()) {}
    
    inline Mat33(double a_11, double a_12, double a_13, double a_21, double a_22, double a_23, double a_31, double a_32, double a_33) :
    m_11(a_11), m_12(a_12), m_13(a_13), m_21(a_21), m_22(a_22), m_23(a_23), m_31(a_31), m_32(a_32), m_33(a_33) {}
    
    inline const Mat33 t() const {
        Mat33 B(m_11, m_21, m_31, m_12, m_22, m_32, m_13, m_23, m_33);
        return B;
    }
    
    inline Mat33& operator=(const Mat33 &rhs) {
        
        // Only do assignment if RHS is a different object from this.
        if (this != &rhs) {
            m_11 = rhs.m_11; m_12 = rhs.m_12; m_13 = rhs.m_13;
            m_21 = rhs.m_21; m_22 = rhs.m_22; m_23 = rhs.m_23;
            m_31 = rhs.m_31; m_32 = rhs.m_32; m_33 = rhs.m_33;
      }
        
        return *this;
    }

    inline const Mat33 inv() const {
        double det = m_11*m_22*m_33 + m_12*m_23*m_31 + m_13*m_21*m_32;
        det -= m_11 * m_23 * m_32 + m_12 * m_21 * m_33 + m_13 * m_22 * m_31;
        
        double A_11 = m_22 * m_33 - m_23*m_32; 
		A_11 = A_11/det;
        
		double A_12 = m_12 * m_33 - m_13* m_32;
		A_12 = -A_12/det;
        
		double A_13 = m_12*m_23 - m_13*m_22;
		A_13 = A_13/det;
        
        double A_21 = m_21*m_33 - m_23*m_31;
		A_21 = - A_21/det;
        
		double A_22 = m_11 * m_33 - m_13*m_31;
		A_22 = A_22/det;
        
		double A_23 = m_11*m_23 - m_13*m_21;
		A_23 = -A_23 /det;
        
        double A_31 = m_21 *m_32 - m_22*m_31;
		A_31 = A_31/det;
        
		double A_32 = m_11*m_32 - m_12*m_31;
		A_32 = -A_32/det;
		
        double A_33 = m_11 *m_22 - m_12*m_21;
		A_33 = A_33/det;
    
    Mat33 A (A_11, A_12, A_13, A_21, A_22, A_23, A_31, A_32, A_33);
    return A;
    
    }
    
    inline Vec3 row1() const {
        return Vec3(m_11, m_12, m_13);
    }
    
    inline Vec3 row2() const {
        return Vec3(m_21, m_22, m_23);
    }
    
    inline Vec3 row3() const {
        return Vec3(m_31, m_32, m_33);
    }
	inline Vec3 row(int j) const {
	    if (j==0) return Vec3(m_11, m_12, m_13);
		if (j==1) return Vec3(m_21, m_22, m_23);
		if (j==2) return Vec3(m_31, m_32, m_33);
		die();
		return Vec3(0.0,0.0,0.0);
	};

	inline const Vec3 operator*( const Vec3 &v) const {
		Vec3 result;
		result.x() = m_11*v.x() + m_12*v.y() + m_13*v.z();
		result.y() = m_21*v.x() + m_22*v.y() + m_23*v.z();
		result.z() = m_31*v.x() + m_32*v.y() + m_33*v.z();
		
		return result;
	}
	
private:
	double m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33;

};

inline Vec3 crop(Vec3 v, Vec3 min, Vec3 max) {
	
	Vec3 ret;
	ret = v;
	ret.x() = std::max( std::min( ret.x(), max.x() ), min.x());
	ret.y() = std::max( std::min( ret.y(), max.y() ), min.y());
	ret.z() = std::max( std::min( ret.z(), max.z() ), min.z());
	
	return ret;
	
	
}

#endif