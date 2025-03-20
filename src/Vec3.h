/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 3, 2022                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include <cassert>
#include <iosfwd>


class Vec3 {
public:
    /**
     * Create a Vec3 whose elements are all 0.
     */
    Vec3() {
        data[0] = data[1] = data[2] = 0.0;
    }
    /**
     * Create a Vec3 with specified x, y, and z components.
     */
    Vec3(double x, double y, double z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }
    double& operator[](int index) {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    bool operator==(const Vec3& rhs) const {
        return (data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
    }

    bool operator!=(const Vec3& rhs) const {
        return (data[0] != rhs[0] || data[1] != rhs[1] || data[2] != rhs[2]);
    }
    
    // Arithmetic operators
    
    // unary plus
    Vec3 operator+() const {
        return Vec3(*this);
    }
    
    // plus
    Vec3 operator+(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return Vec3(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }
    
    Vec3& operator+=(const Vec3& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    // unary minus
    Vec3 operator-() const {
        const Vec3& lhs = *this;
        return Vec3(-lhs[0], -lhs[1], -lhs[2]);
    }
    
    // minus
    Vec3 operator-(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return Vec3(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
    }

    Vec3& operator-=(const Vec3& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

    // scalar product
    Vec3 operator*(double rhs) const {
        const Vec3& lhs = *this;
        return Vec3(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
    }

    Vec3& operator*=(double rhs) {
        data[0] *= rhs;
        data[1] *= rhs;
        data[2] *= rhs;
        return *this;
    }

    // scalar division
    Vec3 operator/(double rhs) const {
        const Vec3& lhs = *this;
        double scale = 1.0/rhs;
        return Vec3(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale);
    }

    Vec3& operator/=(double rhs) {
        double scale = 1.0/rhs;
        data[0] *= scale;
        data[1] *= scale;
        data[2] *= scale;
        return *this;
    }
    
    // dot product
    double dot(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
    }

    // cross product
    Vec3 cross(const Vec3& rhs) const {
        return Vec3(data[1]*rhs[2]-data[2]*rhs[1], data[2]*rhs[0]-data[0]*rhs[2], data[0]*rhs[1]-data[1]*rhs[0]);
    }
    
private:
    double data[3];
};

static Vec3 operator*(double lhs, Vec3 rhs) {
    return Vec3(rhs[0]*lhs, rhs[1]*lhs, rhs[2]*lhs);
}

template <class CHAR, class TRAITS>
std::basic_ostream<CHAR,TRAITS>& operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Vec3& v) {
    o<<'['<<v[0]<<", "<<v[1]<<", "<<v[2]<<']';
    return o;
}

