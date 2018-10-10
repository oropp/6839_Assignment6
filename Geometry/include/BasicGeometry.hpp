#pragma once
#include <Eigen/Dense>
#include <vector>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    // the plane is represented by (x - _p) /dot _normal = 0
    template <typename T>
    class Plane {
    public:
        Plane(Vector3<T> p, Vector3<T> normal) {
            _p = p;
            _normal = normal;
            _normal.normalize();
        }

        Vector3<T>& p() { return _p; }
        Vector3<T>& normal() { return _normal; }
        
        // return if the point is on plane
        // also fill parameter dist field as the signed distance from point to plane
        bool onPlane(Vector3<T> point, T& dist) {
            dist = (point - _p).dot(_normal);
            if (std::fabs(dist) < 1e-6) {
                return true;
            } else {
                return false;
            }
        }

    private:
        Vector3<T> _p;
        Vector3<T> _normal;
    };

    template <typename T>
    class Triangle {
    public:
        Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }

        std::vector<Vector3<T>> IntersectPlane(Plane<T> p) {}

        // Assignment 2: implement ray-triangle intersection.
        // The ray is defined as r(t) = origin + t * dir.
        // You should return a scalar t such that r(t) is the intersection point. Which value
        // to return for the case of no intersections is up to you. You can also change the
        // signature of this function as you see fit.
        const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const {
            /* Assignment 2. */
            /* Implement your code here */

            Vector3<T> A = _vertices[0];
            Vector3<T> B = _vertices[1];
            Vector3<T> C = _vertices[2];

            /* Find the intersection of the ray and the plane */
            Vector3<T> BA = B - A;
            Vector3<T> CA = C - A;
            Vector3<T> n = (BA.cross(CA))/((BA.cross(CA)).norm());
            T D_plane = n.dot(A);
            T t = (D_plane - n.dot(origin))/(n.dot(dir));

            /* Check if the intersection point is in the triangle */
            Vector3<T> intersectionPoint = origin + t*dir;

            Vector3<T> AB = A - B;
            Vector3<T> AIntersectionPoint = A - intersectionPoint;
            Vector3<T> crossA = AB.cross(AIntersectionPoint);

            Vector3<T> BC = B - C;
            Vector3<T> BIntersectionPoint = B - intersectionPoint;
            Vector3<T> crossB = BC.cross(BIntersectionPoint);

            Vector3<T> CIntersectionPoint = C - intersectionPoint;
            Vector3<T> crossC = CA.cross(CIntersectionPoint);

            if (n.dot(crossA) >= 0 && n.dot(crossB) >= 0 && n.dot(crossC) >= 0) {
                return t;
            } else if (n.dot(crossA) < 0 && n.dot(crossB) < 0 && n.dot(crossC) < 0) {
                return t;
            }

            return -1.0;
        }

    private:
        Vector3<T> _vertices[3];
    };
}
