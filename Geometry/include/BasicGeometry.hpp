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

        /* Implement triangle plane intersection.
            Input is a plane, output is a list of intersection points
            Hints:
                1. Take some considerations: handle duplicate intersections (when plane intersect with 
                    triangle at corner but you find them by edge-plane intersection).
                2. You can modify the function input/return type to any format you want.
        */
        std::vector<Vector3<T>> IntersectPlane(Plane<T> p) {
            /* Implement your code here */
	        std::vector<Vector3<T>> edge;

	        for (int i = 0; i < 3; ++i) {
	            Vector3<T> v_1 = _vertices[i];
	            T dist_1;
	            if (p.onPlane(v_1, dist_1)) {
	                edge.push_back(v_1);
	            }
	        }
	        for (int i = 0; i < 3; ++i) {
	            Vector3<T> v_1 = _vertices[i];
	            T dist_1;
	            if (p.onPlane(v_1, dist_1)) {
                    continue;
                }
	            int j = (i + 1)%3;
                Vector3<T> v_2 = _vertices[j];
                T dist_2;
                if (p.onPlane(v_2, dist_2)) {
                    continue;
                }
                if (i == j) {
                    continue;
                } else if (dist_1 * dist_2 < 0) {
                    if (dist_1 < 0) {
                        dist_1 = -1 * dist_1;
                    }
                    if (dist_2 < 0) {
                        dist_2 = -1 * dist_2;
                    }
                    Vector3<T> triangle_edge = v_2 - v_1;
                    Vector3<T> plane_intersection = v_1 + triangle_edge*(dist_1/(dist_1 + dist_2));

                    edge.push_back(plane_intersection);
                }
	        }
	        return edge;

//	        int numOnPlane = 0;
//
//	        for (int i = 0; i < 3; ++i) {
//		        Vector3<T> v_1 = _vertices[i];
//		        T dist_1;
//		        if (p.onPlane(v_1, dist_1)) {
//		            edge.push_back(v_1);
//		            numOnPlane += 1;
//		        } else {
//		            for (int j = 0; j < 3; ++j) {
//			            Vector3<T> v_2 = _vertices[j];
//			            T dist_2;
//			            p.onPlane(v_2, dist_2);
//		                if (i == j) {
//			                continue;
//		                } else if (dist_1 * dist_2 < 0) {
//			                if (dist_1 < 0) {
//				                dist_1 = -1 * dist_1;
//			                }
//			                if (dist_2 < 0) {
//				                dist_2 = -1 * dist_2;
//			                }
//			                Vector3<T> triangle_edge = v_1 - v_2;
//			                Vector3<T> plane_intersection = triangle_edge*(dist_1/(dist_1 + dist_2));
//
//			                edge.push_back(plane_intersection);
//		                }
//		            }
//		            if (edge.size() > 0) {
//		                return edge;
//		            }
//		        }
//		        if (numOnPlane == 2) {
//		            return edge;
//		        }
//	        }
        }
        
    private:
        Vector3<T> _vertices[3];
    };
}
