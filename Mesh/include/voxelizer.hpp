#pragma once

#include "read_stl.hpp"
#include "BasicGeometry.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <unordered_set>
#include <fstream>
#include <ctime>

namespace mesh {

    template<typename T>
    class Voxelizer {
    public:
        Voxelizer(const std::string& stl_file_name, const T dx)
            : _dx(dx) {
            // Randomness.
            srand(static_cast<unsigned>(time(0)));
            // Load triangles from the stl file.
            std::vector<Vector3<T>> normals;
            if (!ReadSTL(stl_file_name, _triangles, normals)) {
                std::cout << "ERROR: cannot read " << stl_file_name << std::endl;
                return;
            }
            // Compute the bounding box of _triangle and save the results into _pmin.
            _pmin = _triangles[0][0];
            Vector3<T> pmax = _triangles[0][0];
            for (const auto& triangle : _triangles)
                for (const auto& v : triangle) {
                    _pmin = _pmin.cwiseMin(v);
                    pmax = pmax.cwiseMax(v);
                }
            for (int i = 0; i < 3; ++i) {
                _pmin[i] -= _dx;
                pmax[i] += _dx;
            }
            // Compute the number of voxels along each direction.
            for (int i = 0; i < 3; ++i)
                _nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;
            // Initialize the voxel array.
            _voxels = std::vector<std::vector<std::vector<bool>>>(_nvoxel.x(),
                std::vector<std::vector<bool>>(_nvoxel.y(),
                    std::vector<bool>(_nvoxel.z(), false)));
        }

        const Vector3<T> pmin() const { return _pmin; }
        const T dx() const { return _dx; }
        const Vector3<T> pmax() const { return _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx; }
        const Vector3<int> voxel_num() const { return _nvoxel; }

        void BasicVoxelization() {
            /* Assignment 2, Part 2.1. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];

            /* Find the intersections by checking each ray with every triangle. */
            std::vector<std::vector<std::vector<T>>> intersections;
            for (int i = 0; i < nx; ++i) {
                std::vector<std::vector<T>> everythingAlongY;
                for (int j = 0; j < ny; ++j) {
                    std::vector<T> everythingAlongZ;
                    for (int tri_index = 0; tri_index < _triangles.size(); ++tri_index) {
                        Vector3<T> origin = Vector3<T>(pmin()[0] + i*dx() + dx()/2, pmin()[1] + j*dx() + dx()/2, pmin()[2]);
                        Vector3<T> dir = Vector3<T>(0, 0, 1);
                        geometry::Triangle<T> triangle = geometry::Triangle<T>(_triangles[tri_index][0], _triangles[tri_index][1], _triangles[tri_index][2]);
                        T t = triangle.IntersectRay(origin, dir);
                        if (t >= 0) {
                            everythingAlongZ.push_back(t);
                        }
                    }
                    everythingAlongY.push_back(everythingAlongZ);
                }
                intersections.push_back(everythingAlongY);
            }

            /* For each x, y coordinate, check along the ray in the z direction and decide for each voxel whether it should
             * be filled or not. This is based on the intersections found in the part above. The first filled voxel is
             * the one above the first intersection, then all of them are filled until another intersection is hit.*/
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    if (intersections[i][j].size() > 0) {
                        std::sort(intersections[i][j].begin(), intersections[i][j].end());
                        T fill = false;
                        T curr_index = 0;
                        T curr_value = intersections[i][j][curr_index];
                        if (curr_value < 1e-6) {
                            fill = true;
                        }
                        for (int k = 0; k < nz; ++k) {
                            if ((k + 0.5) * dx() > curr_value) {
                                curr_index += 1;
                                if (curr_index >= intersections[i][j].size()) {
                                    fill = false;
                                } else {
                                    curr_value = intersections[i][j][curr_index];
                                    if (fill) {
                                        fill = false;
                                    } else {
                                        fill = true;
                                    }
                                }
                            }
                            _voxels[i][j][k] = fill;
                        }
                    } else {
                        for (int k = 0; k < nz; ++k) {
                            _voxels[i][j][k] = false;
                        }
                    }
                }
            }
        }

        void AdvancedVoxelization() {
            /* Assignment 2, Part 2.2. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];

            /* Set up the datastructure to store the intersections. */
            std::vector<std::vector<std::vector<T>>> intersections;
            for (int i = 0; i < nx; ++i) {
                std::vector<std::vector<T>> everythingAlongY;
                for (int j = 0; j < ny; ++j) {
                    std::vector<T> everythingAlongZ;
                    everythingAlongY.push_back(everythingAlongZ);
                }
                intersections.push_back(everythingAlongY);
            }

            /* For each triangle, shoot rays from only the voxels that the projection of the triangle on the xy plane
             * overlaps with.*/
            for (int tri_index; tri_index < _triangles.size(); ++tri_index) {
                Vector3<T> vertex_0 = _triangles[tri_index][0];
                Vector3<T> vertex_1 = _triangles[tri_index][1];
                Vector3<T> vertex_2 = _triangles[tri_index][2];
                geometry::Triangle<T> triangle = geometry::Triangle<T>(vertex_0, vertex_1, vertex_2);
                T x_min = std::min(vertex_0[0], std::min(vertex_1[0], vertex_2[0]));
                T x_max = std::max(vertex_0[0], std::max(vertex_1[0], vertex_2[0]));
                T y_min = std::min(vertex_0[1], std::min(vertex_1[1], vertex_2[1]));
                T y_max = std::max(vertex_0[1], std::max(vertex_1[1], vertex_2[1]));

                for (int i = floor((x_min - pmin()[0]) / dx() - 0.5); i < ceil((x_max - pmin()[0]) / dx() + 0.5); ++i) {
                    for (int j = floor((y_min - pmin()[1]) / dx() - 0.5); j < ceil((y_max - pmin()[1]) / dx() + 0.5); ++j) {
                        Vector3<T> origin = Vector3<T>(pmin()[0] + i*dx() + dx()/2, pmin()[1] + j*dx() + dx()/2, pmin()[2]);
                        Vector3<T> dir = Vector3<T>(0, 0, 1);
                        T t = triangle.IntersectRay(origin, dir);
                        if (t >= 0 ) {
                            intersections[i][j].push_back(t);
                        }
                    }
                }
            }

            /* Decide which voxels to fill in the same way as in the Basic Voxelization. */
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    if (intersections[i][j].size() > 0) {
                        std::sort(intersections[i][j].begin(), intersections[i][j].end());
                        T fill = false;
                        T curr_index = 0;
                        T curr_value = intersections[i][j][curr_index];
                        if (curr_value < 1e-6) {
                            fill = true;
                        }
                        for (int k = 0; k < nz; ++k) {
                            if ((k + 0.5) * dx() > curr_value) {
                                curr_index += 1;
                                if (curr_index >= intersections[i][j].size()) {
                                    fill = false;
                                } else {
                                    curr_value = intersections[i][j][curr_index];
                                    if (fill) {
                                        fill = false;
                                    } else {
                                        fill = true;
                                    }
                                }
                            }
                            _voxels[i][j][k] = fill;
                        }
                    } else {
                        for (int k = 0; k < nz; ++k) {
                            _voxels[i][j][k] = false;
                        }
                    }
                }
            }
        }

        void AdvancedVoxelizationWithApproximation() {
            /* Assignment 2, Part 2.3. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];

            int numRays = 11;

            /* This will be where the different ray directions are stored. */
            std::vector<Vector3<T>> rays;

            /* Set up data structure to store all of the intersections for each ray. */
            std::vector<std::vector<std::vector<std::vector<T>>>> intersections;

            /* Where the rotation matrices for each ray will be stored. This is the rotation from the ray to vertical. */
            std::vector<Eigen::Quaternion<T>> rotationMatrices;

            /* Store the dimensions of the voxel grid for each of the rotations. Each set of dimensions is stored as
             * a Vector3<T> where the first number is the number of voxels in the new x direction, and same with y and z.*/
            std::vector<Vector3<T>> voxelDimensionsAllRays;

            /* This call sets up each of the above data structures. */
            PopulateDataStructures(numRays, rays, intersections, rotationMatrices, voxelDimensionsAllRays);

            /* Find the intersections for each triangle with each ray. */
            for (int tri_index; tri_index < _triangles.size(); ++tri_index) {
                std::vector<Vector3<T>> triangleVector = _triangles[tri_index];
                for (int rayIndex = 0; rayIndex < numRays; ++rayIndex) {
                    RayIntersections(rayIndex, triangleVector, intersections, rays[rayIndex]);
                }
            }

            /* Set up data structure to store whether each voxel should be full or not for each ray. */
            std::vector<std::vector<std::vector<std::vector<T>>>> voxelsForAllRays;
            for (int i = 0; i < nx; ++i) {
                std::vector<std::vector<std::vector<T>>> everythingAlongY;
                for (int j = 0; j < ny; ++j) {
                    std::vector<std::vector<T>> everythingAlongZ;
                    for (int k = 0; k < nz; ++k){
                        std::vector<T> forAllRays;
                        everythingAlongZ.push_back(forAllRays);
                    }
                    everythingAlongY.push_back(everythingAlongZ);
                }
                voxelsForAllRays.push_back(everythingAlongY);
            }

            /* Decide for each voxel if it should be filled or not. */
            for (int rayIndex = 0; rayIndex < numRays; ++rayIndex) {
                CheckRayIntersections(rayIndex,
                                      rays,
                                      intersections,
                                      rotationMatrices,
                                      voxelDimensionsAllRays,
                                      voxelsForAllRays);
            }

            /* Decide which voxels to fill and not fill. This is done by checking the result for all of the different
             * rays. If more of the rays say that the voxel should be filled, then the voxel should be filled. If a
             * majority say it should not be filled, then it will not be filled. If there are an equal number, then
             * the default is to not be filled.*/
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        std::vector<T> allFills = voxelsForAllRays[i][j][k];
                        int numTrue = 0;
                        int numFalse = 0;
                        for (int b = 0; b < allFills.size(); ++b) {
                            if (allFills[b]) {
                                numTrue += 1;
                            }
                        }
                        if (numTrue > numFalse) {
                            _voxels[i][j][k] = true;
                        } else {
                            _voxels[i][j][k] = false;
                        }
                    }
                }
            }
        }

        /* This method sets up the datastructures to be the correct size so that they can be indexed into and
         * populated later. */
        void PopulateDataStructures(int numRays,
                                    std::vector<Vector3<T>>& rays,
                                    std::vector<std::vector<std::vector<std::vector<T>>>>& intersections,
                                    std::vector<Eigen::Quaternion<T>>& rotationMatrices,
                                    std::vector<Vector3<T>>& voxelDimensionsAllRays) {
            for (int ray = 0; ray < numRays; ++ray) {
                Vector3<T> randomRay = Vector3<T>(rand(), rand(), rand());
                rays.push_back(randomRay);

                Vector3<T> verticalRay = Vector3<T>(0, 0, 1);

                Eigen::Quaternion<T> q = Eigen::Quaternion<T>::FromTwoVectors(randomRay, verticalRay);
                q.matrix();
                rotationMatrices.push_back(q);

                Vector3<T> newMin = q * pmin();
                Vector3<T> newMax = q * pmax();

                Vector3<T> voxelDimensions = Vector3<T>(ceil((newMax[0] - newMin[0])/dx()),
                                                        ceil((newMax[1] - newMin[1])/dx()),
                                                        ceil((newMax[2] - newMin[2])/dx()));
                voxelDimensionsAllRays.push_back(voxelDimensions);

                std::vector<std::vector<std::vector<T>>> everythingForThisRay;
                for (int i = 0; i < ceil((newMax[0] - newMin[0])/dx()) + 1; ++i) {
                    std::vector<std::vector<T>> everythingAlongY;
                    for (int j = 0; j < ceil((newMax[1] - newMin[1])/dx()) + 1; ++j) {
                        std::vector<T> everythingAlongZ;
                        everythingAlongY.push_back(everythingAlongZ);

                    }
                    everythingForThisRay.push_back(everythingAlongY);
                }
                intersections.push_back(everythingForThisRay);
            }
        }

        /* This method determines whether a voxel should be filled or not filled in the rotated reference frame and
         * translates it into the original reference frame.
         * This is done based on the intersections already found in the same way done in the Advanced Voxelization.*/
        void CheckRayIntersections(int rayIndex,
                                   std::vector<Vector3<T>>& rays,
                                   std::vector<std::vector<std::vector<std::vector<T>>>>& intersections,
                                   std::vector<Eigen::Quaternion<T>>& rotationMatrices,
                                   std::vector<Vector3<T>>& voxelDimensionsAllRays,
                                   std::vector<std::vector<std::vector<std::vector<T>>>>& voxelsForAllRays) {
            Vector3<T> ray = rays[rayIndex];

            Eigen::Quaternion<T> q = rotationMatrices[rayIndex];
            Vector3<T> newMin = q * pmin();
            Vector3<T> newMax = q * pmax();

            Vector3<T> verticalRay = Vector3<T>(0, 0, 1);
            Eigen::Quaternion<T> rotateBack = Eigen::Quaternion<T>::FromTwoVectors(verticalRay, ray);
            rotateBack.matrix();

            for (int i = 0; i < voxelDimensionsAllRays[rayIndex][0]; ++i) {
                for (int j = 0; j < voxelDimensionsAllRays[rayIndex][1]; ++j) {
                    if (intersections[i][j].size() > 0) {
                        std::sort(intersections[rayIndex][i][j].begin(), intersections[rayIndex][i][j].end());

                        T fill = false;
                        T curr_index = 0;
                        T curr_value = intersections[rayIndex][i][j][curr_index];
                        if (curr_value < 1e-6) {
                            fill = true;
                        }
                        for (int k = 0; k < voxelDimensionsAllRays[rayIndex][2]; ++k) {

                            if ((k + 0.5) * dx() > curr_value) {
                                curr_index += 1;
                                if (curr_index >= intersections[rayIndex][i][j].size()) {
                                    fill = false;
                                } else {
                                    curr_value = intersections[rayIndex][i][j][curr_index];
                                    if (fill) {
                                        fill = false;
                                    } else {
                                        fill = true;
                                    }
                                }
                            }

                            /* The below section is to determine which voxel in the original reference frame
                             * the voxel currently looked at relates to. */
                            Vector3<T> center = Vector3<T>(newMin[0] + (i + 0.5)*dx(),
                                                           newMin[1] + (j + 0.5) * dx(),
                                                           newMin[2] + (k + 0.5) * dx());

                            Vector3<T> setVoxel = PointFallsInVoxel(center);
                            voxelsForAllRays[setVoxel[0]][setVoxel[1]][setVoxel[2]].push_back(fill);
                        }
                    } else {
                        for (int k = 0; k < voxelDimensionsAllRays[rayIndex][2]; ++k) {
                            /* The below section is to determine which voxel in the original reference frame
                             * the voxel currently looked at relates to. */
                            Vector3<T> center = Vector3<T>(newMin[0] + (i + 0.5)*dx(),
                                                           newMin[1] + (j + 0.5) * dx(),
                                                           newMin[2] + (k + 0.5) * dx());

                            Vector3<T> setVoxel = PointFallsInVoxel(center);
                            voxelsForAllRays[setVoxel[0]][setVoxel[1]][setVoxel[2]].push_back(false);
                        }
                    }
                }
            }
        }

        /* This method takes a point representing the center of a voxel in the rotated frame and maps it to
         * a voxel in the original reference frame. The indexes into the original voxel in _voxels is returned as
         * a Vector3<T>  */
        Vector3<T> PointFallsInVoxel(Vector3<T>& point) {
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            bool xSet = false;
            T xFinal;
            bool ySet = false;
            T yFinal;
            bool kSet = false;
            T kFinal;

            for (int i = 0; i < nx; ++i) {
                T x = pmin()[0] + (i + 0.5)*dx();
                if (x < point[0]) {
                    continue;
                } else if (!xSet) {
                    xFinal = x;
                }
                for (int j = 0; j < ny; ++j) {
                    T y = pmin()[1] + (j + 0.5)*dx();
                    if (y < point[1]) {
                        continue;
                    } else if (!ySet) {
                        yFinal = y;
                    }
                    for (int k = 0; k < nz; ++k) {
                        T y = pmin()[2] + (k + 0.5)*dx();
                        if (k < point[2]) {
                            continue;
                        } else if (!kSet) {
                            kFinal = k;
                            return Vector3<T>(xFinal, yFinal, kFinal);
                        }
                    }
                }
            }
        }

        /* This method does the calculations for intersecting rays with a single given triangle. */
        void RayIntersections(int rayIndex,
                              std::vector<Vector3<T>>& triangleVector,
                              std::vector<std::vector<std::vector<std::vector<T>>>>& intersections,
                              Vector3<T>& ray) {
            Vector3<T> verticalRay = Vector3<T>(0, 0, 1);

            Vector3<T> vertex_0 = triangleVector[0];
            Vector3<T> vertex_1 = triangleVector[1];
            Vector3<T> vertex_2 = triangleVector[2];

            Eigen::Quaternion<T> q = Eigen::Quaternion<T>::FromTwoVectors(ray, verticalRay);
            q.matrix();

            Vector3<T> vertex_0_rotated = q * vertex_0;
            Vector3<T> vertex_1_rotated = q * vertex_1;
            Vector3<T> vertex_2_rotated = q * vertex_2;

            geometry::Triangle<T> triangle = geometry::Triangle<T>(vertex_0_rotated, vertex_1_rotated, vertex_2_rotated);
            T x_min = std::min(vertex_0_rotated[0], std::min(vertex_1_rotated[0], vertex_2_rotated[0]));
            T x_max = std::max(vertex_0_rotated[0], std::max(vertex_1_rotated[0], vertex_2_rotated[0]));
            T y_min = std::min(vertex_0_rotated[1], std::min(vertex_1_rotated[1], vertex_2_rotated[1]));
            T y_max = std::max(vertex_0_rotated[1], std::max(vertex_1_rotated[1], vertex_2_rotated[1]));

            Vector3<T> min = q * pmin();

            for (int i = floor((x_min - min[0]) / dx() - 0.5); i < ceil((x_max - min[0]) / dx() + 0.5); ++i) {
                for (int j = floor((y_min - min[1]) / dx() - 0.5); j < ceil((y_max - min[1]) / dx() + 0.5); ++j) {
                    Vector3<T> origin = Vector3<T>(min[0] + i*dx() + dx()/2, min[1] + j*dx() + dx()/2, min[2]);
                    Vector3<T> dir = Vector3<T>(0, 0, 1);
                    T t = triangle.IntersectRay(origin, dir);

                    if (t >= 0 ) {
                        /* The segmentation fault is happening from the line below. Unfortunately I can't seem to fix it. */
                        intersections[rayIndex][i][j].push_back(t);
                    }
                }
            }
        }

        void WriteVoxelToMesh(const std::string& stl_file_name) const {
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            std::vector<std::vector<Vector3<int>>> faces;
            std::vector<Vector3<int>> corners({
                Vector3<int>(0, 0, 0),
                Vector3<int>(0, 0, 1),
                Vector3<int>(0, 1, 0),
                Vector3<int>(0, 1, 1),
                Vector3<int>(1, 0, 0),
                Vector3<int>(1, 0, 1),
                Vector3<int>(1, 1, 0),
                Vector3<int>(1, 1, 1)
            });
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k) {
                        if (!_voxels[i][j][k]) continue;
                        // Check -x direction.
                        Vector3<int> cmin(i, j, k);
                        if (i == 0 || !_voxels[i - 1][j][k]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[1], cmin + corners[3] });
                            faces.push_back({ cmin + corners[0], cmin + corners[3], cmin + corners[2] });
                        }
                        if (i == nx - 1 || !_voxels[i + 1][j][k]) {
                            faces.push_back({ cmin + corners[4], cmin + corners[6], cmin + corners[7] });
                            faces.push_back({ cmin + corners[4], cmin + corners[7], cmin + corners[5] });
                        }
                        if (j == 0 || !_voxels[i][j - 1][k]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[4], cmin + corners[5] });
                            faces.push_back({ cmin + corners[0], cmin + corners[5], cmin + corners[1] });
                        }
                        if (j == ny - 1 || !_voxels[i][j + 1][k]) {
                            faces.push_back({ cmin + corners[2], cmin + corners[3], cmin + corners[7] });
                            faces.push_back({ cmin + corners[2], cmin + corners[7], cmin + corners[6] });
                        }
                        if (k == 0 || !_voxels[i][j][k - 1]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[2], cmin + corners[6] });
                            faces.push_back({ cmin + corners[0], cmin + corners[6], cmin + corners[4] });
                        }
                        if (k == nz - 1 || !_voxels[i][j][k + 1]) {
                            faces.push_back({ cmin + corners[5], cmin + corners[7], cmin + corners[3] });
                            faces.push_back({ cmin + corners[5], cmin + corners[3], cmin + corners[1] });
                        }
                    }
            std::ofstream fout(stl_file_name);
            fout << "solid vcg" << std::endl;
            for (const auto& f : faces) {
                std::vector<Vector3<T>> p;
                for (const auto& fi : f) {
                    Vector3<T> v = _pmin + fi.cast<T>() * _dx;
                    p.push_back(v);
                }
                const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
                fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
                fout << "    outer loop" << std::endl;
                for (const auto& v : p) {
                    fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
                }
                fout << "    endloop" << std::endl;
                fout << "  endfacet" << std::endl;
            }
            fout << "endsolid vcg" << std::endl;
        }

    private:
        std::vector<std::vector<Vector3<T>>> _triangles;
        T _dx;  // The size of each voxel.
        Vector3<T> _pmin;    // The min and max corner of the bounding box.
        Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.
        std::vector<std::vector<std::vector<bool>>> _voxels;   // True <-> voxel is occupied.
    };

}
