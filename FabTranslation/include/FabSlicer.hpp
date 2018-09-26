#pragma once
#include "tri_mesh.hpp"
#include "BasicGeometry.hpp"
#include "IntervalTree.hpp"
#include "cinolib/meshes/meshes.h"
#include <ctime>

namespace fab_translation {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    template <typename T>
    class FabSlicer {
        
    public:
        FabSlicer(mesh::TriMesh<T> tri_mesh, T bottom, T top, T dx,
            T infill_dx)
            : _tri_mesh(tri_mesh), _bottom(bottom), _top(top), _dx(dx), _infill_dx(infill_dx) {

            /* Implement your code here */
            /* 1. Initialize your variables
               2. Build interval tree */

	        for (T h = bottom; h < top; h += dx) {
		        z_map.push_back(std::vector<std::vector<Vector3<T>>>());
	        }

	        for (int i = 0; i < tri_mesh.elements().size(); ++i) {
		        std::vector<Vector3<T>> triangle;
		        for (int j = 0; j < 3; j++) {
		            triangle.push_back(tri_mesh.vertices(tri_mesh.elements(i)[j]));
		        }
		        float min_z = std::min(triangle[0][2], std::min(triangle[1][2], triangle[2][2]));
		        float max_z = std::max(triangle[0][2], std::max(triangle[1][2], triangle[2][2]));

		        int bottom_layer = floor((min_z - bottom)/dx);
		        int top_layer = ceil((max_z - bottom)/dx);

		        for (int j = bottom_layer; j <= top_layer; ++j) {
		            z_map[j].push_back(triangle);
		        }
	        }
        }

        /* Main entrance for FabSlicer
            return contour and infill_edges, which can be directed sent to corresponding visualization function to visualize in MeshLab
            1. each contour contains a list of point in loop order, each layer can have multiple contours, so the contour is vector<vector<vector<Point>>>.
            2. infill edges is a edge soup for each layer, thus it's vector<vector<Edge>>
        */
        void RunTranslation(std::vector<std::vector<std::vector<Vector3<T>>>>& contour,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {

            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> intersection_edges;
            int t_0 = clock();
            Slicing_bruteforce(_tri_mesh, intersection_edges);
            int t_1 = clock();


            Slicing_accelerated(_tri_mesh, intersection_edges);
            int t_2 = clock();

            printf("bruteforce: %d\n", t_1 - t_0);
            printf("accelerated: %d\n", t_2 - t_1);
            CreateContour(_tri_mesh, intersection_edges, contour);
            VisualizeContour("/home/computationalfabrication/CompFab/ComputationalFabrication/data/bunny-contour.ply", 0.001, contour);
            Infill(contour, infill_edges);
            VisualizeInfill("/home/computationalfabrication/CompFab/ComputationalFabrication/data/bunny-infill.ply", 0.001, infill_edges);
        }

        /* Slicing algorithms
            goal: slice the triangle mesh by a set of parallel planes, 
                  output an intersection edge soup for each layer */
        void Slicing_bruteforce(mesh::TriMesh<T>& tri_mesh, 
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> &intersection_edges) {

            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();
            std::vector<Eigen::Vector3i>& edges = tri_mesh.edges();

            intersection_edges.clear();

            for (T h = _bottom; h <= _top; h += _dx) {
                std::vector<std::pair<Vector3<T>, Vector3<T>>> intersections_one_plane;
                intersections_one_plane.clear();

                geometry::Plane<T> plane(Vector3<T>(0, 0, h), Vector3<T>(0, 0, 1));
                for (int i = 0;i < elements.size();++i) {
                    geometry::Triangle<T> triangle(vertices[elements[i](0)], vertices[elements[i](1)], vertices[elements[i](2)]);
                    std::vector<Vector3<T>> intersections = triangle.IntersectPlane(plane);

                    /* Implement your code here */
                    /* What kinds of intersections should be added into intersection edge list? */
		            std::pair<Vector3<T>, Vector3<T>> edge;
		            int size = intersections.size();
		            if (size == 2) {
		                edge = std::make_pair(intersections[0], intersections[1]);
		                intersections_one_plane.push_back(edge);
		            }
                }

                intersection_edges.push_back(intersections_one_plane);
            }
        }

        void Slicing_accelerated(mesh::TriMesh<T>& tri_mesh,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> &intersection_edges) {
            
            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();
            std::vector<Eigen::Vector3i>& edges = tri_mesh.edges();

            intersection_edges.clear();
            int i = 0;
            for (T h = _bottom; h <= _top; h += _dx) {
                /* Implement your code here */
                /* Retrieve candidate triangle list */
                std::vector<std::vector<Vector3<T>>>& candidates = z_map[i];

                std::vector<std::pair<Vector3<T>, Vector3<T>>> intersections_one_plane;
                intersections_one_plane.clear();

                geometry::Plane<T> plane(Vector3<T>(0, 0, h), Vector3<T>(0, 0, 1));
                for (int ii = 0;ii < candidates.size();++ii) {

                    std::vector<Vector3<T>>& tri = candidates[ii];
                    geometry::Triangle<T> triangle(tri[0], tri[1], tri[2]);
                    std::vector<Vector3<T>> intersections = triangle.IntersectPlane(plane);
                    
                    /* Implement your code here */
                    /* What kinds of intersections should be added into intersection edge list? */
                    std::pair<Vector3<T>, Vector3<T>> edge;
                    int size = intersections.size();
                    if (size == 2) {
                        edge = std::make_pair(intersections[0], intersections[1]);
                        intersections_one_plane.push_back(edge);
                    }
                }

                intersection_edges.push_back(intersections_one_plane);
                i += 1;
            }
        }

        /* Find contours
            Goal: Given an intersetion edge soup for each layer, link those edges one by one to form contour loops.
                  Each layer probably has several disjoint contour loops. */
        void CreateContour(mesh::TriMesh<T>& tri_mesh,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> &intersection_edges,
            std::vector<std::vector<std::vector<Vector3<T>>>>& contours) {
            
            /* Implement your code here */
            /* Input is a edge soup, your task is to generate the loop by linking those edges one by one.
               Thinking about how to find two edge sharing a same end point. set a threshold? or a more clever way? */
	        for (int i = 0; i < intersection_edges.size(); ++i) {
		        std::vector<std::pair<Vector3<T>, Vector3<T>>> intersection_edges_plane = intersection_edges[i];

		        std::vector<std::vector<Vector3<T>>> contours_plane;

	    	    while (intersection_edges_plane.size() > 0) {
		            std::vector<Vector3<T>> contour;
	    	        contour.push_back(intersection_edges_plane[0].first);
	    	        contour.push_back(intersection_edges_plane[0].second);
                    intersection_edges_plane.erase(intersection_edges_plane.begin(), intersection_edges_plane.begin() + 1);
	
	    	        bool contour_complete = false;

	    	        while (!contour_complete) {
		                AddToContour(contour, intersection_edges_plane, contour_complete);
	   	            }
	    	        contours_plane.push_back(contour);
	            }
		        contours.push_back(contours_plane);
	        }
        }

	void AddToContour(std::vector<Vector3<T>>& contour,
	                  std::vector<std::pair<Vector3<T>,
	                  Vector3<T>>>& intersection_edges_plane,
	                  bool& contour_complete) {

	    T closest_distance = 99999999.0f;
	    Vector3<T> closest_vertex;
	    Vector3<T> connected_vertex;
	    int closest_index = -1;
        if (intersection_edges_plane.size() == 0) {
            contour_complete = true;
            return;
        }
	    for (int i = 0; i < intersection_edges_plane.size(); ++i) {
		    Vector3<T> first_v = intersection_edges_plane[i].first;
		    Vector3<T> second_v = intersection_edges_plane[i].second;

		    T dist_1 = (first_v - contour[contour.size() - 1]).norm();
	    	T dist_2 = (second_v - contour[contour.size() - 1]).norm();

		    if (dist_1 < closest_distance) {
		        closest_distance = dist_1;
		        closest_vertex = first_v;
		        connected_vertex = second_v;
		        closest_index = i;
		    }
		    if (dist_2 < closest_distance) {
		        closest_distance = dist_2;
		        closest_vertex = second_v;
		        connected_vertex = first_v;
		        closest_index = i;
		    }
	    }
	    T dist_first = (contour[0] - contour[contour.size() - 1]).norm();
	    if (dist_first < 1e-6) {
		    contour_complete = true;

		    closest_distance = dist_first;
		    contour.push_back(contour[0]);

	    } else {
            contour.push_back(connected_vertex);
            intersection_edges_plane.erase(intersection_edges_plane.begin() + closest_index, intersection_edges_plane.begin() + closest_index + 1);
	    }
	}

        /* Generate infill pattern
           Goal: Given the contours at each layer, this function aims to infill the internal part by a pattern which
                 can be procedurally built. (e.g. grid, honey comb, Fermat spiral) 
           The code is for grid pattern, you can rewrite the whole function based on your need. */
        void Infill(std::vector<std::vector<std::vector<Vector3<T>>>>& contours,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            
            infill_edges.clear();

            T min_x = 999999999.0;
            T max_x = -999999999.0;
            T min_y = 999999999.0;
            T max_y = -999999999.0;
            for (int i = 0; i < contours.size(); ++i) {
                std::vector<std::vector<Vector3<T>>> contour_plane = contours[i];
                for (int j = 0; j < contour_plane.size(); ++j ) {
                    std::vector<Vector3<T>> contour = contour_plane[j];
                    for (int k = 0; k < contour.size(); ++k) {
                        Vector3<T> vertex = contour[k];
                        if (vertex[0] < min_x) {
                            min_x = vertex[0];
                        }
                        if (vertex[0] > max_x) {
                            max_x = vertex[0];
                        }
                        if (vertex[1] < min_y) {
                            min_y = vertex[1];
                        }
                        if (vertex[1] > max_y) {
                            max_y = vertex[1];
                        }
                    }
                }
            }

            for (int i = 0;i < contours.size();++i) {
                std::vector<std::pair<Vector3<T>, Vector3<T>>> infill_edges_one_layer;
                infill_edges_one_layer.clear();

                /* Implement your code here */
                /* 1. find all intersections between contours and your infill pattern 
                2. infill internal space with desired pattern */
                T z = i * _dx;
                std::vector<std::vector<Vector3<T>>> contour_plane = contours[i];
                std::vector<std::vector<T>> x_intersections_plane;
                std::vector<std::vector<T>> y_intersections_plane;
                for (int j = 0; j < contour_plane.size(); ++j) {
                    std::vector<Vector3<T>> contour = contour_plane[j];
                    T x_index = 0;
                    for (T x = min_x; x <= max_x; x += _infill_dx) {
                        x_intersections_plane.push_back(std::vector<T>());
                        Vector3<T> prev_vertex = contour[0];
                        for (T k = 0; k < contour.size(); ++k) {
                            Vector3<T> curr_vertex = contour[k];
                            T prev_dist = x - prev_vertex[0];
                            T curr_dist = x - curr_vertex[0];
                            if ((abs(curr_vertex[0] - x) < 1e-6) && (k != 0)) {
                                x_intersections_plane[x_index].push_back(curr_vertex[1]);
                                continue;
                            }
                            if (prev_dist * curr_dist < 0) {
                                if (prev_dist < 0) {
                                    prev_dist = -1 * prev_dist;
                                }
                                if (curr_dist < 0) {
                                    curr_dist = -1 * (curr_dist);
                                }

                                Vector3<T> edge = curr_vertex - prev_vertex;
                                Vector3<T> edge_scaled = prev_vertex + edge*(prev_dist/(prev_dist + curr_dist));

                                x_intersections_plane[x_index].push_back(edge_scaled[1]);
                            }
                            prev_vertex = curr_vertex;
                        }
                        x_index += 1;
                    }

                    T y_index = 0;
                    for (T y = min_y; y <= max_y; y += _infill_dx) {
                        y_intersections_plane.push_back(std::vector<T>());
                        Vector3<T> prev_vertex = contour[0];
                        for (T k = 0; k < contour.size(); ++k) {
                            Vector3<T> curr_vertex = contour[k];
                            T prev_dist = y - prev_vertex[1];
                            T curr_dist = y - curr_vertex[1];
                            if ((abs(curr_vertex[1] - y) < 1e-6) && (k != 0)) {
                                y_intersections_plane[y_index].push_back(curr_vertex[0]);
                                continue;
                            }
                            if (prev_dist * curr_dist < 0) {
                                if (prev_dist < 0) {
                                    prev_dist = -1 * prev_dist;
                                }
                                if (curr_dist < 0) {
                                    curr_dist = -1 * curr_dist;
                                }
                                Vector3<T> edge = curr_vertex - prev_vertex;
                                Vector3<T> edge_scaled = prev_vertex + edge*(prev_dist/(prev_dist + curr_dist));
                                y_intersections_plane[y_index].push_back(edge_scaled[0]);
                            }
                            prev_vertex = curr_vertex;
                        }
                        y_index += 1;
                    }

                }

                for (int k = 0; k < x_intersections_plane.size(); ++k) {
                    std::sort(x_intersections_plane[k].begin(), x_intersections_plane[k].end());
                    T x = k * _infill_dx + min_x;
                    for (int l = 0; l < x_intersections_plane[k].size(); l += 2) {
                        Vector3<T> start = Vector3<T>(x, x_intersections_plane[k][l], z);
                        Vector3<T> stop = Vector3<T>(x, x_intersections_plane[k][l + 1], z);

                        infill_edges_one_layer.push_back(std::make_pair(start, stop));
                    }
                }
                for (int k = 0; k < y_intersections_plane.size(); ++k) {
                    std::sort(y_intersections_plane[k].begin(), y_intersections_plane[k].end());
                    T y = k * _infill_dx + min_y;
                    for (int l = 0; l < y_intersections_plane[k].size(); l += 2) {
                        Vector3<T> start = Vector3<T>(y_intersections_plane[k][l], y, z);
                        Vector3<T> stop = Vector3<T>(y_intersections_plane[k][l + 1], y, z);

                        infill_edges_one_layer.push_back(std::make_pair(start, stop));
                    }
                }

                infill_edges.push_back(infill_edges_one_layer);
            }
        }  

        /* Point cloud visualization for each layer's slicing
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            intersection_edges: edge soup for each layer's slicing */
        void VisualizeSlicing(std::string file_name, 
            T point_density,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> intersection_edges) {
            
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < intersection_edges.size();++i)
                for (int j = 0;j < intersection_edges[i].size();++j) {
                    Vector3<T> s_pos = intersection_edges[i][j].first;
                    Vector3<T> t_pos = intersection_edges[i][j].second;
                    int num_steps = (int)((t_pos - s_pos).norm() / point_density) + 1;
                    for (int step = 0;step <= num_steps;++step) {
                        Vector3<T> pos = s_pos * ((T)step / num_steps) + t_pos * ((T)1.0 - (T)step / num_steps);
                        points.push_back(pos);
                    }
                }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

        /* Point cloud visualization for each layer's contour
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            contour: each layer's contour list, each contour is a list a point in loop order */
        void VisualizeContour(std::string file_name,
            T point_density, 
            std::vector<std::vector<std::vector<Vector3<T>>>>& contour) {
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < contour.size();++i)
                for (int j = 0;j < contour[i].size();++j) 
                    for (int k = 0;k < contour[i][j].size();++k) {
                        Vector3<T> s_pos = contour[i][j][k];
                        Vector3<T> t_pos = contour[i][j][(k + 1) % contour[i][j].size()];
                        int num_steps = (int)((t_pos - s_pos).norm() / point_density) + 1;
                        for (int step = 0;step <= num_steps;++step) {
                            Vector3<T> pos = s_pos * ((T)step / num_steps) + t_pos * ((T)1.0 - (T)step / num_steps);
                            points.push_back(pos);
                        }
                    }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

        /* Point cloud visualization for each layer's slicing
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            infill_edges: edge soup for each layer's slicing */
        void VisualizeInfill(std::string file_name,
            T point_density,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < infill_edges.size();++i)
                for (int j = 0;j < infill_edges[i].size();++j) {
                    int num_steps = (int)((infill_edges[i][j].first - infill_edges[i][j].second).norm() / point_density) + 1;
                    for (int k = 0;k <= num_steps;++k)
                        points.push_back(infill_edges[i][j].first + (infill_edges[i][j].second - infill_edges[i][j].first) * (T)k / (T)num_steps);
                }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

    private:
        mesh::TriMesh<T> _tri_mesh;

        /* Variables for slicing */
        T _bottom, _top, _dx;

        /* Variables for infill algorithm */
        T _infill_dx;                                   // infill pattern will be equal-length grid
        T _infill_x_lower_bound, _infill_x_upper_bound;
        T _infill_y_lower_bound, _infill_y_upper_bound;

        /* accelerated data structure */
        data_structure::IntervalTree<T> _interval_tree;

        std::vector<std::vector<std::vector<Vector3<T>>>> z_map;
    };
}
