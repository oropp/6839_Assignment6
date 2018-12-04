#pragma once
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <iostream>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;
    template <typename T>
    using Vector2 = Eigen::Matrix<T, 2, 1>;
    template <typename T>
    using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    /*
     * A method to see if the new point (p3), requires turning counterclockwise from the first 2 points.
     */
    template <typename T>
    T CCW(Vector2<T> p1, Vector2<T> p2, Vector2<T> p3) {
        return (p2[0] - p1[0])*(p3[1] - p1[1]) - (p2[1] - p1[1])*(p3[0] - p1[0]);
    }


    // Q1: implement Graham Scan for 2D convex hull
    // The input is a vector of 2D points
    // The output should be the 2D points on the convex hull
    // Remove duplicates and sort your return vector
    // for comparing with the reference result
    template <typename T>
    std::vector<Vector2<T>> ConvexHull2D(const std::vector<Vector2<T>> &points) {
        // copy the points passed in so that the order can be changed
        std::vector<Vector2<T>> copiedPoints = points;

        // the index of the point with the lowest y value, starts at 0
        int lowest_y_index = 0;

        // Loop through all of the points and find the point with the lowest y value
        for (int i = 0; i < copiedPoints.size(); ++i) {
            T y_coord = copiedPoints[i][1];
            if (y_coord < copiedPoints[lowest_y_index][1]) {
                lowest_y_index = i;
            }
        }

        // Swap the point with the lowest y index with the first point in the list.
        std::swap(copiedPoints[lowest_y_index], copiedPoints[0]);

        // Store the first point in the list so it can be used in the sort.
        Vector2<T> point0 = copiedPoints[0];

        // Sort the points by their polar angle as related to the first point (point0). This first point should remain
        // first, so the sort starts at the second point.
        std::sort(copiedPoints.begin() + 1, copiedPoints.end(), [point0](const Vector2<T> &i, const Vector2<T> &j) {
                return atan2((i[1] - point0[1]), (i[0] - point0[0])) < atan2((j[1] - point0[1]), (j[0] - point0[0]));
        });

        // Set up the stack as a an empty vector
        std::vector<Vector2<T>> pile;
        // add the first three points to the stack
        pile.push_back(copiedPoints[0]);
        pile.push_back(copiedPoints[1]);
        pile.push_back(copiedPoints[2]);

        // For each of the other points
        for (int i = 3; i < copiedPoints.size(); ++i) {
            // calculate the counterclockwise
            double ccw = CCW(pile[pile.size() - 2], pile[pile.size() - 1], copiedPoints[i]);
            // If the value is negative, then it is counterclockwise, and the previous point was wrong.
            while (ccw <= 0.0) {
                // Remove the previous point from the stack
                pile.pop_back();
                // calculate the counter clockwise again
                ccw = CCW(pile[pile.size() - 2], pile[pile.size() - 1], copiedPoints[i]);
            }
            // Now that it is clockwise, add the point to the stack
            pile.push_back(copiedPoints[i]);
        }

        // Sort the stack by the x coordinate so that it can be compared to the results
        std::sort(pile.begin(), pile.end(), [](const Vector2<T> &i, const Vector2<T> &j) {
            return i[0] < j[0];
        });

        // return the stack
        return pile;
    }


    // Q2: implement brute-force method for Nd Pareto front
    // The input is a vector of Nd points
    // The output should be the Nd points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<VectorX<T>> ParetoFrontNdNaive(const std::vector<VectorX<T>> &points) {
        // Get the dimensions of the points that are passed in
        int dimension = points[0].size();

        // Set up an empty vector that will result in the pareto front
        std::vector<VectorX<T>> paretoFirstRule;

        // Loop through all the points
        for (int i = 0; i < points.size(); ++i) {
            // assume the point is on the pareto front
            bool isPareto = true;

            // Loop through all of the points again
            for (int j = 0; j < points.size(); ++j) {
                // If the two points are the same point, it is irrelevant
                if (i == j) {
                    continue;
                }

                // Assume that the first point i is worse than the second point j
                bool isWorse = true;

                // Loop through all of the dimensions of the point
                for (int k = 0; k < dimension; ++k) {
                    // If the value of point i is less than point j in this dimension, then i is not worse than j
                    if (points[i][k] < points[j][k]) {
                        isWorse = false;
                    }
                }

                // If point i is worse than point j, then it is not on the pareto front
                if (isWorse) {
                    isPareto = false;
                }
            }

            // If point i is on the pareto front, then add it to the pareto vector
            if (isPareto) {
                paretoFirstRule.push_back(points[i]);
            }
        }

        // Sort the points on the pareto front by the first dimension
        std::sort(paretoFirstRule.begin(), paretoFirstRule.end(), [](const VectorX<T> &i, const VectorX<T> &j) {
            return i[0] < j[0];
        });

        // Return the pareto front
        return paretoFirstRule;
    }

    // Q3: implement nlog(n) method for 2d Pareto front
    // The input is a vector of 2d points
    // The output should be the 2d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector2<T>> ParetoFront2D(std::vector<Vector2<T>> &points) {
        // Sort all of the points by the y value
        std::sort(points.begin(), points.end(), [](const Vector2<T> &point1, const Vector2<T> &point2) {
            return point1[0] > point2[0];
        });

        // Initialize an empty vector to store the pareto front
        std::vector<Vector2<T>> pareto;

        // Loop through all of the points
        for (int i = 0; i < points.size(); ++i) {
            // If it's the first point, then add it to the pareto front
            if (i == 0) {
                pareto.push_back(points[i]);
            } else {
                // the most recent point added to the pareto front
                Vector2<T> lastPoint = pareto[pareto.size() - 1];
                // the current point we're looking at
                Vector2<T> currentPoint = points[i];

                // while the y value of the current point is smaller than the y value of the last point on the pareto
                // front, remove that last point on the pareto front.
                while (currentPoint[1] <= lastPoint[1]) {
                    pareto.pop_back();
                    lastPoint = pareto[pareto.size() - 1];
                }

                // add the current point onto the pareto front
                pareto.push_back(points[i]);
            }
        }

        // Sort all the points on the pareto front by their x coordinate so result can be compared to given result
        std::sort(pareto.begin(), pareto.end(), [](const Vector2<T> &i, const Vector2<T> &j) {
            return i[0] < j[0];
        });

        // return the points on the pareto front
        return pareto;
    }



    // bonus question: implement nlog(n) method for 3d Pareto front
    // The input is a vector of 3d points
    // The output should be the 3d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector3<T>> ParetoFront3D(const std::vector<Vector3<T>> &points) {
        return points;
    }
}
