// Tao Du
// taodu@csail.mit.edu
// Oct 12, 2016
#pragma once
#include "deformable_body.hpp"
#include "hexahedral_mesh.hpp"
#include "typedefs.hpp"

namespace materials {

    template<typename T>
    class HexDeformableBody : public DeformableBody<3, T> {
    public:
        HexDeformableBody(const Material<3, T>& material,
                          const Matrix3X<T>& initial_vertex_position,
                          const T density, const HexahedralMesh<T>& undeformed_hex_mesh)
                          : DeformableBody<3, T>(material, initial_vertex_position, undeformed_hex_mesh),
                            hex_size_((undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(0)) - undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(1))).norm()) {}






        HexDeformableBody(const std::vector<std::reference_wrapper<const Material<3, T>>>& materials,
                          const std::vector<int>& material_id,
                          const Matrix3X<T>& initial_vertex_position,
                          const T density, const HexahedralMesh<T>& undeformed_hex_mesh)
                          : DeformableBody<3, T>(materials, material_id, initial_vertex_position, undeformed_hex_mesh),
                            hex_size_((undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(0)) - undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(1))).norm()) {}

        HexDeformableBody(const HexDeformableBody& hex_body) : DeformableBody<3, T>(hex_body),
                                                               hex_size_(hex_body.hex_size_) {}

        ~HexDeformableBody() {}


        //TODO: Students should fill this out
        const Eigen::SparseMatrix<T> ComputeStiffnessMatrix(
                const Matrix3X<T>& vertices) const{
            std::vector<Eigen::Triplet<T>> triplet_list;
            const int vertex_num = static_cast<int>(this->vertex_position_.size() / 3);
           
            Eigen::SparseMatrix<T> K(vertex_num * 3, vertex_num * 3);
            K.setFromTriplets(triplet_list.begin(), triplet_list.end());

            // Get the undeformed mesh
            PolyMesh<3, T> undeformed_mesh = this -> undeformed_mesh_;

            /*
             * Loop through all of the voxels in the mesh. The voxels are stored in the elements matrix of the
             * undeformed mesh.
             *
             * The elements matrix is an 8xn matrix where each column is a single voxel and n is the number of voxels
             */
            for (int voxel_index = 0; voxel_index < undeformed_mesh.NumOfElement(); ++voxel_index) {
                /*
                 * Because the voxels start as a cube, using one side of the undeformed voxel you can calculate the
                 * volume
                 */
                T volume = pow(undeformed_mesh.vertex(1)[1] - undeformed_mesh.vertex(0)[1], 3);

                // Initializing data structures
                Eigen::Matrix<T, 24, 24> k_voxel;
                Eigen::Matrix<T, 3, 8> undeformed_cube;
                Eigen::Matrix<T, 3, 8> deformed_cube;

                /*
                 * Loop through the 8 vertices in the voxel and fill the undeformed_cube and deformed_cube
                 * data structures
                 */
                for (int i = 0; i < 8; ++i) {
                    Eigen::Vector3d vertex = undeformed_mesh.vertex(undeformed_mesh.element(voxel_index)[i]);
                    undeformed_cube(i, 0) = vertex[0];
                    undeformed_cube(i, 1) = vertex[1];
                    undeformed_cube(i, 2) = vertex[2];

                    Eigen::Vector3d deformed_vertex = vertices.col(undeformed_mesh.element(voxel_index)[i]);
                    deformed_cube(i, 0) = deformed_vertex[0];
                    deformed_cube(i, 1) = deformed_vertex[1];
                    deformed_cube(i, 2) = deformed_vertex[2];
                }

                // Loop through all the vertices in the voxel
                for (int v = 0; v < 8; ++v) {
                    // Get the actual x, y, and z coordinates of the vertex we are currently on
                    Vector3<T> vertex = undeformed_mesh.vertex(undeformed_mesh.element(voxel_index)[v]);

                    // Get the dF/dx matrix using DeformationGradientPartialx
                    Eigen::Matrix<T, 9, 24> df_dx = DeformationGradientPartialx(vertex, undeformed_cube, deformed_cube);

                    // Get the dP/dF matrrix from StressDifferential. In order to do this we need the forces F, which we
                    // get from the material.
                    int material_index = this->material_id_[voxel_index];
                    auto material = this -> materials_[material_index];
                    Eigen::Matrix<T, 3, 3> F = DeformationGradient(vertex, undeformed_cube, deformed_cube);
                    Eigen::Matrix<T, 9, 9> dp_df = material.get().StressDifferential(F);

                    // Calculate the K matrix for this voxel and vertex combination, using the equation described in the
                    // write up.
                    Eigen::Matrix<T, 24, 24> k_vertex = (1/8)*(df_dx.transpose()*dp_df*df_dx)*volume;

                    /*
                     * Looping through all of the vertices in the voxel twice to get all the combinations in order to
                     * store the information in the larger K matrix
                     *
                     * For each ij index I store the 3x3 matrix in the correct indices in the larger K matrix
                     */
                    for (int i = 0; i < 8; ++i) {
                        int row_index = undeformed_mesh.element(voxel_index)[i];
                        for (int j = 0; j < 8; ++j) {
                            int col_index = undeformed_mesh.element(voxel_index)[j];
                            triplet_list.push_back(Eigen::Triplet<T>(row_index, col_index, k_vertex(3*i, 3*j)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index, col_index + 1, k_vertex(3*i, 3*j + 1)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index, col_index + 2, k_vertex(3*i, 3*j + 2)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index + 1, col_index, k_vertex(3*i + 1, 3*j)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index + 1, col_index + 1, k_vertex(3*i + 1, 3*j + 1)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index + 1, col_index + 2, k_vertex(3*i + 1, 3*j + 2)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index + 2, col_index, k_vertex(3*i + 2, 3*j)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index + 2, col_index + 1, k_vertex(3*i + 2, 3*j + 1)));
                            triplet_list.push_back(Eigen::Triplet<T>(row_index + 2, col_index + 2, k_vertex(3*i + 2, 3*j + 2)));
                        }
                    }
                }
            }

            // Make sure K is symmetric.
            K = (K + Eigen::SparseMatrix<T>(K.transpose())) / 2.0;
            return K;
        }

        //return dphi (the deformation gradient) for a given voxel:
        //undeformed_vertex is a point in material space.
        //undeformed_cube are the vertices of an undeformed voxel
        //deformed_cube are the vertices of a the deformed voxel.
        static const Eigen::Matrix<T, 3, 3> DeformationGradient(
                const Vector3<T>& undeformed_vertex,
                const Eigen::Matrix<T, 3, 8>& undeformed_cube,
                const Eigen::Matrix<T, 3, 8>& deformed_cube){
            // Rename variables.
            const Vector3<T>& X = undeformed_vertex;
            const Eigen::Matrix<T, 3, 8>& X0 = undeformed_cube;
            const Eigen::Matrix<T, 3, 8>& x0 = deformed_cube;
            const T dx = X0(0, 4) - X0(0, 0);
            const T inv_dx = 1.0 / dx;
            const Vector3<T> offset = X - X0.col(0);
            const T rx = offset.x() / dx;
            const T ry = offset.y() / dx;
            const T rz = offset.z() / dx;
            Eigen::Matrix<T, 3, 3> F = Eigen::Matrix<T, 3, 3>::Zero();
            const T x_factor[2] = {1 - rx, rx};
            const T y_factor[2] = {1 - ry, ry};
            const T z_factor[2] = {1 - rz, rz};
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        F.col(0) += x0.col(4 * i + 2 * j + k)
                                    * (i == 0 ? -inv_dx : inv_dx) * y_factor[j] * z_factor[k];
                        F.col(1) += x0.col(4 * i + 2 * j + k)
                                    * x_factor[i] * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        F.col(2) += x0.col(4 * i + 2 * j + k)
                                    * x_factor[i] * y_factor[j] * (k == 0 ? -inv_dx : inv_dx);
                    }
            return F;
        }




        //return dphi/dx for a given voxel:
        //undeformed_vertex is a point in material space.
        //undeformed_cube are the vertices of an undeformed voxel
        //deformed_cube are the vertices of a the deformed voxel.
        static const Eigen::Matrix<T, 9, 24> DeformationGradientPartialx(
                const Vector3<T>& undeformed_vertex,
                const Eigen::Matrix<T, 3, 8>& undeformed_cube,
                const Eigen::Matrix<T, 3, 8>& deformed_cube){

            const Vector3<T>& X = undeformed_vertex;
            const Eigen::Matrix<T, 3, 8>& X0 = undeformed_cube;
            Eigen::Matrix<T, 9, 24> Jacobian = MatrixX<T>::Zero(9, 24);
            const T dx = X0(0, 4) - X0(0, 0);
            const T inv_dx = 1.0 / dx;
            const Vector3<T> offset = X - X0.col(0);
            const T rx = offset.x() / dx;
            const T ry = offset.y() / dx;
            const T rz = offset.z() / dx;
            Eigen::Matrix<T, 3, 3> F = Eigen::Matrix<T, 3, 3>::Zero();
            const T x_factor[2] = { 1 - rx, rx };
            const T y_factor[2] = { 1 - ry, ry };
            const T z_factor[2] = { 1 - rz, rz };
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        const int index = 4 * i + 2 * j + k;
                        const T scale_first_column = (i == 0 ? -inv_dx : inv_dx)
                                                          * y_factor[j] * z_factor[k];
                        Jacobian(0, 3 * index) += scale_first_column;
                        Jacobian(1, 3 * index + 1) += scale_first_column;
                        Jacobian(2, 3 * index + 2) += scale_first_column;
                        const T scale_second_column = x_factor[i]
                                                           * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        Jacobian(3, 3 * index) += scale_second_column;
                        Jacobian(4, 3 * index + 1) += scale_second_column;
                        Jacobian(5, 3 * index + 2) += scale_second_column;
                        const T scale_third_column = x_factor[i] * y_factor[j]
                                                          * (k == 0 ? -inv_dx : inv_dx);
                        Jacobian(6, 3 * index) += scale_third_column;
                        Jacobian(7, 3 * index + 1) += scale_third_column;
                        Jacobian(8, 3 * index + 2) += scale_third_column;
                    }
            return Jacobian;



        }

        virtual void getInitialNodes(Matrix3X<T>& initial_nodes){
            initial_nodes = this->undeformed_mesh_.vertex();
        }

    private:
        HexDeformableBody& operator=(const HexDeformableBody&);


        //TODO: Studnets fill this function out
        static const Eigen::Matrix<T, 8, 8> GaussIntegrationFactor() {
            // \int_{-1}^{1} f(x) dx \approx f(-1/sqrt(3)) + f(1/sqrt(3)).
            Eigen::Matrix<T, 8, 8> X0_coeff = Eigen::MatrixXd::Zero(8, 8);


            std::vector<Vector3<T>> quadrature_pts;

            // The 8 quadrature points used to calculate
            quadrature_pts.push_back(Vector3<T>(-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(-1/sqrt(3), -1/sqrt(3), 1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(-1/sqrt(3), 1/sqrt(3), -1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(1/sqrt(3), -1/sqrt(3), -1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(1/sqrt(3), -1/sqrt(3), 1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(1/sqrt(3), 1/sqrt(3), -1/sqrt(3)));
            quadrature_pts.push_back(Vector3<T>(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)));

            // Get all combinations of i, j, and k that will give us the nodal rules
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int index = 4*i + 2*j + k;

                        // For each nodal rule, calculate the nodal value for each quadrature point
                        for (int e = 0; e < 8; ++e) {
                            Vector3<T> quad = quadrature_pts[e];
                            T N = (1/8)*((1 + pow(-1, i)*quad[0])*(1 + pow(-1, j)*quad[1])*(1 + pow(-1, k)*quad[2]));
                            X0_coeff(e, index) = N;
                        }

                    }
                }
            }
            return X0_coeff;
        }

        const T hex_size_;
    };

}
