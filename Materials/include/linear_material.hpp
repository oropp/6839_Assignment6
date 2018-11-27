// Tao Du
// taodu@csail.mit.edu
// Sept 22, 2016
#pragma once
#include "material.hpp"

namespace materials {

    template<int dim, typename T>
    class LinearElasticityMaterial : public Material<dim, T> {
    public:
        LinearElasticityMaterial(const T young_modulus,
                                 const T poisson_ratio) : Material<dim, T>(young_modulus, poisson_ratio) {}

        LinearElasticityMaterial(const LinearElasticityMaterial<dim, T>& material)  : Material<dim, T>(material) {}

        ~LinearElasticityMaterial() {}


        const T EnergyDensity(const typename Material<dim, T>::MatrixDimT& F) const {
            const typename Material<dim, T>::MatrixDimT small_strain_tensor =
                    0.5 * (F + F.transpose()) - Material<dim, T>::MatrixDimT::Identity();
            const T trace = small_strain_tensor.trace();
            return Material<dim, T>::mu() * small_strain_tensor.array().square().sum()
                   + Material<dim, T>::lambda() / 2.0 * trace * trace;
        }



        const typename Material<dim, T>::MatrixDimT StressTensor(
                const typename Material<dim, T>::MatrixDimT& F) const {
            const typename Material<dim, T>::MatrixDimT I = Material<dim, T>::MatrixDimT::Identity();
            return Material<dim, T>::mu() * (F + F.transpose() - 2 * I) +
                   Material<dim, T>::lambda() * (F.trace() - dim) * I;
        }


        const typename Material<dim, T>::MatrixDimT StressDifferential(
                const typename Material<dim, T>::MatrixDimT& F,
                const typename Material<dim, T>::MatrixDimT& dF) const {

	    /*Note: This is wrong.  normally, this would return the differential dP/dF * dF.  
	      We have deleted the code (since it would give away the answer to the other StressDifferential).  
              You don't need this function to complete the assignment.
            */
            const typename Material<dim, T>::MatrixDimT I = Material<dim, T>::MatrixDimT::Identity();
            return I;

        };

         //TODO: Students fill this out (and change the return value)
         //F is the deformation gradient.
        const typename Material<dim, T>::MatrixDim2T StressDifferential(
                const typename Material<dim, T>::MatrixDimT& F) const {

            const typename Material<dim, T>::MatrixDim2T zero = Material<dim, T>::MatrixDim2T::Zero();

            typename Material<dim, T>::MatrixDim2T diff_p;

            // Loop through all rows of P
            for (int i = 0; i < 3; ++i) {
                // Loop through all columns of P
                for (int j = 0; j < 3; ++j) {
                    // Loop through all rows of F
                    for (int k = 0; k < 3; ++k) {
                        // Loop through all columns of F
                        for (int l = 0; l < 3; ++l) {
                            /*
                             * Case 1:
                             * the element in P is the transpose location of the element in F.
                             * i = l and j = k
                             */
                            if (i == l && j == k) {
                                diff_p(3*i + j, 3*k + l) = Material<dim, T>::mu();
                            }
                            /*
                             * Case 2:
                             * The differential of an element of P along the diagonal with respect to an element of F
                             * along the diagonal in the same location.
                             * i = j = k = l
                             */
                            else if (i == j == k == l) {
                                diff_p(3*i + j, 3*k + l) = 2*Material<dim, T>::mu() + Material<dim, T>::lambda();
                            }
                            /*
                             * Case 3:
                             * The differential of an element of P along the diagonal with respect to an element of F
                             * along the diagonal, but they are not at the same location.
                             * i = j and k = l, but i /= k
                             */
                            else if (i ==j && k == l) {
                                diff_p(3*i + j, 3*k + l) = Material<dim, T>::lambda();
                            }
                            /*
                             * Case 4:
                             * The differential of an element of P with respect to the same element in F.
                             * i = k and j = l, but none of the previous cases apply
                             */
                            else if (i == k && j == l) {
                                diff_p(3*i + j, 3*k + l) = Material<dim, T>::mu();
                            }
                            /*
                             * Case 5:
                             * If none of the previous cases apply, it falls into this case.
                             */
                            else {
                                diff_p(3*i + j, 3*k + l) = 0;
                            }
                        }
                    }
                }
            }

            return diff_p;

        }

    private:
        LinearElasticityMaterial<dim, T>& operator=(
                const LinearElasticityMaterial<dim, T>&);



    };

}
