// Reference Implementation by: Alexandre Kaspar <akaspar@mit.edu>
#pragma once

#include <algorithm>
#include <cassert>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>
#include <igl/parallel_for.h>
#include <igl/polar_svd3x3.h>

// edge indices
typedef typename Eigen::Matrix<int, 3, 2> Matrix32i;
const Matrix32i edges = (Matrix32i() << 1, 2, 2, 0, 0, 1).finished();

/**
 * Problem to minimize
 *
 *   sum_k sum_{ij in F(k)} c_ij | (v_i - v_j) - R_k (v_i' - v_j')|^2
 *
 */

/**
 * Compute entries of the cotangent matrix C
 * (c_ij) for i in F, and j in {0,1,2} = triangle edge index
 *
 * @param V (nx3) input vertices
 * @param F (mx3) input triangle faces
 * @param C (mx3) output cotangents
 */
template <typename TypeV, typename TypeF, typename TypeC>
void cotmatrix_entries(
  const TypeV & V,
  const TypeF & F,
  TypeC & C)
{
  typedef typename TypeV::Scalar Scalar;

  const int m = F.rows();
  C.resize(m, 3);
  // compute per face (m of them)
  igl::parallel_for( m, [&V, &F, &C](const int i){
    // squared edge lengths
    Scalar l[3], e[3];
    for(int j = 0; j < 3; ++j){
      const auto &v0 = V.row(F(i, edges(j, 0)));
      const auto &v1 = V.row(F(i, edges(j, 1)));
      Scalar l2 = (v0 - v1).squaredNorm();
      l[j] = std::sqrt(l2);
      e[j] = l2;
    }
    // face area using Heron's formula
    // @see https://en.wikipedia.org/wiki/Heron%27s_formula
    // improved stability using Karan's formula
    // @see Miscalculating Area and Angles of a Needle-like Triangle, Kahan 2014
    /* unstable version:
    Scalar s = (l[0] + l[1] + l[2]) * 0.5;
    Scalar dblA = s * (s - l[0]) * (s - l[1]) * (s - l[2]);
    */
    std::sort(l, l + 3); // ensures l[0] <= l[1] <= l[2]
    // stable version:
    Scalar arg  = (l[0] + (l[1] + l[2]))
                * (l[2] - (l[0] - l[1]))
                * (l[2] + (l[0] - l[1]))
                * (l[0] + (l[1] - l[2]));
    Scalar dblA = 2.0 * 0.25 * std::sqrt(arg);
    if(dblA != dblA){
      // NaN
      dblA = 0.0; // C(i, j) = infinity
    }
    
    // cotangent and diagonal entries
    // @see https://www.cs.toronto.edu/~jacobson/images/alec-jacobson-thesis-2013-compressed.pdf
    // @see 2.1.2 Cotangent weights
    // @see Appendix 4.8: A cotangent Laplacian for images as surfaces
    for(int j = 0; j < 3; ++j){
      Scalar e0 = e[edges(j, 0)];
      Scalar e1 = e[edges(j, 1)];
      Scalar e2 = e[3 - edges(j, 0) - edges(j, 1)];
      C(i, j) = (e0 + e1 - e2) / dblA / 4.0;
    }

  }, 1000);
}

/**
 * Compute cotangent laplacian sparse matrix L
 * (l_ii) = 1/2 sum_ik (cot \alpha_k + cot \beta_k)
 * (l_ij) = -1/2 (cot \alpha_j + cot \beta_j)
 * where ij and ik are edges.
 *
 * @param V (nx3) input vertices
 * @param F (mx3) input triangle faces
 * @param C (mx3) input cotangent weights
 * @param L (nxn) output laplacian matrix
 * @see https://graphics.stanford.edu/courses/cs468-13-spring/assets/lecture12-lu.pdf
 */
template <typename TypeV, typename TypeF, typename TypeC, typename TypeL>
void cotmatrix(
  const TypeV & V,
  const TypeF & F,
  const TypeC & C,
  TypeL & L)
{
  typedef typename TypeV::Scalar Scalar;
  typedef typename Eigen::Triplet<Scalar> Entry;

  // allocate laplacian
  const int n = V.rows();
  const int m = F.rows();
  L.resize(n, n);
  L.reserve(10 * n);

  // create sparse matrix entries
  std::vector<Entry> entries;
  entries.reserve(m * 3 * 4);

  /*
   * After discussions with Andrew Spielberg, Alexandre Kaspar, and Andy Wang, I changed
   * my approach from looping over all of the vertices to looping over all of the faces.
   *
   * I believe that this approach will still find all of the weights necessary, because
   * the weights not on the faces should be zero.
   */

  // for all faces
  for (int f = 0; f < m; ++f) {
      // for all half edges
      for (int e = 0; e < 3; ++e) {
          // Get the two vertices at either end of the edge
          auto v0 = edges(e, 0);
          auto v1 = edges(e, 1);
          int i = F(f, v0);
          int j = F(f, v1);

          // Get the corresponding cotangent weight for this edge.
          double val = C(f, e);

          // Add all of the combinations of i and j to entries, making sure that the non-diagonal values are negative
          entries.push_back(Eigen::Triplet<Scalar>(i, j, -val));
          entries.push_back(Eigen::Triplet<Scalar>(j, i, -val));
          entries.push_back(Eigen::Triplet<Scalar>(i, i, val));
          entries.push_back(Eigen::Triplet<Scalar>(j, j, val));
      }
  }
  L.setFromTriplets(entries.begin(), entries.end());
}

/**
 * Create the entries of a block of the ARAP K matrix
 *
 * @param V (nx3) input vertices
 * @param F (mx3) input faces
 * @param C (nxn) input cotangent weight matrix
 * @param d (1)   input dimension (0, 1, or 2) corresponding to (x, y or z)
 * @param fromRow input starting entry row offset
 * @param fromCol input starting entry column offset
 * @param entries output list of triplets to which new entries are appended
 * @param loop    input whether to use the loop implementation from IGL (slower)
 */
template <typename MatV, typename MatF, typename MatC, typename Scalar>
void arap_linear_block(
  const MatV & V,
  const MatF & F,
  const MatC & C,
  const int d,
  const int fromRow,
  const int fromCol,
  std::vector<Eigen::Triplet<Scalar> > &entries,
  bool loop = false)
{
  typedef typename Eigen::Triplet<Scalar> Entry;

  // spokes+rims energy implementation
  const int n = V.rows();
  const int m = F.rows(); // number of faces
  
  // workspace
  entries.reserve(entries.size() + 7 * n);

  // for all faces
  for(int i = 0; i < m; ++i){
    // for each half-edge
    for(int e = 0; e < 3; ++e){
        // Get each of the three vertices from the face
      int src = F(i, edges(e, 0));
      int trg = F(i, edges(e, 1));
      int oth = F(i, 3 - edges(e, 0) - edges(e, 1));

      // Calculate the value for this edge from the C matrix
      auto delta = V.row(src) - V.row(trg);
      double val = C(i, e) * delta(d)/3.0;

      // This is the offset to account for the different sub matrices of K, Kx, Ky, and Kz
      int offset = d*n;

      // Add all of the combinations of the vertices to entries, accoutning for the correct sign of value.
      // Don't need oth, oth because these will be considered in other values of e.
      entries.push_back(Eigen::Triplet<Scalar>(src, trg + offset, val));
      entries.push_back(Eigen::Triplet<Scalar>(trg, src + offset, -val));
      entries.push_back(Eigen::Triplet<Scalar>(src, src + offset, val));
      entries.push_back(Eigen::Triplet<Scalar>(trg, trg + offset, -val));
      entries.push_back(Eigen::Triplet<Scalar>(src, oth + offset, val));
      entries.push_back(Eigen::Triplet<Scalar>(trg, oth + offset, -val));
    }
  }
}


/**
 * Precomputations for the ARAP quadratic problem using min_quad_with_fixed.
 *
 * min_quad_with_fixed minimizes
 *
 *    trace( 0.5*Z'*A*Z + Z'*B + constant )
 *
 * subject to
 *
 *    Z(known,:) = Y, and
 *    Aeq*Z = Beq
 * 
 * The call is:
 * 
 *    min_quad_with_fixed_precompute(A, known, Aeq, pd, data)
 * 
 * where pd = simplification when A is positive definite
 *
 *
 * @param V (nx3)   input vertices of rest shape
 * @param F (mx3)   input faces pointing to vertices
 * @param b (bx1)   input indices of handle vertices
 * @param data      output datastructure for min_quad precomputation
 * @param K (nx3n)  output K matrix
 */
template <typename TypeV, typename TypeF, typename Typeb, typename Scalar>
void arap_precompute(
  const TypeV & V,
  const TypeF & F,
  const Typeb & b,
  igl::min_quad_with_fixed_data<Scalar> & data,
  Eigen::SparseMatrix<Scalar> & K)
{
  using namespace Eigen;

  const int n = V.rows();
  // const int m = F.rows();

  // compute cotangent weight matrix c_ij
  Matrix<Scalar, Dynamic, Dynamic> C;
  cotmatrix_entries(V, F, C);

  // get Laplacian matrix
  SparseMatrix<Scalar> L;
  cotmatrix(V, F, C, L);

  // compute K = [Kx, Ky, Kz]
  std::vector<Triplet<Scalar> > entries;

  // For all values of d
  for (int d = 0; d < 3; ++d) {
      int fromRow = 0;
      // Fill th entries with the values for one chunk of the K matrix
      arap_linear_block(V, F, C, d, fromRow*n, d*n, entries);
  }
  // Transform the entries into the K matrix
  K.resize(n, 3*n);
  K.setFromTriplets(entries.begin(), entries.end());

  // Precompute the min quad so that we can iterate more efficiently
  igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<Scalar>(), true, data);
}

/**
 * Find best rotation matrix R that minimizes
 * the following linear trace problem:
 *
 *   trace(C * R)
 *
 * where C = V^T K = [X Y Z]^T K = [Cx Cy Cz]^T
 * with  V in |R^{nx3}  and K in |R^{nx3n}
 *       C in |R^{3x3n} and R in |R^{3nx3} = [Rx Ry Rz]^T
 *
 * Notes:
 *   Rx, Ry, Rz in |R^{nx3}
 * with
 *   Rd = [Rd_1^T Rd_2^T ... Rd_n^T]^T
 * where
 *   Rd_i in |R^{1x3}
 * which can be assembled into individual rotations
 *   R_i = [Rx_i^T Ry_i^T Rz_i^T]^T in |R^{3x3}
 *
 * The sub-rotations R_i are all disconnected.
 * Thus they are solved in parallel independently using SVD.
 *
 * @param C (3x3n) input coefficient matrix encapsulating V and K
 * @param R (3nx3) output rotation matrix
 */
template<typename MatS, typename MatR>
void solve_rotations(
  const Eigen::PlainObjectBase<MatS> & C,
  Eigen::PlainObjectBase<MatR> & R)
{
  typedef typename MatR::Scalar Scalar;
  typedef Eigen::Matrix<Scalar, 3, 3> Mat3;
  
  const int dim = C.rows();
  assert(dim == 3 && "Only supporting 3D rotations");
  const int nr = C.cols() / dim;
  assert(nr * dim == C.cols() && "Dimension of C must use multiple of dimension");

  R.resize(3 * nr, 3);
  // for each rotation

  // For all vertices
  for (int i = 0; i < nr; i += 3) {
      Mat3 C_i, R_i;

      // Fill the C_i matrix from the C matrix that's passed in
      for (int j = 0; j < 3; ++j) {
          C_i(0, j) = C(j, i);
          C_i(1, j) = C(j, i + nr);
          C_i(2, j) = C(j, i + 2*nr);
      }


      // Use the SVD to get the R_i matrix from C_i
      igl::polar_svd3x3(C_i, R_i);

      // Decompose the R_i matrix into the larger R matrix
      R_i.transpose();
      R.row(i) = R_i.row(0);
      R.row(i + nr) = R_i.row(1);
      R.row(i + 2*nr) = R_i.row(2);
  }
}

template<typename Scalar, typename Typebc, typename TypeU>
void arap_single_iteration(
  const igl::min_quad_with_fixed_data<Scalar> & data,
  const Eigen::SparseMatrix<Scalar> & K,
  const Typebc & bc,
  TypeU & U)
{
  using namespace Eigen;
  
  // number of vertices
  const int n = U.rows();

  // enforce boundary conditions exactly
  for(int bi = 0; bi < bc.rows(); bi++)
    U.row(data.known(bi)) = bc.row(bi);

  // matrix of rotations
  Matrix<Scalar, Dynamic, Dynamic> R(n * 3, 3);

  // local solve: fix vertices, find rotations
local: {
    // Create the C matrix needed to pass into solve rotations.
    Matrix<Scalar, Dynamic, Dynamic> C(3, 3*n);
    C = U.transpose()*K;

    solve_rotations(C, R);
  }

  // global solve: fix rotations, find vertices
global: {
    typedef Matrix<Scalar, Dynamic, 1> Vector;

    // Create the B and Beq matrices needed to call igl::min_quad_with_fixed_solve
    auto B = -(K*R).eval();
    Matrix<Scalar, Dynamic, Dynamic> Beq;

    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
    }
}
