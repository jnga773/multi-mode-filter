!----------------------------------!
! < b^{\dagger}_{j} a_{k} \sigma > !
!----------------------------------!
! Set the diagonal matrix elements for M
Mat = Mat_OG
DO x = 1, N_mat
  Mat(x, x) = Mat(x, x) - (kappa(a) + kappa(b) - i * (wl(j, b) - wl(k, a)))
END DO

! Set the non-homogeneous vector
B_vec = 0.0d0
B_vec(1) = -CONJG(gkl(j, b)) * cavsig2_ss(k, a, eg, 1) + &
         & -gkl(k, a) * cavsig2_ss(j, at, ge, 2)
B_vec(2) = -CONJG(gkl(j, b)) * cavsig2_ss(k, a, ee, 1) + &
         & -gkl(k, a) * xi * cavsig2_ss(j, at, gf, 2)
B_vec(3) = -CONJG(gkl(j, b)) * xi * cavsig2_ss(k, a, fg, 1) + &
         & -gkl(k, a) * cavsig2_ss(j, at, ee, 2)
B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, ata, 2, 1) + &
         & -CONJG(gkl(j, b)) * xi * cavsig2_ss(k, a, fe, 1) + &
         & -gkl(k, a) * xi * cavsig2_ss(j, at, ef, 2)
B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, ata, 2, 1) + &
         & CONJG(gkl(j, b)) * xi * cavsig2_ss(k, a, gg, 1) + &
         & CONJG(gkl(j, b)) * xi * cavsig2_ss(k, a, ee, 1) + &
         & -CONJG(gkl(j, b)) * xi * cav1_ss(k, a, 1)
B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, ata, 2, 1) + &
         & gkl(k, a) * xi * cavsig2_ss(j, at, gg, 2) + &
         & gkl(k, a) * xi * cavsig2_ss(j, at, ee, 2) + &
         & -gkl(k, a) * xi * cav1_ss(j, at, 2)
B_vec(7) = -CONJG(gkl(j, b)) * cavsig2_ss(k, a, ef, 1)
B_vec(8) = - gkl(k, a) * cavsig2_ss(j, at, fe, 2)

! Set inverse matrix
Mat_inv = Mat
! Perform LU-factorization of matrix
CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
! Invert Matrix (Optimal LWORK = 8)
LWORK = 8
CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

! Calculate steady state
cavsig3_ss(j, k, bta, :) = -MATMUL(Mat_inv, B_vec)
