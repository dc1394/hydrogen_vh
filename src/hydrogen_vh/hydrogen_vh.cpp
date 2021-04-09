/*! \file hydrogen_vh.cpp
    \brief 水素原子のHartreeポテンシャルを求めるクラスの実装
    Copyright © 2020 @dc1394 All Rights Reserved.
    (but this is originally adapted by sunsetyuhi for fem1d_poisson.py from https://github.com/sunsetyuhi/fem_py/blob/master/fem1d_poisson )

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "hydrogen_vh.h"
#include <cmath>                        // for std::pow, std::sqrt
#include <cstdio>                       // for FILE, std::fclose, std::fopen, std::fprintf
#include <memory>                       // for std::unique_ptr
#include <boost/assert.hpp>             // for boost::assert        
#include <Eigen/SparseLU>               // for Eigen::SparseLU

namespace hydrogen_vh {
    // #region コンストラクタ・デストラクタ

    using namespace hydrogen_fem;

    Hydrogen_Vh::Hydrogen_Vh()
        :   gl_(INTEGTABLENUM),
            mat_A_ele_(boost::extents[hydrogen_fem::Hydrogen_FEM::ELE_TOTAL][2][2]),
            mat_A_glo_(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL, hydrogen_fem::Hydrogen_FEM::NODE_TOTAL),
            u_(Eigen::VectorXd::Zero(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL)),
            vec_b_ele_(boost::extents[hydrogen_fem::Hydrogen_FEM::ELE_TOTAL][2]),
            vec_b_glo_(Eigen::VectorXd::Zero(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL))
    {
        hfem_.do_run();
    }

    // #endregion コンストラクタ・デストラクタ

    // #region publicメンバ関数 

    void Hydrogen_Vh::do_run()
    {
        // 要素行列とLocal節点ベクトルを生成
        make_element_matrix_and_vector();

        // 全体行列と全体ベクトルを生成
        make_global_matrix_and_vector();

        // 境界条件処理を行う
        boundary_conditions();

        // 連立方程式を解く
        Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        solver.analyzePattern(mat_A_glo_);
        solver.factorize(mat_A_glo_);

        u_ = solver.solve(vec_b_glo_);
    }

    void Hydrogen_Vh::save_result() const
    {
        std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(Hydrogen_Vh::RESULT_FILENAME, "w"), std::fclose);

        auto const dr = (Hydrogen_FEM::R_MAX - Hydrogen_FEM::R_MIN) / static_cast<double>(hydrogen_fem::Hydrogen_FEM::ELE_TOTAL);
        for (auto i = 1; i < hydrogen_fem::Hydrogen_FEM::NODE_TOTAL; i++) {
            auto const r = static_cast<double>(i) * dr;
            // 厳密な結果と比較
            std::fprintf(fp.get(), "%.14f, %.14f, %.14f\n", r, u_[i] / r, - (1.0 + 1.0 / r) * exp(-2.0 * r) + 1.0 / r);
        }
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void Hydrogen_Vh::boundary_conditions()
    {
        auto const a = 0.0;
        mat_A_glo_.coeffRef(0, 0) = 1.0;
        vec_b_glo_(0) = a;
        vec_b_glo_(1) -= a * mat_A_glo_.coeff(0, 1);
        mat_A_glo_.coeffRef(0, 1) = 0.0;
        mat_A_glo_.coeffRef(1, 0) = 0.0;

        auto const b = 1.0;
        mat_A_glo_.coeffRef(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 1, hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 1) = 1.0;
        vec_b_glo_(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 1) = b;
        vec_b_glo_(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 2) -= b * mat_A_glo_.coeff(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 2, hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 1);
        mat_A_glo_.coeffRef(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 2, hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 1) = 0.0;
        mat_A_glo_.coeffRef(hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 1, hydrogen_fem::Hydrogen_FEM::NODE_TOTAL - 2) = 0.0;
    }

    void Hydrogen_Vh::make_element_matrix_and_vector()
    {
        // 要素行列とLocal節点ベクトルの各成分を計算
        for (auto e = 0; e < hydrogen_fem::Hydrogen_FEM::ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    mat_A_ele_[e][i][j] = (std::pow(-1, i + 1) * std::pow((-1), (j + 1)) / hfem_.Length()[e]);
                }

                switch (i)
                {
                case 0:
                    vec_b_ele_[e][i] = gl_.qgauss([this, e](double r) { return r * hfem_.rho(r) * (hfem_.Node_r_ele()[e][1] - r) / hfem_.Length()[e]; },
                        hfem_.Node_r_ele()[e][0],
                        hfem_.Node_r_ele()[e][1]);
                    break;

                case 1:
                    vec_b_ele_[e][i] = gl_.qgauss([this, e](double r) { return r * hfem_.rho(r) * (r - hfem_.Node_r_ele()[e][0]) / hfem_.Length()[e]; },
                        hfem_.Node_r_ele()[e][0],
                        hfem_.Node_r_ele()[e][1]);

                    break;

                default:
                    BOOST_ASSERT(!"switch文のdefaultに来てしまった！");
                    break;
                }
            }
        }
    }

    void Hydrogen_Vh::make_global_matrix_and_vector()
    {
        // 全体行列と全体ベクトルを生成
        for (auto e = 0; e < hydrogen_fem::Hydrogen_FEM::ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    mat_A_glo_.coeffRef(hfem_.Node_num_seg()[e][i], hfem_.Node_num_seg()[e][j]) = mat_A_glo_.coeff(hfem_.Node_num_seg()[e][i], hfem_.Node_num_seg()[e][j]) + mat_A_ele_[e][i][j];
                }
                vec_b_glo_(hfem_.Node_num_seg()[e][i]) += vec_b_ele_[e][i];
            }
        }
    }

    // #endregion privateメンバ関数
}
