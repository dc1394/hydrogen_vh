/*! \file hydrogen_vh.h
    \brief 水素原子のHartreeポテンシャルを求めるクラスの宣言
    Copyright © 2020 @dc1394 All Rights Reserved.

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

#ifndef _HYDROGEN_VH_H_
#define _HYDROGEN_VH_H_

#pragma once

#include "gausslegendre/gausslegendre.h"
#include "../hydrogen_fem/hydrogen_fem.h"
#include <cstdint>                      // for std::int32_t
#include <valarray>                     // for std::valarray
#include <boost/multi_array.hpp>        // for boost::multi_array
#include <Eigen/Core>                   // for Eigen::MatrixXd, Eigen::VectorXd
#include <Eigen/Sparse>                 // for Eigen::SparseMatrix

namespace hydrogen_vh {
    //! A class.
    /*!
        VWN-LDAを用い、Kohn-Sham法でヘリウム原子のエネルギーを計算するクラス
    */
    class Hydrogen_Vh final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
        */
        Hydrogen_Vh();

        //! A destructor.
        /*!
            デストラクタ
        */
        ~Hydrogen_Vh() = default;

        // #region publicメンバ関数

        //! A public member function.
        /*!
            実際に計算を行う
        */
        void do_run();

        //! A public member function.
        /*!
            計算結果をファイルに出力する
        */
        void save_result() const;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            境界条件処理を行う
        */
        void boundary_conditions();

        //! A private member function.
        /*!
            要素行列とLocal節点ベクトルの各成分を計算する
        */
        void make_element_matrix_and_vector();

        //! A private member function.
        /*!
            全体行列と全体ベクトルの各成分を計算する
        */
        void make_global_matrix_and_vector();

        // #endregion privateメンバ関数

        // #region メンバ変数

public:
        //! A public member variable (constant expression).
        /*!
            出力ファイル名
        */
        static auto constexpr RESULT_FILENAME = "result.csv";

private:
        //! A private member variable (constant expression).
        /*!
            Gauss-Legendre積分の分点
        */
        static auto constexpr INTEGTABLENUM = 100;

        //! A private member variable.
        /*!
            Gauss-Legendre積分用オブジェクト
        */
        gausslegendre::Gauss_Legendre gl_;

        //! A private member variable.
        /*!
            水素原子の有限要素法の解法オブジェクト
        */
        hydrogen_fem::Hydrogen_FEM hfem_;

        //! A private member variable.
        /*!
            左辺の係数行列
        */
        boost::multi_array<double, 3> mat_A_ele_;

        //! A private member variable.
        /*!
            左辺の全体行列
        */
        Eigen::SparseMatrix<double> mat_A_glo_;

        //! A private member variable.
        /*!
            連立方程式の解
        */
        Eigen::VectorXd u_;

        //! A private member variable.
        /*!
            右辺のLocal節点ベクトル
        */
        boost::multi_array<double, 2> vec_b_ele_;

        //! A private member variable.
        /*!
            右辺の全体ベクトル
        */
        Eigen::VectorXd vec_b_glo_;

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A public copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        Hydrogen_Vh(Hydrogen_Vh const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Hydrogen_Vh & operator=(Hydrogen_Vh const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
        
}

#endif  // _HYDROGEN_VH_H_
