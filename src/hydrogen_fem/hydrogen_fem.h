﻿/*! \file hydrogen_fem.h
    \brief FEMで水素原子に対するSchrödinger方程式を解くクラスの宣言

    Copyright © 2019 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.    
*/

#ifndef _HYDROGEN_FEM_H_
#define _HYDROGEN_FEM_H_

#pragma once

#include "../utility/property.h"
#include <cstdint>                  // for std::int32_t
#include <valarray>                 // for std::valarray
#include <Eigen/Core>               // for Eigen::MatrixXd
#include <boost/multi_array.hpp>    // for boost::multi_array

namespace hydrogen_fem {
    //! A class.
    /*!
        FEMで水素原子に対するSchrödinger方程式を解くクラス
    */
    class Hydrogen_FEM final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
        */
        Hydrogen_FEM();

        //! A destructor.
        /*!
            デストラクタ
        */
        ~Hydrogen_FEM() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function.
        /*!
            実際に計算を行い、水素原子の基底状態のエネルギー固有値を返す
            \return 水素原子の基底状態のエネルギー固有値
        */
        void do_run();

        //! A public member function (const).
        /*!
            一次要素でρ(r)を計算して返す
            \param r rの値
            \return ρ(r)の値
        */
        double rho(double r) const;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            境界条件を設定する
        */
        void boundary_conditions();

        //! A private member function (const).
        /*!
            左辺の要素行列を計算する
            \param e e番目の要素の添字
            \param le e番目の線分要素の長さ
            \param p 左辺の要素行列の行
            \param q 左辺の要素行列の列
            \return 左辺の要素行列の要素
        */
        [[nodiscard]]
        double get_A_matrix_element(std::int32_t e, double le, std::int32_t p, std::int32_t q) const;

        //! A private member function (const).
        /*!
            右辺の要素行列を計算する
            \param e e番目の要素の添字
            \param le e番目の線分要素の長さ
            \param p 右辺の要素行列の行
            \param q 右辺の要素行列の列
            \return 右辺の要素行列の要素
        */
        [[nodiscard]]
        double get_B_matrix_element(std::int32_t e, double le, std::int32_t p, std::int32_t q) const;
                
        //! A private member function.
        /*!
            全体行列を生成する
        */
        void make_global_matrix();

        //! A private member function.
        /*!
            要素行列を求める
        */
        void make_element_matrix();

        //! A private member function.
        /*!
            各種データを生成する
        */
        void make_data();

        //! A private member function.
        /*!
            固有ベクトル（波動関数）を規格化する
        */
        void normalize();

        // #endregion privateメンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            各要素の長さを返す
            \return 各要素の長さ
        */
        utility::Property<std::valarray<double> const &> const Length;
        
        //! A property.
        /*!
            各要素のGlobal節点番号を返す
            \return 各要素のGlobal節点番号
        */
        utility::Property<boost::multi_array<std::int32_t, 2> const &> const Node_num_seg;

        //! A property.
        /*!
            各要素のLocal節点のr座標を返す
            \return 各要素のLocal節点のr座標
        */
        utility::Property<boost::multi_array<double, 2> const &> const Node_r_ele;

        // #endregion プロパティ

        // #region メンバ変数

    public:
        //! A private member variable (constant expression).
        /*!
            節点数
        */
        static auto constexpr NODE_TOTAL = 5000;

        //! A private member variable (constant expression).
        /*!
            要素数
        */
        static auto constexpr ELE_TOTAL = NODE_TOTAL - 1;

        //! A private member variable (constant expression).
        /*!
            積分区間の上限
        */
        static auto constexpr R_MAX = 50.0;

        //! A private member variable (constant expression).
        /*!
            積分区間の上限
        */
        static auto constexpr R_MIN = 0.0;

    private:
        //! A private member variable.
        /*!
            左辺の全体行列
        */
        Eigen::MatrixXd hg_;

        //! A private member variable.
        /*!
            各要素の長さ
        */
        std::valarray<double> length_;

        //! A private member variable.
        /*!
            左辺要素係数行列
        */
        boost::multi_array<double, 3> mat_A_ele_;

        //! A private member variable.
        /*!
            右辺要素係数行列
        */
        boost::multi_array<double, 3> mat_B_ele_;

        //! A private member variable.
        /*!
            各要素のGlobal節点番号
        */
        boost::multi_array<std::int32_t, 2> node_num_seg_;

        //! A private member variable.
        /*!
            各要素のLocal節点のr座標
        */
        boost::multi_array<double, 2> node_r_ele_;

        //! A private member variable.
        /*!
            Global節点のr座標
        */
        std::valarray<double> node_r_glo_;

        //! A private member variable.
        /*!
            固有ベクトル
        */
        Eigen::VectorXd phi_;

        //! A private member variable.
        /*!
            右辺の全体行列
        */
        Eigen::MatrixXd ug_;

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A public copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        Hydrogen_FEM(Hydrogen_FEM const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Hydrogen_FEM & operator=(Hydrogen_FEM const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
        
}

#endif  // _HYDROGEN_FEM_H_
