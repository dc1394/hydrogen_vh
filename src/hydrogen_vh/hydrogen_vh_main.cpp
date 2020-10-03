/*! \file hydrogen_vh_main.cpp
    \brief 水素原子のHartreeポテンシャルを有限要素法で求める
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

#include "hydrogen_vh.h"
#include <chrono>           // for std::chrono
#include <iostream>         // for std::cout
#include <boost/format.hpp> // for boost::format

int main()
{
    using namespace std::chrono;
    
    auto const start = system_clock::now();
    hydrogen_vh::Hydrogen_Vh hv;
    hv.do_run();
    auto const end = system_clock::now();
    std::cout << boost::format("計算が終わりました．\n計算時間 = %.14f（秒）\n") % duration_cast< duration<double> >(end - start).count();

    hv.save_result();
    std::cout << "計算結果を" << hydrogen_vh::Hydrogen_Vh::RESULT_FILENAME << "に書き込みました．" << std::endl;

    return 0;
}
