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

int main()
{
    hydrogen_vh::Hydrogen_Vh hl;
    hl.do_run();
    hl.save_result();

    return 0;
}
