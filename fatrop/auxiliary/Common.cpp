/*
 * Fatrop - A fast trajectory optimization solver
 *  Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "fatrop/auxiliary/Common.hpp"
#include <limits>
#include <climits>
#include <cstdlib>
#include <algorithm>

namespace fatrop
{
    bool CompareLessEqual(double lhs, double rhs)
    {
        double mach_eps = std::numeric_limits<double>::epsilon();
        return (lhs - rhs <= 10. * mach_eps * std::max(std::abs(rhs), std::abs(lhs)));
    }
    bool CompareLessEqual(double lhs, double rhs, double ref)
    {
        double mach_eps = std::numeric_limits<double>::epsilon();
        return (lhs - rhs <= 10. * mach_eps * std::abs(ref));
    }
}