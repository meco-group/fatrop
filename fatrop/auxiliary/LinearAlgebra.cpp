/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
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
#include "auxiliary/LinearAlgebra.hpp"
namespace fatrop
{
    void FatropMat::print()
    {
        fatrop_int n_rows = nrows();
        fatrop_int n_cols = ncols();
        for (fatrop_int ai = 0; ai < n_rows; ai++)
        {
            for (fatrop_int aj = 0; aj < n_cols; aj++)
            {
                printf("%9.5f ", get_el(ai, aj));
            }
            printf("\n");
        }
    }
    void FatropVec::print()
    {
        fatrop_int n_el = nels();
        for (fatrop_int ai = 0; ai < n_el; ai++)
        {
            printf("%9.5f\n", get_el(ai));
        }
    }
}
