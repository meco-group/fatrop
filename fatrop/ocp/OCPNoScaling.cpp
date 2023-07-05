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
#include "ocp/OCPNoScaling.hpp"
using namespace fatrop;
using namespace std;
OCPNoScaling::OCPNoScaling(const shared_ptr<FatropOptions> &fatrop_params) : OCPScalingMethod(fatrop_params){};


fatrop_int OCPNoScaling::compute_scalings(
    OCPKKTMemory *OCP,
    double &obj_scale,
    FatropVecBF &x_scales,
    FatropVecBF &lam_scales,
    const FatropVecBF &grad_curr) 
{
    obj_scale = 1.0;
    VECSE(x_scales.nels(), 1.0, (VEC *)x_scales, 0);
    VECSE(lam_scales.nels(), 1.0, (VEC *)lam_scales, 0);
    return 0;
};