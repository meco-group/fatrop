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
#ifndef FATROP_INCLUDED
#define FATROP_INCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "ocp/OCPKKT.hpp"
#include "ocp/OCPAdapter.hpp"
#include "ocp/OCPAbstact.hpp"
#include "auxiliary/FatropVector.hpp"
#include "solver/FatropAlg.hpp"
#include "solver/FatropData.hpp"
#include "ocp/OCPScalingMethod.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "solver/AlgStrategy.hpp"
#include "solver/FatropOptions.hpp"
#include "function_evaluation/CasadiCodegen.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/StageOCPApplication.hpp"
#include "ocp/FatropOCPBuilder.hpp"

// #include "SparseSolvers/InterfaceMUMPS.hpp"
#endif //FATROP_INCLUDED