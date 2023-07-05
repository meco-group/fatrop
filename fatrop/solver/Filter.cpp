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
#include "solver/Filter.hpp"
using namespace fatrop;
using namespace std;
Filter::Filter(const fatrop_int size)
{
    filterdata_.reserve(size + 1);
}
void Filter::augment(const FilterData &filterdata)
{
    filterdata_.push_back(filterdata);
}
inline bool Filter::a_dominmates_b(const FilterData &fdin0, const FilterData &fdin1) const
{
    // worse barrier filter and constraint violation -> dominated
    if (fdin0.phi > fdin1.phi && fdin0.theta > fdin1.theta)
    {
        return true;
    }
    return false;
}
bool Filter::is_acceptable(const FilterData &fdin) const
{
    // run over filterdata_ elements
    for (vector<double>::size_type k = 0; k < filterdata_.size(); k++)
    {
        if (a_dominmates_b(fdin, filterdata_.at(k)))
        {
            return false;
        };
    }
    // fdin is not dominated by one of the filterdata_ elements -> acceptable to filter
    return true;
}
void Filter::reset()
{
    filterdata_.resize(0);
}