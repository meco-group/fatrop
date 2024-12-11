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
#ifndef __fatrop_solver_AlgBuilder_hpp__
#define __fatrop_solver_AlgBuilder_hpp__
// #include "NLPAlg.hpp"
#include "FatropAlg.hpp"
#include "fatrop/auxiliary/Common.hpp"
namespace fatrop
{
    class AlgBuilder
    {
    public:
        void build_fatrop_algorithm_objects(const std::shared_ptr<FatropNLP> &nlp, const std::shared_ptr<FatropNLP> &nlp_resto,
                                            const std::shared_ptr<FatropOptions> &fatropoptions,
                                            std::shared_ptr<FatropData> &fatropdata, std::shared_ptr<FatropData> &fatropdata_resto,
                                            std::shared_ptr<Journaller> &journaller)
        {
            if (fatropprinter_ == nullptr)
            {
                fatropprinter_ = std::make_shared<FatropPrinter>();
            }
            fatropdata = std::make_shared<FatropData>(nlp->get_nlp_dims(), fatropoptions, fatropprinter_);
            fatropdata_resto = std::make_shared<FatropData>(nlp_resto->get_nlp_dims(), fatropoptions, fatropprinter_);
            journaller = std::make_shared<Journaller>(fatropoptions->max_iter.get() + 1, fatropprinter_);
            fatropdata_ = fatropdata; // keep this around for building the algorithm
            fatropdata_resto_ = fatropdata_resto;
            journaller_ = journaller;
            nlp_ = nlp;
            nlp_resto_ = nlp_resto;
            fatropoptions_ = fatropoptions;
        }
        void set_printer(const std::shared_ptr<FatropPrinter> &printer)
        {
            fatropprinter_ = printer;
        }
        std::shared_ptr<FatropAlg> build_algorithm()
        {
            filter_ = std::make_shared<Filter>(fatropoptions_->max_iter.get() + 1);
            linesearch_ = std::make_shared<BackTrackingLineSearch>(fatropoptions_, nlp_, fatropdata_, filter_, journaller_, fatropprinter_);
            filter_resto_ = std::make_shared<Filter>(fatropoptions_->max_iter.get() + 1);
            linesearch_resto_ = std::make_shared<BackTrackingLineSearch>(fatropoptions_, nlp_resto_, fatropdata_resto_, filter_resto_, journaller_, fatropprinter_);
            orig_alg_ = std::make_shared<FatropAlg>(nlp_, fatropdata_, fatropoptions_, filter_, linesearch_, journaller_, fatropprinter_, nullptr, nullptr, false);
            resto_alg_ = std::make_shared<FatropAlg>(nlp_resto_, fatropdata_resto_, fatropoptions_, filter_resto_, linesearch_resto_, journaller_, fatropprinter_, orig_alg_, nullptr, true);
            orig_alg_->set_resto_alg(resto_alg_);
            return orig_alg_;
        }

        void update_options(const FatropOptions& options)
        {
            nlp_ -> update_options(options);
            nlp_resto_ -> update_options(options);
            fatropdata_ -> update_options(options);
            fatropdata_resto_ -> update_options(options);
            filter_ -> reserve(options.max_iter.get() + 1);
            linesearch_ -> update_options(options);
            filter_resto_ -> reserve(options.max_iter.get() + 1);
            linesearch_resto_ -> update_options(options);
            orig_alg_ -> update_options(options);
            resto_alg_ -> update_options(options);
            fatropprinter_ -> print_level() = options.print_level.get();
        };

    private:
        std::shared_ptr<FatropNLP> nlp_;
        std::shared_ptr<FatropNLP> nlp_resto_;
        std::shared_ptr<FatropOptions> fatropoptions_;
        std::shared_ptr<FatropData> fatropdata_;
        std::shared_ptr<FatropData> fatropdata_resto_;
        std::shared_ptr<Journaller> journaller_;
        std::shared_ptr<FatropPrinter> fatropprinter_;
        std::shared_ptr<Filter> filter_;
        std::shared_ptr<LineSearch> linesearch_;
        std::shared_ptr<Filter> filter_resto_;
        std::shared_ptr<LineSearch> linesearch_resto_;
        std::shared_ptr<FatropAlg> orig_alg_;
        std::shared_ptr<FatropAlg> resto_alg_;
    };
};

#endif // __fatrop_solver_AlgBuilder_hpp__