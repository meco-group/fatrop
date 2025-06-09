//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_options_hpp__
#define __fatrop_common_options_hpp__

#include "fatrop/context/context.hpp"
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace fatrop
{
    /**
     * @brief Variant type for different option values.
     */
    typedef std::variant<Index, Scalar, bool, std::string> OptionVariant;

    enum class OptionTypeId
    {
        Scalar = 0, ///< Represents a floating-point option.
        Index = 1,  ///< Represents an integer option.
        Bool = 2,   ///< Represents a boolean option.
        String = 3, ///< Represents a string option.
        Unknown = 4 ///< Represents an unknown option type.
    };

    template <typename OptionType> class OptionTypeTraits
    {
    public:
        OptionTypeId get_option_type() const { return OptionTypeId::Unknown; }
    };

    template <> class OptionTypeTraits<Index>
    {
    public:
        OptionTypeId get_option_type() const { return OptionTypeId::Index; }
    };

    template <> class OptionTypeTraits<Scalar>
    {
    public:
        OptionTypeId get_option_type() const { return OptionTypeId::Scalar; }
    };

    template <> class OptionTypeTraits<bool>
    {
    public:
        OptionTypeId get_option_type() const { return OptionTypeId::Bool; }
    };

    template <> class OptionTypeTraits<std::string>
    {
    public:
        OptionTypeId get_option_type() const { return OptionTypeId::String; }
    };

    /**
     * @brief Base class for option setters, implementing the visitor pattern.
     *
     * This class defines the interface for option setters and provides default
     * implementations that throw exceptions for invalid types.
     */
    class OptionSetterBase
    {
    public:
        virtual void operator()(const Index &value)
        {
            throw std::runtime_error("Invalid type for option, got type Index (int).");
        }
        virtual void operator()(const Scalar &value)
        {
            throw std::runtime_error("Invalid type for option, got type Scalar (double).");
        }
        virtual void operator()(const bool &value)
        {
            throw std::runtime_error("Invalid type for option, got type bool.");
        }
        virtual void operator()(const std::string &value)
        {
            throw std::runtime_error("Invalid type for option, got type string.");
        }
        virtual OptionTypeId get_option_type() const { return OptionTypeId::Unknown; }
        virtual ~OptionSetterBase() = default;
    };

    /**
     * @brief Template class for option setters, implementing the visitor pattern.
     *
     * @tparam OptionType The type of the option value.
     * @tparam AlgoType The type of the algorithm class.
     */
    template <typename OptionType, typename AlgoType> class OptionSetter : public OptionSetterBase
    {
    public:
        /**
         * @brief Construct a new Option Setter object.
         *
         * @param set_option Function pointer to the setter method in AlgoType.
         * @param algo Pointer to the algorithm object.
         */
        OptionSetter(void (AlgoType::*set_option)(const OptionType &), AlgoType *algo)
            : set_option(set_option), algo(algo)
        {
        }

        /**
         * @brief Set the option value.
         *
         * @param value The value to set.
         */
        void operator()(const OptionType &value) override { (algo->*set_option)(value); }

        inline OptionTypeId get_option_type() const
        {
            return OptionTypeTraits<OptionType>().get_option_type();
        }

    private:
        void (AlgoType::*set_option)(const OptionType &);
        AlgoType *algo;
    };

    /**
     * @brief Specialization of OptionSetter for bool, implementing the visitor pattern.
     *
     * This specialization allows for setting bool options and converting int to bool.
     *
     * @tparam AlgoType The type of the algorithm class.
     */
    template <typename AlgoType> class OptionSetter<bool, AlgoType> : public OptionSetterBase
    {
    public:
        /**
         * @brief Construct a new Option Setter object for bool options.
         *
         * @param set_option Function pointer to the setter method in AlgoType.
         * @param algo Pointer to the algorithm object.
         */
        OptionSetter(void (AlgoType::*set_option)(const bool &), AlgoType *algo)
            : set_option(set_option), algo(algo)
        {
        }

        /**
         * @brief Set the bool option value.
         *
         * @param value The bool value to set.
         */
        void operator()(const bool &value) override { (algo->*set_option)(value); }

        /**
         * @brief Set the bool option value from an int.
         *
         * @param value The int value to convert to bool and set.
         */
        void operator()(const Index &value) override { (algo->*set_option)(value != 0); }

        /**
         * @brief Set the bool option value from a string (yes/no).
         */

        void operator()(const std::string &value) override
        {
            if (value == "yes")
            {
                (algo->*set_option)(true);
            }
            else if (value == "no")
            {
                (algo->*set_option)(false);
            }
            else
            {
                throw std::runtime_error("Trying to set bool option with string, value should be "
                                         "\"yes\" or \"no\", got " +
                                         value + ".");
            }
        }

        inline OptionTypeId get_option_type() const
        {
            return OptionTypeTraits<bool>().get_option_type();
        }

    private:
        void (AlgoType::*set_option)(const bool &);
        AlgoType *algo;
    };

    /**
     * @brief Registry for managing algorithm options.
     *
     * This class allows for registering and setting options for different algorithm types.
     * Multiple options with the same name can be registered, and all registered setters
     * will be called when the option is set.
     */
    class OptionRegistry
    {
    public:
        friend inline std::ostream &operator<<(std::ostream &os, const OptionRegistry &registry);
        /**
         * @brief Register an option for an algorithm.
         *
         * @tparam OptionType The type of the option value.
         * @tparam AlgoType The type of the algorithm class.
         * @param option_name The name of the option.
         * @param set_option Function pointer to the setter method in AlgoType.
         * @param algo Pointer to the algorithm object.
         */
        template <typename OptionType, typename AlgoType>
        void register_option(const std::string &option_name,
                             void (AlgoType::*set_option)(const OptionType &), AlgoType *algo)
        {
            options[option_name].push_back(
                std::make_unique<OptionSetter<OptionType, AlgoType>>(set_option, algo));
        }

        template <typename AlgoType> void register_options(AlgoType &algo)
        {
            algo.register_options(*this);
        }

        OptionTypeId get_option_type(const std::string &option_name) const
        {
            if (options.find(option_name) == options.end())
                return OptionTypeId::Unknown;
            return options.at(option_name)[0]->get_option_type();
        }

        /**
         * @brief Set the value of a registered option.
         *
         * This method calls all registered setters for the given option name.
         *
         * @tparam ValueType The type of the value to set.
         * @param option_name The name of the option to set.
         * @param value The value to set.
         * @throws std::runtime_error if the option is not found or if there's an error setting the
         * value.
         */
        template <typename ValueType>
        void set_option(const std::string &option_name, const ValueType &value)
        {
            auto it = options.find(option_name);
            if (it == options.end())
            {
                throw std::runtime_error("Option " + option_name + " not found.");
            }
            try
            {
                for (const auto &setter : it->second)
                {
                    std::visit([&](auto &&arg) { (*setter)(arg); }, OptionVariant(value));
                }
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error("Error setting option '" + option_name + "': " + e.what());
            }
        }

    private:
        std::unordered_map<std::string, std::vector<std::unique_ptr<OptionSetterBase>>> options;
    };

    inline std::ostream &operator<<(std::ostream &os, const OptionRegistry &registry)
    {
        os << "Option registry with registered options :\n";
        for (const auto &[option_name, setters] : registry.options)
        {
            os << "  " << option_name << " of type " << (int) registry.get_option_type(option_name) <<"\n";
        }
        return os;
    }

} // namespace fatrop

#endif //__fatrop_common_options_hpp__
