#ifndef FATROPOPTIONSINCLUDED
#define FATROPOPTIONSINCLUDED
#include <string>
namespace fatrop
{
    class FatropOption
    {
    public:
        FatropOption(const std::string &name, const std::string &descr) :name_(name), description_(descr){

                                            };
        std::string name_;
        std::string description_;
    };
    class NumericOption : public FatropOption
    {
    public:
        NumericOption(const std::string &name, const std::string &descr, double default_value) : FatropOption(name, descr), value_(default_value), default_value_(default_value) {}
        virtual void set_value(double value)
        {
            value_ = value;
        }
        double GetValue()
        {
            return value_;
        }

    private:
        double value_;
        double default_value_;
    };
    class NumericOptionBounded : public NumericOption
    {
        // decorator for NumericOption that also checks bounds
    public:
        NumericOptionBounded(const std::string &name, const std::string &descr, double default_value, double lower_bound, double upper_bound) : NumericOption(name, descr, default_value), lower_bound_(lower_bound), upper_bound_(upper_bound){};
        virtual void set_value(double value)
        {
            if ((value > lower_bound_) && (value < upper_bound_))
            {
                NumericOption::set_value(value);
            }
        }

    private:
        double lower_bound_;
        double upper_bound_;
    };
} // namespace fatrop
#endif