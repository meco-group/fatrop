#ifndef FATROPOPTIONSINCLUDED
#define FATROPOPTIONSINCLUDED
#include <string>
using namespace std;
namespace fatrop
{
    class FatropOption
    {
    public:
        FatropOption(const string &descr) : description_(descr){

                                            };
        string description_;
    };
    class NumericOption : public FatropOption
    {
    public:
        NumericOption(const string &descr, double default_value) : FatropOption(descr), value_(default_value), default_value_(default_value) {}
        virtual void SetValue(double value)
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
        NumericOptionBounded(const string &descr, double default_value, double lower_bound, double upper_bound) : NumericOption(descr, default_value), lower_bound_(lower_bound), upper_bound_(upper_bound){};
        virtual void SetValue(double value)
        {
            if ((value > lower_bound_) && (value < upper_bound_))
            {
                NumericOption::SetValue(value);
            }
        }

    private:
        double lower_bound_;
        double upper_bound_;
    };
} // namespace fatrop
#endif