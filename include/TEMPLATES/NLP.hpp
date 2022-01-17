#ifndef NLPINCLUDED
#define NLPINCLUDED
class NLP
{
    virtual int EvalHess() = 0;
    virtual int EvalJac() = 0;
};
#endif // NLPINCLUDED