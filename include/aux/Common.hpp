#ifndef COMMONINCLUDED
#define COMMONINCLUDED
#include <cassert>
#if DEBUG
#define DBGASSERT(assertion) assert(assertion);
#else
#define DBGASSERT(assertion)
#endif
namespace fatrop
{
 
} // namespace fatrop
#endif //  COMMONINCLUDED