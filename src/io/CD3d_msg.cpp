#include "io/CD3d_msg.h"

#define TURNOFF_ERR_MSG 0
#define TURNOFF_WRN_MSG 0
#define TURNOFF_OUT_MSG 0

//output error message
void cd_err(const cd_msg& msg)
{
#if !TURNOFF_ERR_MSG
    cerr<<"! ERROR: "<<msg.m_info<<endl;
#endif
}

//output warning message
void cd_wrn(const cd_msg& msg)
{
#if !TURNOFF_WRN_MSG
    cerr<<"* WARNING: "<<msg.m_info<<endl;
#endif
}

//output normal message
void cd_out(const cd_msg& msg)
{
#if !TURNOFF_OUT_MSG
    cerr<<"- "<<msg.m_info<<endl;
#endif
}
