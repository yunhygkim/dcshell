#ifndef _CD3D_MSG_H_
#define _CD3D_MSG_H_

#include <iostream>
#include <stdio.h>
#include <string>
#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
using namespace mathtool;
using namespace std;

//
//  This class is for global message output
//

class cd_msg{
public:
    cd_msg(){}
    cd_msg(const char * info){m_info=info;}
    cd_msg(const string& info){m_info=info;}
    cd_msg(double value){m_info=to_string(value);}
    cd_msg(int value){m_info=to_string(value);}
    cd_msg(unsigned int value){m_info=to_string(value);}
    // operator +
    cd_msg operator+(const cd_msg& other) const { return cd_msg(m_info+other.m_info); }
    cd_msg operator+(const string& s) const{ return cd_msg(m_info+s); }
    cd_msg operator+(const char * c) const{ return (*this)+string(c); }
    cd_msg operator+(double value) const { return (*this)+cd_msg(value); }
    cd_msg operator+(int value) const { return (*this)+cd_msg(value); }
    cd_msg operator+(unsigned int value) const {return (*this)+cd_msg(value);}
    cd_msg operator+(const Point3d& pt) const {return (*this)+to_string(pt);}
    cd_msg operator+(const Vector3d& vec) const {return (*this)+to_string(vec);}
    // operator +=
    void operator+=(cd_msg& other) {m_info=m_info+other.m_info;}
    void operator+=(const string& s){m_info=m_info+s;}
    void operator+=(const char * c){m_info=m_info+c;}
    void operator+=(double value){m_info=m_info+to_string(value);}
    void operator+=(int value){m_info=m_info+to_string(value);}
    void operator+=(unsigned int value){m_info=m_info+to_string(value);}
    void operator+=(const Point3d& pt) {m_info=m_info+to_string(pt);}
    void operator+=(const Vector3d& vec) {m_info=m_info+to_string(vec);}
    string m_info;
    //func
    string to_string(double value) const { 
        char info[32];
        sprintf(info,"%f",value);
        return info;
    }
    string to_string(int value) const { 
        char info[32];
        sprintf(info,"%d",value);
        return info;
    }
    string to_string(unsigned int value) const { 
        char info[32];
        sprintf(info,"%d",value);
        return info;
    }
    string to_string(const Point3d& pt) const {
        string info("(");
        return info+to_string(pt[0])+","+to_string(pt[1])+","+to_string(pt[2])+")";
    }
    string to_string(const Vector3d& pt) const {
        string info("(");
        return info+to_string(pt[0])+","+to_string(pt[1])+","+to_string(pt[2])+")";
    }
};


//output error message
void cd_err(const cd_msg& msg);

//output warning message
void cd_wrn(const cd_msg& msg);

//output normal message
void cd_out(const cd_msg& msg);

#endif //_CD3D_MSG_H_
