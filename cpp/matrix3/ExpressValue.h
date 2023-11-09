#pragma once

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include <istream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cassert>
using namespace std;

#include "Matrix.h"
#include "libLua/lua.hpp"

class ExpressValue
{
public:
    ExpressValue(const string &expr, int n);
    ~ExpressValue();

public:
    void setStr(const string &expr);
    void setNum(int n);
    void setValue(const string &key, double value) const;
    void setValue(const vector<pair<string, double>> &vp) const;

    double eval() const;

    Matrix F2(const Matrix &beta, const Matrix &X);
    Matrix statjacobian(const Matrix &X, const Matrix &beta);

private:
    lua_State *m_LuaStater;
    string str;
};
