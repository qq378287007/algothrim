#include "ExpressValue.h"

static int l_log10(lua_State *L)
{
    double d = luaL_checknumber(L, 1);
    lua_pushnumber(L, log10(d));
    return 1;
}

ExpressValue::ExpressValue(const string &expr, int n)
{
    m_LuaStater = luaL_newstate();
    luaL_openlibs(m_LuaStater);

    luaL_dostring(m_LuaStater, "sin = math.sin; cos = math.cos;"
                               "tan = math.tan; abs = math.abs;"
                               "acos = math.acos; asin = math.asin;"
                               "atan = math.atan; exp = math.exp;"
                               "log = math.log; "
                               "sqrt = math.sqrt; int = math.floor;"
                               "min = math.min; max = math.max;");

    lua_pushcfunction(m_LuaStater, l_log10);
    lua_setglobal(m_LuaStater, "log10");

    lua_checkstack(m_LuaStater, 100);

    setStr(expr);
    setValue("x", 1.0);
    setNum(n);
}

ExpressValue::~ExpressValue()
{
    lua_close(m_LuaStater);
}

void ExpressValue::setStr(const string &expr)
{
    str = expr;
}

void ExpressValue::setNum(int n)
{
    for (int i = 0; i < n; i++)
        setValue("beta" + to_string(i), 1.0);
}

void ExpressValue::setValue(const string &key, double value) const
{
    const char *name = key.c_str();

    lua_getglobal(m_LuaStater, name);
    lua_pop(m_LuaStater, 1);
    lua_pushnumber(m_LuaStater, value);
    lua_setglobal(m_LuaStater, name);
}

void ExpressValue::setValue(const vector<pair<string, double>> &vp) const
{
    for (const pair<string, double> &p : vp)
        setValue(p.first, p.second);
}

double ExpressValue::eval() const
{
    luaL_dostring(m_LuaStater, ("return 0 + " + str).c_str());
    double value = lua_tonumber(m_LuaStater, -1);
    lua_pop(m_LuaStater, 1);
    return value;
}

Matrix ExpressValue::F2(const Matrix &beta, const Matrix &X)
{
    int m_Row = X.length();
    int m_Col = beta.length();

    Matrix Y(m_Row, 1);

    for (int i = 0; i < m_Col; i++)
        setValue("beta" + to_string(i), beta[i]);

    for (int i = 0; i < m_Row; i++)
    {
        setValue("x", X[i]);
        Y[i] = eval();
    }

    return Y;
}

Matrix ExpressValue::statjacobian(const Matrix &X, const Matrix &beta)
{
    const double eps = 2.220446049250313e-16;
    const double DerivStep = pow(eps, 1.0 / 3.0);
    const double nTheta = sqrt(beta.sum());
    const double DerivStep2 = DerivStep * (nTheta + (nTheta == 0.0));

    int m_Row = X.length();
    int m_Col = beta.length();

    Matrix J(m_Row, m_Col);

    for (int i = 0; i < m_Col; i++)
        setValue("beta" + to_string(i), beta[i]);

    Matrix thetaNew(1, m_Col);
    Matrix dy(m_Row, 1);
    Matrix y0(m_Row, 1);
    Matrix yplus(m_Row, 1);
    for (int i = 0; i < m_Row; i++)
    {
        setValue("x", X[i]);
        y0[i] = eval();
    }

    for (int i = 0; i < m_Col; i++)
    {
        double delta = DerivStep * beta[i];
        if (delta == 0.0)
            delta = DerivStep2;
        thetaNew = beta;
        thetaNew[i] += delta;

        for (int j = 0; j < m_Col; j++)
            setValue("beta" + to_string(j), thetaNew[j]);
            
        for (int j = 0; j < m_Row; j++)
        {
            setValue("x", X[j]);
            yplus[j] = eval();
        }

        dy = 1.0 / delta * (yplus - y0);

        for (int row = 0; row < m_Row; row++)
            J(row, i) = dy[row];
    }

    return J;
}