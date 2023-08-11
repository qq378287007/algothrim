#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <iterator>
using namespace std;

vector<double> values{0.0, 0.1, 0.2, 0.3};
// vector<int> colors{0, 1, 2, 3, 4};

int value2color(double value)
{
    for (int i = 0; i < values.size(); i++)
        if (value <= values[i])
            return i;
    return values.size();
}

vector<int> getValue(double minValue, double maxValue)
{
    vector<int> data;
    for (int i = 0; i < values.size(); i++)
        if (values[i] > minValue && values[i] < maxValue)
            data.push_back(i);
    return data;
}

class Point
{
public:
    double x;
    double y;
    double value;
};

Point getPoint(const Point &p1, const Point &p2, double value)
{
    double x = (value - p1.value) / (p2.value - p1.value) * (p2.x - p1.x) + p1.x;
    double y = (value - p1.value) / (p2.value - p1.value) * (p2.y - p1.y) + p1.y;
    return Point{x, y, value};
}

vector<pair<vector<Point>, int>> getPair(const Point &minP, const Point &midP, const Point &maxP)
{
    double downValue = *values.begin();
    double upValue = *values.end();

    double minValue = minP.value;
    double midValue = midP.value;
    double maxValue = maxP.value;

    vector<pair<vector<Point>, int>> m_P;
    if (minValue >= upValue || maxValue <= downValue || (minValue == midValue && midValue == maxValue))
    {
        return vector<pair<vector<Point>, int>>{make_pair(vector<Point>{minP, midP, maxP}, value2color(maxValue))};
    }
    else if (midValue == minValue)
    {
        vector<int> data = getValue(minValue, maxValue);
        if (data.empty())
            return vector<pair<vector<Point>, int>>{make_pair(vector<Point>{minP, midP, maxP}, value2color(maxValue))};

        vector<pair<vector<Point>, int>> m_P;
        Point p1 = minP;
        Point p2 = midP;
        for (int i = 0; i < data.size(); i++)
        {
            double pValue = values[data[i]];
            Point p3 = getPoint(midP, maxP, pValue);
            Point p4 = getPoint(minP, maxP, pValue);
            vector<Point> ps{p1, p2, p3, p4};
            m_P.push_back(make_pair(ps, value2color(pValue)));
            p1 = p4;
            p2 = p3;
        }
        m_P.push_back(make_pair(vector<Point>{p1, p2, maxP}, value2color(maxValue)));
        return m_P;
    }
    else if (midValue == maxValue)
    {
        vector<int> data = getValue(minValue, maxValue);
        if (data.empty())
            return vector<pair<vector<Point>, int>>{make_pair(vector<Point>{minP, midP, maxP}, value2color(maxValue))};

        vector<pair<vector<Point>, int>> m_P;
        Point p1 = midP;
        Point p2 = maxP;
        for (int i = data.size() - 1; i >= 0; i--)
        {
            double pValue = values[data[i]];
            Point p3 = getPoint(minP, maxP, pValue);
            Point p4 = getPoint(minP, midP, pValue);
            vector<Point> ps{p1, p2, p3, p4};
            m_P.push_back(make_pair(ps, value2color(p1.value)));
            p1 = p4;
            p2 = p3;
        }
        m_P.push_back(make_pair(vector<Point>{p1, p2, minP}, value2color(p1.value)));
        return m_P;
    }
    else
    {
        Point tmpP = getPoint(minP, maxP, midValue);
        auto v1 = getPair(minP, midP, tmpP);
        auto v2 = getPair(tmpP, midP, maxP);

        vector<pair<vector<Point>, int>> v3;
        v3.resize(v1.size() + v2.size());
        for (auto iter = v1.cbegin(); iter != v1.cend(); iter++)
            v3.push_back(*iter);
        for (auto iter = v2.cbegin(); iter != v2.cend(); iter++)
            v3.push_back(*iter);
        return v3;
    }

    return m_P;
}

int main()
{
    return 0;
}
