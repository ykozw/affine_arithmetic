//
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
//
#include <cstdio>
#include <map>
#include <vector>
#include <format>
#include <set>
#include <ranges>

// Interval arithmetic
class IA
{
public:
    //
    IA() = default;
    //
    explicit IA(const float low, const float high)
    {
        if (false)
        {
            low_ = std::nextafter(low_, -std::numeric_limits<float>::infinity());
            high_ = std::nextafter(high_, std::numeric_limits<float>::infinity());
        }
        else
        {
            low_ = low;
            high_ = high;
        }
    }
    //
    IA operator + (const IA& rhs) const
    {
        static_assert(std::numeric_limits<float>::is_iec559);
        return IA((low_ + rhs.low_), (high_ + rhs.high_));
    }
    //
    IA operator + (const float fv) const
    {
        return IA((low_ + fv), (high_ + fv));
    }
    //
    friend IA operator + (const float fv, const IA& ia)
    {
        return ia + fv;
    }
    //
    IA operator - (const IA& rhs) const
    {
        return IA((low_ - rhs.high_), (high_ - rhs.low_));
    }
    //
    friend IA operator-(float f, const IA& ia)
    {
        return IA(f, f) - ia;
    }
    //
    IA operator * (const IA& other) const
    {
        const float x0 = high_ * other.high_;
        const float x1 = high_ * other.low_;
        const float x2 = low_ * other.high_;
        const float x3 = low_ * other.low_;
        return
            IA(
                std::min({ x0, x1,x2,x3 }),
                std::max({ x0, x1,x2,x3 })
            );
    }
    //
    IA operator / (const IA& other) const
    {
        IA r;
        if ((other.low_ <= 0.0f) &&
            (0.0f <= other.high_))
        {
            return
                IA(
                    -std::numeric_limits<float>::infinity(),
                    +std::numeric_limits<float>::infinity());
        }
        else
        {
            const float x0 = high_ / other.high_;
            const float x1 = high_ / other.low_;
            const float x2 = low_ / other.high_;
            const float x3 = low_ / other.low_;
            return
                IA(
                    std::min({ x0, x1, x2, x3 }),
                    std::max({ x0, x1, x2, x3 }));
        }
        return r;
    }
    //
    IA sqrt() const
    {
        return
            IA(
                std::sqrtf(std::max(low_, 0.0f)),
                std::sqrtf(std::max(high_, 0.0f)));
    }
    //
    std::string toString() const
    {
        return std::format("[ {:.2f} : {:.2f} ]", low_, high_);
    }
    //
    float low() const
    {
        return low_;
    }
    //
    float high() const
    {
        return high_;
    }

private:
    float low_ = 0.0f;
    float high_ = 0.0f;
};

// Affine arithmetic
class AA
{
public:
    //
    AA() = default;
    //
    AA(float x)
    {
        cv_ = x;
    }
    //
    AA(const AA& aa)
    {
        cv_ = aa.cv_;
        deviations_ = aa.deviations_;
    }
    //
    AA(const IA& ia)
    {
        cv_ = (ia.low() + ia.high()) * 0.5f;
        deviations_[nextID()] = (ia.high() - ia.low()) * 0.5f;
    }
    //
    AA& operator = (const AA& other)
    {
        cv_ = other.cv_;
        deviations_ = other.deviations_;
        return *this;
    }
    //
    AA operator + (const AA& other) const
    {
        AA ret;
        ret.cv_ = cv_ + other.cv_;
        //
        for (auto& [k, v] : deviations_)
        {
            ret.deviations_[k] += v;
        }
        for (auto& [k, v] : other.deviations_)
        {
            ret.deviations_[k] += v;
        }
        return ret;
    }
    //
    AA operator + (const float fv) const
    {
        return (*this) + AA(fv);
    }
    //
    friend AA operator + (const float fv, AA& aa)
    {
        return aa + fv;
    }
    //
    AA operator - (const AA& other) const
    {
        AA ret;
        ret.cv_ = cv_ - other.cv_;
        //
        for (auto& [k, v] : deviations_)
        {
            ret.deviations_[k] += v;
        }
        for (auto& [k, v] : other.deviations_)
        {
            ret.deviations_[k] -= v;
        }
        return ret;
    }
    //
    friend AA operator -(float fv, const AA& aa)
    {
        return AA(fv) - aa;
    }
    // "Affine Arithmetic and its Applications to Computer Graphics." の実装。
    // aaflibはまた違った実装になっている
    AA operator * (const AA& other) const
    {
        const auto getDeviation = [](const auto& dev, const int id)
        {
            if (auto ite = dev.find(id); ite != dev.end())
            {
                return ite->second;
            }
            else
            {
                return 0.0f;
            }
        };
        AA ret;
        ret.cv_ = cv_ * other.cv_;
        //
        const float x0 = cv_;
        const float y0 = other.cv_;
        for (const auto& id : allIDs(*this, other))
        {
            const float xi = getDeviation(deviations_, id);
            const float yi = getDeviation(other.deviations_, id);
            ret.deviations_[id] = (x0 * yi + xi * y0);
        }
        //
        float U = 0.0f;
        for (auto [k, v] : deviations_)
        {
            U += std::abs(v);
        }
        float V = 0.0f;
        for (auto [k, v] : other.deviations_)
        {
            V += std::abs(v);
        }
        ret.deviations_[nextID()] = U * V;
        return ret;
    }
    //
    friend AA operator * (const AA& aa, const float fv)
    {
        AA ret;
        ret.cv_ = fv * aa.cv_;
        ret.deviations_ = aa.deviations_;
        for (auto& [k, v] : ret.deviations_)
        {
            v *= fv;
        }
        return ret;
    }
    //
    friend AA operator * (const float fv, const AA& aa)
    {
        return aa * fv;
    }
    //
    AA operator / (const AA& other) const
    {
        return
            (this == &other) ?
            AA(1.0) : (*this) * inv(other);
    }
    //
    friend AA operator / (const float fv, const AA& aa)
    {
        return inv(aa) * fv;
    }
    //
    std::string toString() const
    {
        auto ia = toIA();
        return std::format("[{:.2f},{:.2f}]", ia.low(), ia.high());
    }
    //
    float rad() const
    {
        float total = 0.0f;
        for (auto [k, v] : deviations_)
        {
            total += std::fabs(v);
        }
        return total;
    }
    //
    IA toIA() const
    {
        return IA{ cv_ - rad(), cv_ + rad() };
    }
    // sqrtのMinRange実装
    AA sqrt() const
    {
        if (deviations_.empty())
        {
            return AA(std::sqrt(cv_));
        }
        const float r = rad();
        // 区間
        float a = cv_ - r;
        float b = cv_ + r;
        float dneg = 0.0;
        // 区間での値
        float fa;
        if (a < 0.0f)
        {
            dneg = -a / 2.0f;
            a = 0.0f;
            fa = 0.0f;
        }
        else
        {
            fa = std::sqrt(a);
        }
        const float fb = std::sqrt(b);
        const float alpha = 1.0f / (2.0f * fb);
        const float delta = 0.5f * alpha * (a - 2.0f * fb * fa + b);
        const float dzeta = 0.5f * fb - delta;
        //
        AA ret;
        ret.cv_ = alpha * (cv_ + dneg) + dzeta;
        ret.deviations_ = deviations_;
        //
        const float alpha2 = alpha - dneg / rad();
        for (auto& [k, v] : ret.deviations_)
        {
            v *= alpha2;
        }
        ret.deviations_[AA::nextID()] = delta;
        return ret;
    }
    // inverseのMinRange実装
    static AA inv(const AA& aa)
    {
        if (aa.deviations_.empty())
        {
            return AA(1.0f / aa.cv_);
        }
        const float r = aa.rad();
        // 区間
        const float a = aa.cv_ - r;
        const float b = aa.cv_ + r;
        // 区間での値
        const float fa = 1.0f / a;
        const float fb = 1.0f / b;
        // 傾き
        float p;
        // Y軸切片
        float ya, yb;
        //
        if (a > 0.0f)
        {
            p = -fb / b;
            ya = fa - p * a;
            yb = 2.0f * fb;
        }
        else
        {
            p = -fa / a;
            ya = 2.0f * fa;
            yb = fb - p * b;
        }
        //
        AA ret;
        ret.cv_ = p * aa.cv_ + (0.5f * (ya + yb));
        ret.deviations_ = aa.deviations_;
        for (auto& [k, v] : ret.deviations_)
        {
            v *= p;
        }
        ret.deviations_[nextID()] = 0.5f * (ya - yb);
        return ret;
    }
    //
    static std::set<int> allIDs(const AA& aa0, const AA& aa1)
    {
        std::set<int> keys;
        auto k0 = std::views::keys(aa0.deviations_);
        auto k1 = std::views::keys(aa1.deviations_);
        keys.insert(k0.begin(), k0.end());
        keys.insert(k1.begin(), k1.end());
        return keys;
    }
    //
    static int nextID()
    {
        static int cnt = 0;
        return ++cnt;
    }

private:
    // Center Value
    float cv_ = 0.0f;
    // Key=変数ID、Value=その変数の偏微分
    std::map<int, float> deviations_;
};


void test()
{
    {
        const auto x = IA(4.0f, 6.0f);
        const auto y = x * (10.0f - x);
        assert(y.low() == 16.0f);
        assert(y.high() == 36.0f);
    }
    {
        const auto x = AA(IA(4.0f, 6.0f));
        const auto y = x * (10.0f - x);
        const auto iy = y.toIA();
        assert(iy.low() == 24.0f);
        assert(iy.high() == 26.0f);
    }
    {
        const auto x = AA(IA(4.0f, 6.0f));
        const auto y = x / x;
        const auto iy = y.toIA();
        assert(iy.low() == 1.0f);
        assert(iy.high() == 1.0f);
    }
    {
        auto x = AA(IA(4.0f, 6.0f));
        auto y = AA(IA(4.0f, 6.0f));
        auto z = (y * y) / (x * x);
        auto iz = z.toIA();
        assert(iz.low() == -0.0912699699f);
        assert(iz.high() == 2.57142878f);
    }
    // section 4.3 にある例題
    {
        auto x = IA(-2.0f, 2.0f);
        auto r = IA(-1.0f, 1.0f);
        auto s = IA(-1.0f, 1.0f);
        auto z = (10.0f + x + r) * (10.0 - x + s);
        assert(z.low() == 49.0f);
        assert(z.high() == 169.0f);
    }
    {
        auto x = AA(IA(-2.0f, 2.0f));
        auto r = AA(IA(-1.0f, 1.0f));
        auto s = AA(IA(-1.0f, 1.0f));
        auto z = (10.0f + x + r) * (10.0 - x + s);
        auto zi = z.toIA();
        assert(zi.low() == 71.0f);
        assert(zi.high() == 129.0f);
    }
}

//
PYBIND11_MODULE(AA, m)
{
    //
    test();

    // Interval
    namespace py = pybind11;
    py::class_<IA>(m, "IA")
        .def(py::init<float, float>())
        .def(py::self + py::self)
        .def(py::self + float())
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(float() - py::self)
        .def("low", &IA::low)
        .def("high", &IA::high)
        .def("sqrt", &IA::sqrt)
        .def("__repr__", &IA::toString)
        ;

    // Interval
    py::class_<AA>(m, "AA")
        .def(py::init<float>())
        .def(py::init<const AA&>())
        .def(py::init<const IA&>())
        .def(py::self + py::self)
        .def(py::self + float())
        .def(py::self - py::self)
        .def(float() - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(float() / py::self)
        .def("sqrt", &AA::sqrt)
        .def("toIA", &AA::toIA)
        .def("__repr__", &AA::toString)
        ;
}
