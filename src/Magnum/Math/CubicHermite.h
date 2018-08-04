#ifndef Magnum_Math_CubicHermiteSpline_h
#define Magnum_Math_CubicHermiteSpline_h
/*
    This file is part of Magnum.

    Copyright © 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
              Vladimír Vondruš <mosra@centrum.cz>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
*/

/** @file
 * @brief Class @ref Magnum::Math::CubicHermitePoint, alias @ref Magnum::Math::CubicHermitePoint2D, @ref Magnum::Math::CubicHermitePoint3D, function @ref Magnum::Math::select(), @ref Magnum::Math::lerp(), @ref Magnum::Math::splerp()
 */

#include "Magnum/Math/Bezier.h"
#include "Magnum/Math/Functions.h"

namespace Magnum { namespace Math {

/**
@brief Cubic Hermite spline point

Represents a point on a [cubic Hermite spline](https://en.wikipedia.org/wiki/Cubic_Hermite_spline).

Unlike @ref Bezier, which describes a curve segment, this structure describes
a spline point @f$ \boldsymbol{p} @f$, with in-tangent @f$ \boldsymbol{m} @f$
and out-tangent @f$ \boldsymbol{n} @f$. This form is more suitable for
animation keyframe representation. The structure assumes the in/out tangents
to be in their final form, i.e. already normalized by length of their adjacent
segments.

Cubic Hermite splines are fully interchangeable with Bézier curves, see
@ref fromBezier() and @ref Bezier::fromCubicHermite() for more information
about the conversion.

@see @ref CubicHermitePoint2D, @ref CubicHermitePoint3D,
    @ref Magnum::CubicHermitePoint2D, @ref Magnum::CubicHermitePoint2Dd,
    @ref Magnum::CubicHermitePoint3D, @ref Magnum::CubicHermitePoint3Dd,
    @ref CubicBezier
*/
template<UnsignedInt dimensions, class T> class CubicHermitePoint {
    public:
        typedef T Type;             /**< @brief Underlying data type */

        enum: UnsignedInt {
            Dimensions = dimensions /**< Dimensions of control points */
        };

        /**
         * @brief Create cubic Hermite spline point from adjacent Bézier curve segments
         *
         * Given two adjacent cubic Bézier curve segments defined by points
         * @f$ \boldsymbol{a}_i @f$ and @f$ \boldsymbol{b}_i @f$,
         * @f$ i \in \{ 0, 1, 2, 3 \} @f$, the corresponding cubic Hermite
         * spline point @f$ \boldsymbol{p} @f$, in-tangent @f$ \boldsymbol{m} @f$
         * and out-tangent @f$ \boldsymbol{n} @f$ is defined as: @f[
         *      \begin{array}{rcl}
         *          \boldsymbol{m} & = & 3 (\boldsymbol{a}_3 - \boldsymbol{a}_2)
         *                           =   3 (\boldsymbol{b}_0 - \boldsymbol{a}_2) \\
         *          \boldsymbol{p} & = & \boldsymbol{a}_3 = \boldsymbol{b}_0 \\
         *          \boldsymbol{n} & = & 3 (\boldsymbol{b}_1 - \boldsymbol{a}_3)
         *                           =   3 (\boldsymbol{b}_1 - \boldsymbol{b}_0)
         *      \end{array}
         * @f]
         *
         * Expects that the two segments are adjacent (i.e., the endpoint of
         * first segment is the start point of the second). If you need to
         * create a cubic Hermite spline point that's at the beginning or at
         * the end of a curve, simply pass a dummy Bézier segment that
         * satisfies this constraint as the other parameter:
         *
         * @snippet MagnumMath.cpp CubicHermitePoint-fromBezier
         */
        constexpr static CubicHermitePoint<dimensions, T> fromBezier(const CubicBezier<dimensions, T>& a, const CubicBezier<dimensions, T>& b) {
            return CORRADE_CONSTEXPR_ASSERT(a[3] == b[0],
                "Math::CubicHermitePoint::fromBezier(): segments are not adjacent"),
                CubicHermitePoint<dimensions, T>{3*(a[3] - a[2]), a[3], 3*(b[1] - a[3])};
        }

        /**
         * @brief Default constructor
         *
         * Construct cubic Hermite spline point with all control points being
         * zero vectors.
         */
        constexpr /*implicit*/ CubicHermitePoint(ZeroInitT = ZeroInit) noexcept: _inTangent{ZeroInit}, _point{ZeroInit}, _outTangent{ZeroInit} {}

        /** @brief Construct cubic Hermite spline point without initializing its contents */
        explicit CubicHermitePoint(NoInitT) noexcept: _inTangent{NoInit}, _point{NoInit}, _outTangent{NoInit} {}

        /**
         * @brief Construct cubic Hermite spline point with given control points
         * @param inTangent     In-tangent @f$ \boldsymbol{m} @f$
         * @param point         Point @f$ \boldsymbol{p} @f$
         * @param outTangent    Out-tangent @f$ \boldsymbol{n} @f$
         */
        constexpr /*implicit*/ CubicHermitePoint(const Vector<dimensions, T>& inTangent, const Vector<dimensions, T>& point, const Vector<dimensions, T>& outTangent) noexcept: _inTangent{inTangent}, _point{point}, _outTangent{outTangent} {}

        /**
         * @brief Construct subic Hermite spline point from another of different type
         *
         * Performs only default casting on the values, no rounding or
         * anything else.
         */
        template<class U> constexpr explicit CubicHermitePoint(const CubicHermitePoint<dimensions, U>& other) noexcept: _inTangent{other._inTangent}, _point{other._point}, _outTangent{other._outTangent} {}

        /** @brief Equality comparison */
        bool operator==(const CubicHermitePoint<dimensions, T>& other) const {
            return _inTangent == other._inTangent && _point == other._point && _outTangent == other._outTangent;
        }

        /** @brief Non-equality comparison */
        bool operator!=(const CubicHermitePoint<dimensions, T>& other) const {
            return !operator==(other);
        }

        /** @brief In-tangent @f$ \boldsymbol{m} @f$ */
        Vector<dimensions, T>& inTangent() { return _inTangent; }
        /* returns const& so [] operations are also constexpr */
        constexpr const Vector<dimensions, T>& inTangent() const { return _inTangent; } /**< @overload */

        /** @brief Point @f$ \boldsymbol{p} @f$ */
        Vector<dimensions, T>& point() { return _point; }
        /* returns const& so [] operations are also constexpr */
        constexpr const Vector<dimensions, T>& point() const { return _point; } /**< @overload */

        /** @brief Out-tangent @f$ \boldsymbol{n} @f$ */
        Vector<dimensions, T>& outTangent() { return _outTangent; }
        /* returns const& so [] operations are also constexpr */
        constexpr const Vector<dimensions, T>& outTangent() const { return _outTangent; } /**< @overload */

    private:
        template<UnsignedInt, class> friend class CubicHermitePoint;

        Vector<dimensions, T> _inTangent;
        Vector<dimensions, T> _point;
        Vector<dimensions, T> _outTangent;
};

/**
@brief Two-dimensional cubic Hermite spline point

Convenience alternative to @cpp CubicHermitePoint<2, T> @ce. See
@ref CubicHermitePoint for more information.
@see @ref CubicHermitePoint3D, @ref Magnum::CubicHermitePoint2D,
    @ref Magnum::CubicHermitePoint2Dd
*/
#ifndef CORRADE_MSVC2015_COMPATIBILITY /* Multiple definitions still broken */
template<class T> using CubicHermitePoint2D = CubicHermitePoint<2, T>;
#endif

/**
@brief Three-dimensional cubic Hermite spline point

Convenience alternative to @cpp CubicHermitePoint<3, T> @ce. See
@ref CubicHermitePoint for more information.
@see @ref CubicHermitePoint2D, @ref Magnum::CubicHermitePoint3D,
    @ref Magnum::CubicHermitePoint3Dd
*/
#ifndef CORRADE_MSVC2015_COMPATIBILITY /* Multiple definitions still broken */
template<class T> using CubicHermitePoint3D = CubicHermitePoint<3, T>;
#endif

/**
@brief Cubic Hermite spline quaternion

A special case of @ref CubicHermitePoint, used for interpolating
@ref Math::Quaternion using splines

@see @ref Magnum::CubicHermiteQuaternion, @ref Magnum::CubicHermiteQuaterniond
*/
template<class T> class CubicHermiteQuaternion {
    public:
        typedef T Type;             /**< @brief Underlying data type */

        /**
         * @brief Default constructor
         *
         * Construct cubic Hermite spline point with all control points being
         * zero vectors.
         */
        constexpr /*implicit*/ CubicHermiteQuaternion(ZeroInitT = ZeroInit) noexcept: _inTangent{ZeroInit}, _point{ZeroInit}, _outTangent{ZeroInit} {}

        /** @brief Construct cubic Hermite spline point without initializing its contents */
        explicit CubicHermiteQuaternion(NoInitT) noexcept: _inTangent{NoInit}, _point{NoInit}, _outTangent{NoInit} {}

        /**
         * @brief Construct cubic Hermite spline point with given control points
         * @param inTangent     In-tangent @f$ \boldsymbol{m} @f$
         * @param point         Point @f$ \boldsymbol{p} @f$
         * @param outTangent    Out-tangent @f$ \boldsymbol{n} @f$
         */
        constexpr /*implicit*/ CubicHermiteQuaternion(const Quaternion<T>& inTangent, const Quaternion<T>& point, const Quaternion<T>& outTangent) noexcept: _inTangent{inTangent}, _point{point}, _outTangent{outTangent} {}

        /**
         * @brief Construct subic Hermite spline point from another of different type
         *
         * Performs only default casting on the values, no rounding or
         * anything else.
         */
        template<class U> constexpr explicit CubicHermiteQuaternion(const CubicHermiteQuaternion<U>& other) noexcept: _inTangent{other._inTangent}, _point{other._point}, _outTangent{other._outTangent} {}

        /** @brief Equality comparison */
        bool operator==(const CubicHermiteQuaternion<T>& other) const {
            return _inTangent == other._inTangent && _point == other._point && _outTangent == other._outTangent;
        }

        /** @brief Non-equality comparison */
        bool operator!=(const CubicHermiteQuaternion<T>& other) const {
            return !operator==(other);
        }

        /** @brief In-tangent @f$ \boldsymbol{m} @f$ */
        Quaternion<T>& inTangent() { return _inTangent; }
        /* returns const& so [] operations are also constexpr */
        constexpr const Quaternion<T>& inTangent() const { return _inTangent; } /**< @overload */

        /** @brief Point @f$ \boldsymbol{p} @f$ */
        Quaternion<T>& point() { return _point; }
        /* returns const& so [] operations are also constexpr */
        constexpr const Quaternion<T>& point() const { return _point; } /**< @overload */

        /** @brief Out-tangent @f$ \boldsymbol{n} @f$ */
        Quaternion<T>& outTangent() { return _outTangent; }
        /* returns const& so [] operations are also constexpr */
        constexpr const Quaternion<T>& outTangent() const { return _outTangent; } /**< @overload */

    private:
        template<UnsignedInt, class> friend class CubicHermiteQuaternion;

        Quaternion<T> _inTangent;
        Quaternion<T> _point;
        Quaternion<T> _outTangent;
};

/** @debugoperator{CubicHermitePoint} */
template<UnsignedInt dimensions, class T> Corrade::Utility::Debug& operator<<(Corrade::Utility::Debug& debug, const CubicHermitePoint<dimensions, T>& value) {
    debug << "CubicHermitePoint({" << Corrade::Utility::Debug::nospace;
    debug << Corrade::Utility::Debug::nospace << value.inTangent()[0] << Corrade::Utility::Debug::nospace;
    for(UnsignedInt i = 1; i != dimensions; ++i)
        debug << "," << value.inTangent()[i] << Corrade::Utility::Debug::nospace;
    debug << "}, {" << Corrade::Utility::Debug::nospace;
    debug << Corrade::Utility::Debug::nospace << value.point()[0] << Corrade::Utility::Debug::nospace;
    for(UnsignedInt i = 1; i != dimensions; ++i)
        debug << "," << value.point()[i] << Corrade::Utility::Debug::nospace;
    debug << "}, {" << Corrade::Utility::Debug::nospace;
    debug << Corrade::Utility::Debug::nospace << value.outTangent()[0] << Corrade::Utility::Debug::nospace;
    for(UnsignedInt i = 1; i != dimensions; ++i)
        debug << "," << value.outTangent()[i] << Corrade::Utility::Debug::nospace;
    return debug << "})";
}

/** @debugoperator{CubicHermiteQuaternion} */
template<UnsignedInt dimensions, class T> Corrade::Utility::Debug& operator<<(Corrade::Utility::Debug& debug, const CubicHermiteQuaternion<T>& value) {
    debug << "CubicHermiteQuaternion({" << Corrade::Utility::Debug::nospace;
    debug << Corrade::Utility::Debug::nospace << value.inTangent()[0] << Corrade::Utility::Debug::nospace;
    for(UnsignedInt i = 1; i != dimensions; ++i)
        debug << "," << value.inTangent()[i] << Corrade::Utility::Debug::nospace;
    debug << "}, {" << Corrade::Utility::Debug::nospace;
    debug << Corrade::Utility::Debug::nospace << value.point()[0] << Corrade::Utility::Debug::nospace;
    for(UnsignedInt i = 1; i != dimensions; ++i)
        debug << "," << value.point()[i] << Corrade::Utility::Debug::nospace;
    debug << "}, {" << Corrade::Utility::Debug::nospace;
    debug << Corrade::Utility::Debug::nospace << value.outTangent()[0] << Corrade::Utility::Debug::nospace;
    for(UnsignedInt i = 1; i != dimensions; ++i)
        debug << "," << value.outTangent()[i] << Corrade::Utility::Debug::nospace;
    return debug << "})";
}

/* Explicit instantiation for commonly used types */
#ifndef DOXYGEN_GENERATING_OUTPUT
extern template MAGNUM_EXPORT Corrade::Utility::Debug& operator<<(Corrade::Utility::Debug&, const CubicHermitePoint<2, Float>&);
extern template MAGNUM_EXPORT Corrade::Utility::Debug& operator<<(Corrade::Utility::Debug&, const CubicHermitePoint<3, Float>&);
extern template MAGNUM_EXPORT Corrade::Utility::Debug& operator<<(Corrade::Utility::Debug&, const CubicHermitePoint<2, Double>&);
extern template MAGNUM_EXPORT Corrade::Utility::Debug& operator<<(Corrade::Utility::Debug&, const CubicHermitePoint<3, Double>&);
#endif

/** @relatesalso CubicHermitePoint
@brief Constant interpolation of two cubic Hermite spline points
@param a     First value
@param b     Second value
@param t     Interpolation phase

Given segment points @f$ \boldsymbol{p}_i @f$, in-tangents @f$ \boldsymbol{m}_i @f$
and out-tangents @f$ \boldsymbol{n}_i @f$, the interpolated value @f$ \boldsymbol{p} @f$
at phase @f$ t @f$ is: @f[
    \boldsymbol{p}(t) = \begin{cases}
        \boldsymbol{p}_a, & t < 1 \\
        \boldsymbol{p}_b, & t \ge 1
    \end{cases}
@f]

Equivalent to calling @cpp Math::select(a.point(), b.point(), t) @ce.
@see @ref select(const T&, const T&, U),
    @ref lerp(const CubicHermitePoint<dimensions, T>&, const CubicHermitePoint<dimensions, T>&, T),
    @ref splerp(const CubicHermitePoint<dimensions, T>&, const CubicHermitePoint<dimensions, T>&, T)
*/
template<UnsignedInt dimensions, class T> Vector<std::size_t(dimensions), T> select(const CubicHermitePoint<dimensions, T>& a, const CubicHermitePoint<dimensions, T>& b, T t) {
    return select(a.point(), b.point(), t);
}

/** @relatesalso CubicHermiteQuaternion
@brief Constant interpolation of two cubic Hermite spline points
@param a     First value
@param b     Second value
@param t     Interpolation phase

*/
template<class T> Quaternion<T> select(const CubicHermiteQuaternion<T>& a, const CubicHermiteQuaternion<T>& b, T t) {
    return select(a.point(), b.point(), t);
}

/** @relatesalso CubicHermitePoint
@brief Linear interpolation of two cubic Hermite points
@param a     First value
@param b     Second value
@param t     Interpolation phase

Given segment points @f$ \boldsymbol{p}_i @f$, in-tangents @f$ \boldsymbol{m}_i @f$
and out-tangents @f$ \boldsymbol{n}_i @f$, the interpolated value @f$ \boldsymbol{p} @f$
at phase @f$ t @f$ is: @f[
    \boldsymbol{p}(t) = (1 - t) \boldsymbol{p}_a + t \boldsymbol{p}_b
@f]

Equivalent to calling @cpp Math::lerp(a.point(), b.point(), t) @ce.
@see @ref lerp(const T&, const T&, U),
    @ref select(const CubicHermitePoint<dimensions, T>&, const CubicHermitePoint<dimensions, T>&, T),
    @ref splerp(const CubicHermitePoint<dimensions, T>&, const CubicHermitePoint<dimensions, T>&, T)
*/
template<UnsignedInt dimensions, class T> Vector<std::size_t(dimensions), T> lerp(const CubicHermitePoint<dimensions, T>& a, const CubicHermitePoint<dimensions, T>& b, T t) {
    return lerp(a.point(), b.point(), t);
}

/** @relatesalso CubicHermitePoint
@brief Spline interpolation of two cubic Hermite points
@param a     First segment
@param b     Second segment
@param t     Interpolation phase

Given segment points @f$ \boldsymbol{p}_i @f$, in-tangents @f$ \boldsymbol{m}_i @f$
and out-tangents @f$ \boldsymbol{n}_i @f$, the interpolated value @f$ \boldsymbol{p} @f$
at phase @f$ t @f$ is: @f[
    \boldsymbol{p}(t) = (2 t^3 - 3 t^2 + 1) \boldsymbol{p}_a + (t^3 - 2 t^2 + t) \boldsymbol{n}_a + (-2 t^3 + 3 t^2) \boldsymbol{p}_b + (t^3 - t^2)\boldsymbol{m}_b
@f]

@see @ref select(const CubicHermitePoint<dimensions, T>&, const CubicHermitePoint<dimensions, T>&, T),
    @ref lerp(const CubicHermitePoint<dimensions, T>&, const CubicHermitePoint<dimensions, T>&, T)
*/
template<UnsignedInt dimensions, class T> Vector<std::size_t(dimensions), T> splerp(const CubicHermitePoint<dimensions, T>& a, const CubicHermitePoint<dimensions, T>& b, T t) {
    return (T(2)*pow<3>(t) - T(3)*pow<2>(t) + 1)*a.point() +
        (pow<3>(t) - T(2)*pow<2>(t) + t)*a.outTangent() +
        (T(-2)*pow<3>(t) + T(3)*pow<2>(t))*b.point() +
        (pow<3>(t) - pow<2>(t))*b.inTangent();
}

}}

#endif
