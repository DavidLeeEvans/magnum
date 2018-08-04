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

#include <Corrade/TestSuite/Tester.h>

#include "Magnum/Math/CubicHermite.h"
#include "Magnum/Math/Vector2.h"

namespace Magnum { namespace Math { namespace Test {

struct CubicHermiteTest: Corrade::TestSuite::Tester {
    explicit CubicHermiteTest();

    void construct();
    void constructDefault();
    void constructNoInit();
    void constructConversion();
    void constructFromBezier();
    void constructCopy();
    void convert();

    void data();

    void compare();

    void select();
    void lerp();

    void debug();
    void configuration();
};

CubicHermiteTest::CubicHermiteTest() {
    addTests({&CubicHermiteTest::construct,
              &CubicHermiteTest::constructDefault,
              &CubicHermiteTest::constructNoInit,
              &CubicHermiteTest::constructConversion,
              &CubicHermiteTest::constructFromBezier,
              &CubicHermiteTest::constructCopy,
              &CubicHermiteTest::convert,

              &CubicHermiteTest::data,

              &CubicHermiteTest::compare,

              &CubicHermiteTest::select,
              &CubicHermiteTest::lerp,

              &CubicHermiteTest::debug,
              &CubicHermiteTest::configuration});
}

typedef Math::Vector2<Float> Vector2;
typedef Math::CubicBezier2D<Float> CubicBezier2D;
typedef Math::CubicHermitePoint2D<Float> CubicHermitePoint2D;

void CubicHermiteTest::construct() {
}

void CubicHermiteTest::constructDefault() {
}

void CubicHermiteTest::constructNoInit() {
}

void CubicHermiteTest::constructConversion() {
}

void CubicHermiteTest::constructFromBezier() {
}

void CubicHermiteTest::constructCopy() {
}

void CubicHermiteTest::convert() {
}

void CubicHermiteTest::data() {
}

void CubicHermiteTest::compare() {
}

void CubicHermiteTest::select() {
}

void CubicHermiteTest::lerp() {
    CubicBezier2D bezier{Vector2{0.0f, 0.0f}, Vector2{10.0f, 15.0f}, Vector2{20.0f, 4.0f}, Vector2{5.0f, -20.0f}};
    /** @todo uhhhhh */
    auto a = CubicHermitePoint2D::fromBezier({Vector2{}, Vector2{}, Vector2{}, bezier[0]}, bezier);
    auto b = CubicHermitePoint2D::fromBezier(bezier, {bezier[3], Vector2{}, Vector2{}, Vector2{}});

    CORRADE_COMPARE(bezier.value(0.0f), (Vector2{0.0f, 0.0f}));
    CORRADE_COMPARE(Math::lerp(a, b, 0.0f), (Vector2{0.0f, 0.0f}));

    CORRADE_COMPARE(bezier.value(0.2f), (Vector2{5.8f, 5.984f}));
    CORRADE_COMPARE(Math::lerp(a, b, 0.2f), (Vector2{5.8f, 5.984f}));

    CORRADE_COMPARE(bezier.value(0.5f), (Vector2{11.875f, 4.625f}));
    CORRADE_COMPARE(Math::lerp(a, b, 0.5f), (Vector2{11.875f, 4.625f}));

    CORRADE_COMPARE(bezier.value(1.0f), (Vector2{5.0f, -20.0f}));
    CORRADE_COMPARE(Math::lerp(a, b, 1.0f), (Vector2{5.0f, -20.0f}));
}

void CubicHermiteTest::debug() {
}

void CubicHermiteTest::configuration() {
}

}}}

CORRADE_TEST_MAIN(Magnum::Math::Test::CubicHermiteTest)
