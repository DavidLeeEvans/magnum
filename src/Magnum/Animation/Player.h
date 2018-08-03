#ifndef Magnum_Animation_Player_h
#define Magnum_Animation_Player_h
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
 * @brief Class @ref Magnum::Animation::Player, enum @ref Magnum::Animation::State
 */

#include <chrono>
#include <vector>

#include "Magnum/Animation/Track.h"
#include "Magnum/Math/Range.h"

namespace Magnum { namespace Animation {

/**
@brief Player state

@see @ref Player
@experimental
*/
enum class State: UnsignedByte {
    /**
     * The animation clip is currently playing. Setting the state to
     * @ref State::Playing does nothing.
     */
    Playing,

    /**
     * The animation clip is currently paused. Setting the state to
     * @ref State::Playing starts playing from where it left, setting the state
     * to @ref State::Stopped stops the animation, setting the state to
     * @ref State::Paused does nothing.
     */
    Paused,

    /**
     * The animation clip is currently stopped. Setting the state to
     * @ref State::Playing starts playing from the beginning, attempting to set
     * the state to @ref State::Paused will retain the @ref State::Stopped
     * state, setting the state to @ref State::Stopped does nothing.
     */
    Stopped
};

/** @debugoperatorenum{State} */
MAGNUM_EXPORT Debug& operator<<(Debug& debug, State value);

namespace Implementation {
    template<class, class> struct DefaultScaler;
}

/**
@brief Animation player
@tparam T   Time type
@tparam K   Key type

Provides a generic way for querying interpolated results from multiple
animation tracks of distinct types from a single place, together with managing
their running state.

@section Animation-Player-setup Setting up

The @ref Player class is used by adding tracks to it and specifying what should
be done with interpolation result values. The simplest option is specifying a
destination location when adding the track using @ref add() --- that'll mean
you get a fresh set of animated values at your disposal after every iteration:

@snippet MagnumAnimation.cpp Player-usage

The Player stores just @ref TrackView instances, for every @ref Track instance
you have to ensure that it stays alive for the whole lifetime of the player
instance.

In case you need to apply the animated values using a setter, it's possible
to fire a callback every iteration. Note that the @ref addWithCallback()
function has also a typeless version taking just @cpp void* @ce user pointer
instead of a reference to a concrete type. Below is an example of animating
@ref SceneGraph object transformations using the
@ref SceneGraph::TranslationRotationScalingTransformation3D transformation
implementation:

@snippet MagnumAnimation.cpp Player-usage-callback

The @ref addWithCallbackOnChange() variant will fire the callback only if the
interpolated value changes, which is useful for triggering other events. See
@ref Animation-Player-higher-order "below" for an example. Lastly, there is
@ref addRawCallback() that allows for greater control and further performance
optimizations. See its documentation for a usage example code snippet.

By default, the @ref duration() of an animation is calculated implicitly from
all added tracks. You can use @ref setDuration() to specify a custom duration
--- if it extends beyond the keyframe values, values of begin/end keyframes
will be extrapolated according to @ref Extrapolation specified for every track;
if it will be shorter, only a slice of the animation will be played. The
animation is implicitly played only once, use @ref setPlayCount() to set a
number of repeats or make it repeat indefinitely.

@section Animation-Player-playback Animation playback

The @ref Player class doesn't access any global timer functionality, but
instead requires you to call its APIs with explicit time values. That allows
for greater flexibility and control over animation playback, among other
things.

By default, the player is in a @ref State::Stopped state. Call @ref play() with
a time value denoting the moment at which the animation should start. After
that, the @ref advance() function is meant to be called every frame with
a current time value. As long as the animation is playing, the @ref advance()
function will update track result destination locations with interpolated
values and/or fire user-defined callbacks described above.

Once the animation playback is finished (exhausing the whole @ref duration() of
all @ref playCount() iterations), the @ref advance() will update the
destination locations and/or fire user-defined callbacks with values that
correspond to @ref duration() end time. This is guaranteed to be always the
case in order to correctly "park" the animations --- even if your app would
freeze for a while and @ref advance() would get called later, the result values
will never be calculated from a key value that's outside @ref duration().

Calling @ref stop() immediately transfers @ref state() to @ref State::Stopped
and the next @ref advance() iteration will give out interpolated values
corresponding to the begin time of @ref duration(), again to "park" the
animation back to its initial state. After that, no more updates are done until
the animation is started again. Compared to when the animation stops by itself,
this will park it at the beginning, not at the end.

Calling @ref pause() while the animation is running immediately transfers the
animation state to @ref State::Paused and the next @ref advance() iteration
will give out interpolated values corresponding to a time that was passed to
the @ref pause() function. After that, no more updates are done until the
animation is resumed again with @ref play() or stopped with @ref stop().

The callbacks are only ever fired from within the @ref advance() function,
never from @ref pause(), @ref stop() or any other API.

For managing global application you can use @ref Timeline, @ref std::chrono
APIs or any other type that supports basic arithmetic. The time doesn't have to
be monotonic or have constant speed, but note that non-continuous and backward
time jumps may have worse performance than going monotonically forward. See
@ref Animation-Player-time-type "below" for more information about using
different time types.

@snippet MagnumAnimation.cpp Player-usage-playback

@section Animation-Player-time-type Using custom time/key types

In long-running apps it's not desirable to use @ref Magnum::Float "Float" for
global application time, since its precision will deteriorate over time. Even
after one hour the precision loss might start to get noticeable. To overcome
this problem, it's possible to specify a type for time values that's different
from type used for animation track keys. In contrast, using
@ref Magnum::Float "Float" for animation track key values is usually good
enough, as the tracks are never too long for this to become a problem --- and
if the tracks *are* long, you can always use a different key type for them as
well. A good choice is @ref std::chrono::nanoseconds as a time type and keeping
track key values as @ref Magnum::Float "Float" seconds:

@snippet MagnumAnimation.cpp Player-usage-chrono

While there's a builtin support for the above, you are free to use any other
type combination --- for that you need to provide a *scaler* function that will
take care of converting a time difference to play iteration index and key value
inside given iteration. The types should be implicitly constructible, and have
basic arithmetic and comparison operators. In order to reduce header size, the
@ref Player implementation is in a separate @ref Player.hpp file that you need
to include to get all needed template function definitions. See also
@ref compilation-speedup-hpp for more information.

@snippet MagnumAnimation-custom.cpp Player-usage-custom

@section Animation-Player-higher-order Higher-order players, animating time

Sometimes you might want to control multiple players at the same time or
animate player state. That's doable by creating specialized tracks that control
given player via a state change callback. By adding more tracks you can control
multiple players from a central location.

@snippet MagnumAnimation.cpp Player-higher-order

Besides state, you can also animate @ref setDuration() and @ref setPlayCount(),
but be aware that setting those while the animation is playing might cause
unwanted jumps and abrupt stops. Time is also completely in your control and
you can employ another @ref Player instance to speed it up or slow it down for
a particular animation:

@snippet MagnumAnimation.cpp Player-higher-order-animated-time

@section Animation-Player-explicit-specializations Explicit template specializations

The following specializations are explicitly compiled into the @ref Animation
library. For other specializations (e.g. using an integer key type) you have to
use the @ref Player.hpp implementation file to avoid linker errors. See also
@ref compilation-speedup-hpp for more information.

-   @ref Player "Player<Float, Float>"
-   @ref Player "Player<std::chrono::nanoseconds, Float>"

@experimental
*/
template<class T, class K
    #ifdef DOXYGEN_GENERATING_OUTPUT
    = T
    #endif
> class Player {
    public:
        /** @brief Time type */
        typedef T TimeType;

        /** @brief Key type */
        typedef K KeyType;

        /**
         * @brief Scaler function type
         *
         * The function gets time from when the animation started and combined
         * duration of all tracks; returns play iteration index and key value
         * inside given iteration.
         */
        typedef std::pair<UnsignedInt, K>(*Scaler)(T, K);

        /** @brief Constructor */
        explicit Player();

        /**
         * @brief Construct with a custom scaler function
         * @param scaler    Scaler function
         */
        explicit Player(Scaler scaler);

        /** @brief Copying is not allowed */
        Player(const Player<T, K>&) = delete;

        /** @brief Move constructor */
        /* The noexcept specifier was removed because it causes *extreme*
           issues on various Clang versions (iOS, Android, Emscripten)
           currently used on Travis CI. Not on my stock 6.0 or NDK r17 tho. */
        Player(Player<T, K>&&);

        ~Player();

        /** @brief Copying is not allowed */
        Player<T, K>& operator=(const Player<T, K>&) = delete;

        /** @brief Move assignment */
        /* Deliberately not noexcept, see above */
        Player<T, K>& operator=(Player<T, K>&&);

        /** @brief Time-to-key scaler */
        Scaler scaler() const { return _scaler; }

        /**
         * @brief Duration
         *
         * If the duration was not set explicitly using @ref setDuration(),
         * returns value calculated implicitly from all added tracks. If no
         * tracks are added, returns default-constructed value.
         */
        Math::Range1D<K> duration() const { return _duration; }

        /**
         * @brief Set duration
         *
         * The duration is initially a default-constructed value, then
         * calculated implicitly from added tracks. Setting it explicitly will
         * overwrite the implicitly calculated value. Adding a track after the
         * duration was set explicitly will extend the duration to span all
         * track durations.
         *
         * Setting a duration that extends beyond the keyframe values will
         * cause values of begin/end keyframes to be extrapolated according to
         * @ref Extrapolation specified for given track. Setting a shorter
         * duration will cause only a slice of all tracks to be played.
         *
         * Modifying this value while @ref state() is @ref State::Playing may
         * cause the animation to jump or abruptly stop after next call to
         * @ref advance().
         * @see @ref TrackView::duration()
         */
        Player<T, K>& setDuration(const Math::Range1D<K>& duration) {
            _duration = duration;
            return *this;
        }

        /** @brief Play count */
        UnsignedInt playCount() const { return _playCount; }

        /**
         * @brief Set play count
         *
         * By default, play count is set to @cpp 1 @ce, meaning the animation
         * @ref duration() is played once. Value of @cpp 0 @ce means the
         * animation is repeated indefinitely.
         *
         * Modifying this value while @ref state() is @ref State::Playing may
         * cause the animation to jump or abruptly stop after next call to
         * @ref advance().
         */
        Player<T, K>& setPlayCount(UnsignedInt count) {
            _playCount = count;
            return *this;
        }

        /**
         * @brief Whether the player is empty
         *
         * @see @ref size(), @ref add(), @ref addWithCallback(),
         *      @ref addWithCallbackOnChange(), @ref addRawCallback()
         */
        bool isEmpty() const;

        /**
         * @brief Count of tracks managed by this player
         *
         * @see @ref isEmpty(), @ref add(), @ref addWithCallback(),
         *      @ref addWithCallbackOnChange(), @ref addRawCallback()
         */
        std::size_t size() const;

        /**
         * @brief Track at given position
         *
         * Due to the type-erased nature of the player implementation, it's not
         * possible to know the exact track type.
         */
        const TrackViewStorage<K>& track(std::size_t i) const;

        /**
         * @brief Add a track with a result destination
         *
         * The @p destination is updated with new value after each call to
         * @ref advance() as long as the animation is playing.
         */
        template<class V, class R> Player<T, K>& add(const TrackView<K, V, R>& track, R& destination);

        /** @overload
         *
         * Note that track ownership is *not* transferred to the @ref Player
         * and you have to ensure that it's kept in scope for the whole
         * lifetime of the @ref Player instance.
         */
        template<class V, class R> Player<T, K>& add(const Track<K, V, R>& track, R& destination) {
            return add(TrackView<K, V, R>{track}, destination);
        }

        /**
         * @brief Add a track with a result callback
         *
         * The @p callback is called with current key value, interpolated
         * result value and the @p userData pointer after each call to
         * @ref advance() as long as the animation is playing. The key value is
         * guaranteed to never be outside of the @ref duration() ranage, with
         * the interpolated result always corresponding to that key value.
         *
         * See the overload below for a more convenient type-safe way to pass
         * user data.
         * @see @ref addRawCallback()
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R> Player<T, K>& addWithCallback(const TrackView<K, V, R>& track, void(*callback)(const K&, const R&, void*), void* userData = nullptr);
        #else
        /* Otherwise the user would be forced to use the + operator to convert
           a lambda to a function pointer and (besides being weird and
           annoying) it's also not portable because it doesn't work on MSVC
           2015 and older versions of MSVC 2017. OTOH, putting this in the docs
           wouldn say nothing about how the callback signature should look. */
        template<class V, class R, class Callback> Player<T, K>& addWithCallback(const TrackView<K, V, R>& track, Callback callback, void* userData = nullptr);
        #endif

        /** @overload
         *
         * Note that the track ownership is *not* transferred to the
         * @ref Player and you have to ensure that it's kept in scope for the
         * whole lifetime of the @ref Player instance.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R> Player<T, K>& addWithCallback(const Track<K, V, R>& track, void(*callback)(const K&, const R&, void*), void* userData = nullptr);
        #else /* See above why */
        template<class V, class R, class Callback> Player<T, K>& addWithCallback(const Track<K, V, R>& track, Callback callback, void* userData = nullptr) {
            return addWithCallback(TrackView<K, V, R>{track}, callback, userData);
        }
        #endif

        /**
         * @brief Add a track with a result callback
         *
         * Equivalent to calling the above with a lambda wrapper that casts
         * @cpp void* @ce back to @cpp T* @ce and dereferences it in order to
         * pass it to @p callback. There is no additional overhead compared to
         * the overload taking the @cpp void* @ce pointer, however see
         * @ref addRawCallback() for optimization possibilities.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R, class U> Player<T, K>& addWithCallback(const TrackView<K, V, R>& track, void(*callback)(const K&, const R&, U&), U& userData);
        #else /* See above why */
        template<class V, class R, class U, class Callback> Player<T, K>& addWithCallback(const TrackView<K, V, R>& track, Callback callback, U& userData);
        #endif

        /** @overload
         *
         * Note that the track ownership is *not* transferred to the
         * @ref Player and you have to ensure that it's kept in scope for the
         * whole lifetime of the @ref Player instance.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R, class U> Player<T, K>& addWithCallback(const Track<K, V, R>& track, void(*callback)(const K&, const R&, U&), U& userData);
        #else /* See above why */
        template<class V, class R, class U, class Callback> Player<T, K>& addWithCallback(const Track<K, V, R>& track, Callback callback, U& userData) {
            return addWithCallback(TrackView<K, V, R>{track}, callback, userData);
        }
        #endif

        /**
         * @brief Add a track with a result callback that's called on change
         *
         * A combination of @ref add() and @ref addWithCallback() --- during
         * each call to @ref advance(), as long as the animation is playing,
         * the new value is compared to @p destination. If the new value is
         * different from the stored one, @p callback is called and
         * @p destination is updated. Note that in order to keep the memory
         * management inside the player class simple, the value can't be cached
         * inside and you are required to provide the @p destination location.
         *
         * See the overload below for a more convenient type-safe way to pass
         * user data.
         * @see @ref addRawCallback()
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R> Player<T, K>& addWithCallbackOnChange(const TrackView<K, V, R>& track, void(*callback)(const K&, const R&, void*), R& destination, void* userData = nullptr);
        #else /* See above why */
        template<class V, class R, class Callback> Player<T, K>& addWithCallbackOnChange(const TrackView<K, V, R>& track, Callback callback, R& destination, void* userData = nullptr);
        #endif

        /** @overload
         *
         * Note that the track ownership is *not* transferred to the
         * @ref Player and you have to ensure that it's kept in scope for the
         * whole lifetime of the @ref Player instance.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R> Player<T, K>& addWithCallbackOnChange(const Track<K, V, R>& track, void(*callback)(const K&, const R&, void*), R& destination, void* userData = nullptr);
        #else /* See above why */
        template<class V, class R, class Callback> Player<T, K>& addWithCallbackOnChange(const Track<K, V, R>& track, Callback callback, R& destination, void* userData = nullptr) {
            return addWithCallbackOnChange(TrackView<K, V, R>{track}, callback, destination, userData);
        }
        #endif

        /**
         * @brief Add a track with a result callback that's called on change
         *
         * Equivalent to calling the above with a lambda wrapper that casts
         * @cpp void* @ce back to @cpp T* @ce and dereferences it in order to
         * pass it to @p callback. There is no additional overhead compared to
         * the overload taking the @cpp void* @ce pointer, however see
         * @ref addRawCallback() for optimization possibilities.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R, class U> Player<T, K>& addWithCallbackOnChange(const TrackView<K, V, R>& track, void(*callback)(const K&, const R&, void*), R& destination, U& userData);
        #else /* See above why */
        template<class V, class R, class U, class Callback> Player<T, K>& addWithCallbackOnChange(const TrackView<K, V, R>& track, Callback callback, R& destination, U& userData);
        #endif

        /** @overload
         *
         * Note that the track ownership is *not* transferred to the
         * @ref Player and you have to ensure that it's kept in scope for the
         * whole lifetime of the @ref Player instance.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R, class U> Player<T, K>& addWithCallbackOnChange(const Track<K, V, R>& track, void(*callback)(const K&, const R&, void*), R& destination, U& userData);
        #else
        template<class V, class R, class U, class Callback> Player<T, K>& addWithCallbackOnChange(const Track<K, V, R>& track, Callback callback, R& destination, U& userData) {
            return addWithCallbackOnChange(TrackView<K, V, R>{track}, callback, destination, userData);
        }
        #endif

        /**
         * @brief Add a track with a raw callback
         *
         * This is a low-level function meant to be used if you want to avoid
         * the extra overhead of an additional callback in @ref addWithCallback()
         * or @ref addWithCallbackOnChange(), want more flexibility in the user
         * callback or want to control the track interpolation directly --- for
         * example taking advantage of @ref TrackView::atStrict() or passing an
         * inlineable interpolator function instead of using the saved
         * interpolator function pointer.
         *
         * The callback takes the raw @ref TrackViewStorage reference (which
         * you need to cast to a correct type), the interpolated key and hint
         * that's meant to be passed to @ref TrackView::at(), the destination
         * pointer (equivalent to the one passed to @ref add()), user callback
         * pointer (which again needs to be cast to a correct type) and user
         * data pointer. The following code snippet shows implementation of the
         * @ref addWithCallbackOnChange() API using this function, using a
         * custom callback to add a value to a vector if it changes:
         *
         * @snippet MagnumAnimation.cpp Player-addRawCallback
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R, class Callback> Player<T, K>& addRawCallback(const TrackView<K, V, R>& track, void(*callback)(const TrackViewStorage<K>&, K, std::size_t&, void*, void(*)(), void*), void* destination, void(*userCallback)(), void* userData);
        #else
        template<class V, class R, class Callback> Player<T, K>& addRawCallback(const TrackView<K, V, R>& track, Callback callback, void* destination, void(*userCallback)(), void* userData);
        #endif

        /** @overload
         *
         * Note that the track ownership is *not* transferred to the
         * @ref Player and you have to ensure that it's kept in scope for the
         * whole lifetime of the @ref Player instance.
         */
        #ifdef DOXYGEN_GENERATING_OUTPUT
        template<class V, class R> Player<T, K>& addRawCallback(const Track<K, V, R>& track, void(*callback)(const TrackViewStorage<K>&, K, std::size_t&, void*, void(*)(), void*), void* destination, void(*userCallback)(), void* userData);
        #else
        template<class V, class R, class Callback> Player<T, K>& addRawCallback(const Track<K, V, R>& track, Callback callback, void* destination, void(*userCallback)(), void* userData) {
            return addRawCallback(TrackView<K, V, R>{track}, callback, destination, userCallback, userData);
        }
        #endif

        /**
         * @brief State
         *
         * The player is @ref State::Stopped by default.
         * @see @ref play(), @ref pause(), @ref stop(), @ref setState()
         */
        State state() const { return _state; }

        /**
         * @brief Play
         *
         * Starts playing all tracks added to the player at given @p startTime.
         * If @ref state() is already @ref State::Playing, the animation is
         * restarted from the beginning at @p startTiime. If @ref state() is
         * @ref State::Paused, the animation continues from the time that was
         * passed to @ref pause().
         *
         * If @p startTime is in the future (that is, time passed to the next
         * @ref advance() iteration will be less than @p startTime),
         * @ref advance() will do nothing until given point in the future.
         * Setting time to such a particular value can be used to synchronize
         * playback of multiple independent animation clips.
         * @see @ref setState()
         */
        Player<T, K>& play(T startTime);

        /**
         * @brief Pause
         *
         * Pauses the currently playing animation at given @p pauseTime. If
         * @ref state() is not @ref State::Playing, the function does nothing.
         * See @ref advance() for a detailed description of behavior when the
         * animation gets paused.
         * @see @ref setState()
         */
        Player<T, K>& pause(T pauseTime);

        /**
         * @brief Stop
         *
         * Stops the currently playing animation. If @ref state() is
         * @ref State::Paused, discard the pause information. If @ref state()
         * is already @ref State::Stopped, the function does nothing. See
         * @ref advance() for a detailed description of behavior when the
         * animation gets stopped.
         * @see @ref setState()
         */
        Player<T, K>& stop();

        /**
         * @brief Set state
         *
         * Convenience function that calls @ref play(), @ref pause() or
         * @ref stop() based on @p state. See documentation of these functions
         * for detailed description. The @p time parameter is used only when
         * @p state is @ref State::Playing or @ref State::Paused, it's ignored
         * for @ref State::Stopped.
         */
        Player<T, K>& setState(State state, T time);

        /**
         * @brief Advance the animation
         *
         * As long as @ref state() is @ref State::Playing, goes through all
         * tracks added with @ref add(), @ref addWithCallback() or
         * @ref addWithCallbackOnChange() in order they were added and updates
         * the destination locations and/or fires the callbacks with
         * interpolation results.
         *
         * If @ref state() is @ref State::Paused or @ref State::Stopped, the
         * function does nothing. If @p time is less than time that was passed
         * to @ref play(), the function does nothing. If @p time is large
         * enough that @ref duration() times @ref playCount() got exhausted,
         * the function will update destination locations and/or fire
         * user-defined callback with key and result values corresponding to
         * the end time of @ref duration() in order to correctly "park" the
         * animation. The state then becomes @ref State::Stopped and no more
         * updates are done until the animation is started again.
         *
         * If @ref pause() was called right before a particular @ref advance()
         * iteration, the function will update destination locations and/or
         * fire user-defined callbacks with key and result values corresponding
         * to the time passed to the @ref pause() call before to correctly
         * "park" the animation. After that, no more updates are done until the
         * animation is started again.
         *
         * If @ref stop() was called right before a particular @ref advance()
         * iteration, the function will update destination locations and/or
         * fire user-defined callbacks with key and result values corresponding
         * to the begin time of @ref duration() to correctly "park" the
         * animation back to its initial state. After that, no more updates are
         * done until the animation is started again.
         */
        Player<T, K>& advance(T time);

    private:
        struct Track;

        Player<T, K>& addInternal(const TrackViewStorage<K>& track, void (*advancer)(const TrackViewStorage<K>&, K, std::size_t&, void*, void(*)(), void*), void* destination, void(*userCallback)(), void* userCallbackData);

        std::vector<Track> _tracks;
        Math::Range1D<K> _duration;
        UnsignedInt _playCount{1};
        State _state{State::Stopped};
        T _startTime{}, _pauseTime{};
        Scaler _scaler;
};

template<class T, class K> template<class V, class R> Player<T, K>& Player<T, K>::add(const TrackView<K, V, R>& track, R& destination) {
    return addInternal(track,
        [](const TrackViewStorage<K>& track, K key, std::size_t& hint, void* destination, void(*)(), void*) {
            *static_cast<R*>(destination) = static_cast<const TrackView<K, V, R>&>(track).at(key, hint);
        }, &destination, nullptr, nullptr);
}

#ifndef DOXYGEN_GENERATING_OUTPUT
template<class T, class K> template<class V, class R, class Callback> Player<T, K>& Player<T, K>::addWithCallback(const TrackView<K, V, R>& track, Callback callback, void* userData) {
    auto callbackPtr = static_cast<void(*)(const K&, const R&, void*)>(callback);
    return addInternal(track,
        [](const TrackViewStorage<K>& track, K key, std::size_t& hint, void*, void(*callback)(), void* userData) {
            reinterpret_cast<void(*)(const K&, const R&, void*)>(callback)(key, static_cast<const TrackView<K, V, R>&>(track).at(key, hint), userData);
        }, nullptr, reinterpret_cast<void(*)()>(callbackPtr), userData);
}

template<class T, class K> template<class V, class R, class U, class Callback> Player<T, K>& Player<T, K>::addWithCallback(const TrackView<K, V, R>& track, Callback callback, U& userData) {
    auto callbackPtr = static_cast<void(*)(const K&, const R&, U&)>(callback);
    return addInternal(track,
        [](const TrackViewStorage<K>& track, K key, std::size_t& hint, void*, void(*callback)(), void* userData) {
            reinterpret_cast<void(*)(const K&, const R&, U&)>(callback)(key, static_cast<const TrackView<K, V, R>&>(track).at(key, hint), *static_cast<U*>(userData));
        }, nullptr, reinterpret_cast<void(*)()>(callbackPtr), &userData);
}

template<class T, class K> template<class V, class R, class Callback> Player<T, K>& Player<T, K>::addWithCallbackOnChange(const TrackView<K, V, R>& track, Callback callback, R& destination, void* userData) {
    auto callbackPtr = static_cast<void(*)(const K&, const R&, void*)>(callback);
    return addInternal(track,
        [](const TrackViewStorage<K>& track, K key, std::size_t& hint, void* destination, void(*callback)(), void* userData) {
            R result = static_cast<const TrackView<K, V, R>&>(track).at(key, hint);
            if(result == *static_cast<R*>(destination)) return;
            reinterpret_cast<void(*)(const K&, const R&, void*)>(callback)(key, result, userData);
            *static_cast<R*>(destination) = result;
        }, &destination, reinterpret_cast<void(*)()>(callbackPtr), userData);
}

template<class T, class K> template<class V, class R, class U, class Callback> Player<T, K>& Player<T, K>::addWithCallbackOnChange(const TrackView<K, V, R>& track, Callback callback, R& destination, U& userData) {
    auto callbackPtr = static_cast<void(*)(const K&, const R&, U&)>(callback);
    return addInternal(track,
        [](const TrackViewStorage<K>& track, K key, std::size_t& hint, void* destination, void(*callback)(), void* userData) {
            R result = static_cast<const TrackView<K, V, R>&>(track).at(key, hint);
            if(result == *static_cast<R*>(destination)) return;
            reinterpret_cast<void(*)(const K&, const R&, U&)>(callback)(key, result, *static_cast<U*>(userData));
            *static_cast<R*>(destination) = result;
        }, &destination, reinterpret_cast<void(*)()>(callbackPtr), &userData);
}

template<class T, class K> template<class V, class R, class Callback> Player<T, K>& Player<T, K>::addRawCallback(const TrackView<K, V, R>& track, Callback callback, void* destination, void(*userCallback)(), void* userData) {
    auto callbackPtr = static_cast<void(*)(const TrackViewStorage<K>&, K, std::size_t&, void*, void(*)(), void*)>(callback);
    return addInternal(track, callbackPtr, destination, userCallback, userData);
}
#endif

#if defined(CORRADE_TARGET_WINDOWS) && !defined(__MINGW32__)
extern template class MAGNUM_EXPORT Player<Float, Float>;
extern template class MAGNUM_EXPORT Player<std::chrono::nanoseconds, Float>;
#endif

}}

#endif