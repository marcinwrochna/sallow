#pragma once
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <sstream>


/** Printf to stderr with `indentation`-many tabs. */
__attribute__((format(printf, 2, 3)))  // Treat second argument as printf-like format string.
static void log([[maybe_unused]] int indentation, [[maybe_unused]] const char* format, ...) {
    if (indentation > -1)
        return;
    va_list argptr;
    va_start(argptr, format);
    fprintf(stderr, "%3d: ", indentation);
    for (int i = 0; i < indentation; ++i)
        fprintf(stderr, "\t");
    vfprintf(stderr, format, argptr);
    va_end(argptr);
}


/** Convert vector to string (using `<<` operator). */
template <typename T>
std::string toString(const std::vector<T>& v) {
    std::ostringstream stream;
    bool first = true;
    for (const T& element : v) {
        if (!first)
            stream << " ";
        first = false;
        stream << element;
    }
    return stream.str();
}


/** Return ceil(log_2(n+1)), or length of binary representation of n. */
static constexpr int bitlength(int n) {
    int result = 0;
    while (n) {
        n /= 2;
        ++result;
    }
    return result;
}


/** Return whether a range includes a value. */
template <typename It, typename V>
bool contains(const It& b, const It& e, V value) {
    for (It i = b; i != e; ++i)
        if (*i == value)
            return true;
    return false;
}

/** Return whether a vector includes a value. */
template <typename V>
bool contains(const std::vector<V> vector, V value) {
    return contains(vector.begin(), vector.end(), value);
}


/** The Predicate is only taken by reference, so FilterView becomes invalid as soon as the Predicate goes out of scope. */
template <typename BaseView, typename Predicate>
class FilterView
{
  public:
    class Sentinel {};
    class const_iterator {
     public:
        const_iterator() = default;
        constexpr explicit const_iterator(const FilterView& view, typename BaseView::const_iterator it) noexcept
            : _it(it), _end(view._base.end()), _predicate(&view._predicate) {}
        constexpr auto operator*() const noexcept { return *_it; }
        constexpr const_iterator& operator++() noexcept { do { ++_it; } while(_it != _end && !(*_predicate)(*_it)); return *this; }
        constexpr void operator++(int) noexcept { do { ++_it; } while(_it != _end && !(*_predicate)(*_it)); }
        constexpr bool operator!=(const Sentinel&) const noexcept { return _it != _end; }
     private:
        typename BaseView::const_iterator _it;
        typename BaseView::const_iterator _end;
        const Predicate* _predicate = nullptr;
    };

    constexpr explicit FilterView(BaseView base, const Predicate& predicate) noexcept
        : _base(base), _predicate(predicate)  {}
    constexpr const_iterator begin() const noexcept {
        auto it = _base.begin();
        while (it != _base.end() && !_predicate(*it))
            ++it;
        return const_iterator{*this, it};
    }
    constexpr Sentinel end() const noexcept { return Sentinel{}; }
 private:
    const BaseView _base;
    const Predicate& _predicate;
};
