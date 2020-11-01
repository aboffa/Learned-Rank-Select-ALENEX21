//
// Created by boffa on 20/02/20.
//

#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <sdsl/sd_vector.hpp>
#include <sdsl/rrr_vector.hpp>

template<typename T>
class array {
protected:
    std::vector<T> a;

public:
    array(std::vector<T> &a0) : a(a0) {}

    inline T operator[](size_t i) const {
        assert(i >= 0 && i < a.size());
        return a[i];
    }

    inline T select(size_t i) const {
        assert(i > 0 && i <= a.size());
        return a[i - 1];
    }

    inline double bits_per_element() const {
        return static_cast<double>(sizeof(T) * 8);
    };

    inline size_t rank(T x) const {
        const T *base = a.data();
        size_t n = a.size();
        while (n > 1) {
            size_t half = n / 2;
            base = (base[half] < x) ? base + half : base;
            n -= half;
        }
        return (*base < x) + base - a.data();
    }

    inline void decode(T *out) const {
        for (size_t i = 0; i < a.size(); ++i) {
            *(out++) = a[i];
        }
    }

    inline std::vector<T> decode() const {
        std::vector<T> out(a.size());
        decode(out.data());
        return out;
    }
};



#endif //ARRAY_HPP
