#pragma once

#include <chrono>
#include <iostream>
#include "climits"

#include "array.hpp"
#include <sdsl/enc_vector.hpp>
#include <sdsl/dac_vector.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <index_types.hpp>
#include <succinct/mapper.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/io.hpp>
#include "datasets_generation_utils.hpp"
#include "wrapper_sdsl.hpp"

const int start_bpc = 6;
const int end_bpc = 16;
const int step_bpc = 1;

const int start_rrr = 4;
const int end_rrr = 7;
const int step_rrr = 1;

const int start_elias = 3;
const int end_elias = 8;
const int step_elias = 1;


constexpr int TIMES_TEST = 500000;
using timer = std::chrono::high_resolution_clock;

template<class T>
void do_not_optimize(T const &value) {
    asm volatile("" : : "r,m"(value) : "memory");
}


template<typename V, typename Q>
void test_select(const V &v, const Q &queries) {
    auto start = timer::now();
    auto cnt = 0;
    for (auto i = 0; i < queries.size(); ++i) {
        cnt += v.select(queries[i]);
    }
    do_not_optimize(cnt);

    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count();
    std::cout << "," << static_cast<double>(elapsed) / queries.size()
              << "," << v.bits_per_element() << std::flush;
}

template<typename V, typename Q>
void test_rank(const V &v, const Q &queries) {
    auto start = timer::now();
    auto cnt = 0;
    for (auto i = 0; i < queries.size(); ++i) {
        auto res_lower_bound = v.rank(queries[i]);
        if constexpr (std::is_scalar<decltype(res_lower_bound)>::value)
            cnt += res_lower_bound;
        else
            cnt += *res_lower_bound;
    }

    do_not_optimize(cnt);
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(timer::now() - start).count();
    std::cout << "," << static_cast<double>(elapsed) / queries.size() << std::flush;
}

template<typename V, typename D>
void test_linear_decompress(const V &v, const D &data) {
    std::vector<typename D::value_type> out_buffer(data.size());
    auto start = timer::now();
    v.decode(out_buffer.data());
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(timer::now() - start).count();

    if (out_buffer != data) {
        std::cerr << "Decompression failed" << std::endl;
        exit(1);
    }

    std::cout << "," << elapsed << std::flush;
}

template<int First, int Last, int step, typename Lambda>
inline void static_for(Lambda const &f) {
    if constexpr (First <= Last) {
        f(std::integral_constant<size_t, First>{});
        static_for<First + step, Last, step>(f);
    }
}

template<typename Sequence, class T, typename Qsize, typename Qvalue>
void test_ds2i(const std::vector<T> &dataset, const Qsize &rands1, const Qvalue &rands2) {
    ds2i::global_parameters params;
    succinct::bit_vector_builder docs_bits;
    T universe = dataset.back() + 1;
    Sequence::write(docs_bits, dataset.begin(), universe, dataset.size(), params);
    succinct::bit_vector bit_vector(&docs_bits);
    typename Sequence::enumerator enumerator(bit_vector, 0, universe, dataset.size(), params);
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto val = 0;
        for (size_t i = 0; i < TIMES_TEST; ++i) {
            val += enumerator.move(rands1[i]).second;
#ifndef NDEBUG
            if (dataset[rands1[i]] != val)
                exit(1);
#endif
        }
        do_not_optimize(val);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed =
                static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / TIMES_TEST;
        double bpk = static_cast<double>(((succinct::mapper::size_tree_of(bit_vector)->size) * 8)) / dataset.size();
        std::cout << "," << elapsed << "," << bpk << std::flush;
    }
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto val = 0;
        for (size_t i = 0; i < TIMES_TEST; ++i) {
            val += enumerator.next_geq(rands2[i]).second;
#ifndef NDEBUG
            if (dataset[rands2[i]] != val)
                exit(1);
#endif
        }
        do_not_optimize(val);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed =
                static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / TIMES_TEST;
        std::cout << "," << elapsed << std::flush;
    }
}


template<typename D, typename Qsize, typename Qvalue>
void test_one(D &d, Qsize &rands1, Qvalue &rands2) {
    test_select(d, rands1);
    test_rank(d, rands2);
}

template<bool la_vector_epsilon, bool la_vector_opt, bool elias_codes, bool ds2i, typename T>
void test_all(std::vector<T> &a) {
    std::mt19937 mt1(2323);
    std::uniform_int_distribution<size_t> dist1(1, a.size() - 1);
    std::vector<size_t> rands1(TIMES_TEST);
    for (auto i = 0; i < TIMES_TEST; ++i) {
        rands1[i] = (dist1(mt1));
    }
    std::mt19937 mt2(4242);
    std::uniform_int_distribution<T> dist2(a.front(), a.back() - 1);
    std::vector<T> rands2(TIMES_TEST);
    for (auto i = 0; i < TIMES_TEST; ++i) {
        rands2[i] = (dist2(mt2));
    }
    {
        array<T> arr(a);
        test_one(arr, rands1, rands2);
    }
    {
        wrapper_sdsl<T, sdsl::sd_vector<>> ef(a);
        test_one(ef, rands1, rands2);
    }
    static_for<start_rrr, end_rrr, step_rrr>([&a, &rands1, &rands2](auto e) {
        constexpr auto blocksize = (1u << e) - 1;
        wrapper_sdsl<T, sdsl::rrr_vector<blocksize>> rrr(a);
        test_one(rrr, rands1, rands2);
    });
    if constexpr (la_vector_epsilon) {
        static_for<start_bpc, end_bpc, step_bpc>([&](auto bits_per_correction) {
            la_vector<T, bits_per_correction> la_vec(a);
            test_one(la_vec, rands1, rands2);
            std::cout << "," << la_vec.segments_count();
        });
    }
    if constexpr (la_vector_opt) {
        la_vector<T> la_vec_opt(a);
        test_one(la_vec_opt, rands1, rands2);
        std::cout << "," << la_vec_opt.segments_count();
    }
    if constexpr (elias_codes) {
        static_for<start_elias, end_elias, step_elias>([&a, &rands1, &rands2](auto e) {
            constexpr auto dens = 1u << e;
            {
                wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_delta_rank, dens>, std::vector<T>> v(a);
                test_one(v, rands1, rands2);
            }
            {
                wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_gamma_rank, dens>, std::vector<T>> v(a);
                test_one(v, rands1, rands2);
            }
        });
    }
    if constexpr(ds2i) {
        //test_ds2i<ds2i::compact_elias_fano, T>(a, rands1, rands2);
        test_ds2i<ds2i::uniform_partitioned_sequence<>, T>(a, rands1, rands2);
        test_ds2i<ds2i::partitioned_sequence<>, T>(a, rands1, rands2);
    }
}

template<bool la_vector_epsilon, bool la_vector_opt, bool elias_codes, bool ds2i>
std::string get_csv_header() {
    std::string to_return = "i,n,u,ratio,mean_gap,stddev_gap";
    to_return.append(",array_time_select,array_bpk"
                     ",array_time_rank"

                     ",ef_sd_time_select,ef_sd_bpk"
                     ",ef_sd_time_rank"
    );
    static_for<start_rrr, end_rrr, step_rrr>([&to_return](auto e) {
        constexpr auto blocksize = (1u << e) - 1;
        std::string blocksize_str = std::to_string(blocksize);
        to_return.append(",ef_rrr_" + blocksize_str + "_time_select,ef_rrr_" + blocksize_str + "_bpk"
                         ",ef_rrr_" + blocksize_str + "_time_rank"
        );
    });
    if constexpr (la_vector_epsilon) {
        static_for<start_bpc, end_bpc, step_bpc>([&](uint8_t bits_per_correction) {
            std::string bpc = std::to_string(bits_per_correction);
            to_return.append(",la_vector_" + bpc + "_time_select,la_vector_" + bpc + "_bpk"
                                 ",la_vector_" + bpc + "_time_rank"
                                 ",la_vector_" + bpc + "_segments");
        });
    }
    if constexpr (la_vector_opt) {
        to_return.append(",la_vector_opt_time_select,la_vector_opt_bpk"
                         ",la_vector_opt_time_rank"
                         ",la_vector_opt_segments");
    }
    if constexpr (elias_codes) {
        static_for<start_elias, end_elias, step_elias>([&to_return](auto e) {
            auto e_str = std::to_string(e);
            to_return.append(",elias_delta_" + e_str + "_time_select,elias_delta_" + e_str + "_bpk,"
                             "elias_delta_" + e_str + "_time_rank");
            to_return.append(",elias_gamma_" + e_str + "_time_select,elias_gamma_" + e_str + "_bpk,"
                             "elias_gamma_" + e_str + "_time_rank");
        });
    }
    if constexpr(ds2i) {
        //to_return.append(",compact_select"
        //                 ",compact_bpk"
        //                 ",compact_rank");

        to_return.append(",uniform_select"
                         ",uniform_sequence_bpk"
                         ",uniform_rank");

        to_return.append(",partitioned_select"
                         ",partitioned_bpk"
                         ",partitioned_rank");
    }

    return to_return;
}
