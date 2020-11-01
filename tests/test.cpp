#include <random>
#include <datasets_generation_utils.hpp>
#include <algorithm>
#include <cstdlib>

#include "gtest/gtest.h"
#include "array.hpp"
#include "la_vector.hpp"
#include "benchmark_utils.hpp"

TEST(Operations_Test, big) {
    auto const u_pow = 26;
    int n = 10000000;
    std::srand(42210);
    std::mt19937 gen(42210);
    std::uniform_int_distribution<uint32_t> distribution(0, (1ul << u_pow) - 1);
    auto data = generate_unique(distribution, gen, n, true);
    array<uint32_t> arr(data);
    wrapper_sdsl<uint32_t, sdsl::sd_vector<>> ef_sd(data);
    wrapper_sdsl<uint32_t, sdsl::rrr_vector<>> rrr(data);
    la_vector<uint32_t, 7> lav(data);
    la_vector<uint32_t> lav_opt(data);

    constexpr auto dens = 8;

    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_delta_rank, dens>, std::vector<uint32_t>> v_delta(data);
    wrapper_sdsl_enc_vector<sdsl::enc_vector_rank<sdsl::coder::elias_gamma_rank, dens>, std::vector<uint32_t>> v_gamma(data);

    for (int i = 1; i < n; ++i) {
        auto arr_result = arr.select(i);
        auto ef_sd_result = ef_sd.select(i);
        auto ef_rrr_result = rrr.select(i);
        auto enc_vec_gamma_result = v_gamma.select(i);
        auto enc_vec_delta_result = v_delta.select(i);
        auto lav_result = lav.select(i);
        auto lav_result_opt = lav_opt.select(i);

        ASSERT_EQ(arr_result, ef_sd_result) << i;
        ASSERT_EQ(arr_result, ef_rrr_result) << i;
        ASSERT_EQ(arr_result, enc_vec_gamma_result) << i;
        ASSERT_EQ(arr_result, enc_vec_delta_result) << i;
        ASSERT_EQ(lav_result, arr_result) << i;
        ASSERT_EQ(lav_result_opt, lav_result_opt) << i;

        auto arr_result_rank = arr.rank(i);
        auto ef_sd_result_rank = ef_sd.rank(i);
        auto ef_rrr_result_rank = rrr.rank(i);
        auto enc_vec_gamma_result_rank = v_gamma.rank(i);
        auto enc_vec_delta_result_rank = v_delta.rank(i);
        auto lav_result_rank = lav.rank(i);
        auto lav_result_opt_rank = lav_opt.rank(i);

        ASSERT_EQ(arr_result_rank, ef_sd_result_rank) << i;
        ASSERT_EQ(arr_result_rank, enc_vec_gamma_result_rank) << i;
        ASSERT_EQ(arr_result_rank, enc_vec_delta_result_rank) << i;
        ASSERT_EQ(arr_result_rank, ef_rrr_result_rank) << i;
        ASSERT_EQ(lav_result_rank, arr_result_rank) << i;
        ASSERT_EQ(lav_result_opt_rank, arr_result_rank) << i;

    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
