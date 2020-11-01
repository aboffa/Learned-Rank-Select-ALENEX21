#include <iostream>
#include <cstdint>
#include <vector>
#include <random>

#include "la_vector.hpp"
#include "benchmark_utils.hpp"
#include "index_reader.hpp"


int main(int argc, char *argv[]) {
    std::cout << std::scientific;
    if (argc < 2) {
        std::cout << "Usage: "+std::string(argv[0])+" <file_1> <file_2> ... <file_n> \n"
                     "Files must be generated using the function 'write_data_binary' in the file generate_datasets.cpp.\n"
                     "This program compares the select time, rank time and bit per integer of:\n"
                     " - plain 32-bit vector (std::vector)\n"
                     " - Elias-Fano (sdsl::sd_vector)\n"
                     " - rrr_vector varying the block size (sdsl::rrr_vector)\n"
                     " - la_vector varying epsilon (la_vector<c>)\n"
                     " - la_vector space optimized (la_vector_opt)\n"
                     " - gap+encoding + elias delta/elias gamma varying the sample density (sdsl::enc_vector)\n"
                     " - Partitioned Elias Fano (uniform partition ds2i::uniform_partitioned_sequence)\n"
                     " - Partitioned Elias Fano (optimal partition ds2i::partitioned_sequence)"
                     << std::endl;
        return -1;
    }

    constexpr bool la_vector_epsilon = true;
    constexpr bool la_vector_opt = true;
    constexpr bool elias_codes = true;
    constexpr bool ds2i = true;

    using type = uint32_t;

    std::cout << get_csv_header<la_vector_epsilon, la_vector_opt, elias_codes, ds2i>() << std::endl;
    for (int i = 1; i < argc; ++i) {
        std::string filename(argv[i]);
        std::vector<type> data(read_data_binary<type, type>(filename, true, (1u << 31) - 1));
        assert(std::is_sorted(data.begin(), data.end()));
        type u = data.back() + 1;
        type n = data.size();
        double ratio = static_cast<double>(n) / u;
        double mean = 0;
        double stdev = 0;
        {
            std::vector<type> data_gaps;
            data_gaps.reserve(n);
            for (auto j = 0; j < n - 1; ++j) {
                data_gaps.push_back(data[j + 1] - data[j]);
            }
            uint64_t sum = std::accumulate(data_gaps.begin(), data_gaps.end(), (uint64_t) 0);
            mean = static_cast<double>(sum) / data_gaps.size();
            double sq_sum = std::inner_product(data_gaps.begin(), data_gaps.end(), data_gaps.begin(), 0.0);
            stdev = std::sqrt(sq_sum / data_gaps.size() - mean * mean);
            data_gaps.clear();
        }
        std::cout << filename << "," << n << "," << u << "," << ratio << "," << mean << "," << stdev << std::flush;
        test_all<la_vector_epsilon, la_vector_opt, elias_codes, ds2i>(data);
        std::cout << std::endl;
    }
    return 0;
}
