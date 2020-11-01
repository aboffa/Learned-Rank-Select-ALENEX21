#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <climits>

struct WT_char_positions {
    char c;
    std::vector<uint32_t> *vec;

    WT_char_positions() {}

    WT_char_positions(std::vector<uint32_t> *_vec, char _c) {
        c = _c;
        vec = _vec;
    }
};

template<typename Dist, typename Gen>
std::vector<typename Dist::result_type> generate_data(Dist &distribution, Gen &generator, size_t n, bool sorted) {
    using T = typename Dist::result_type;
    std::vector<T> out(n);
    std::generate(out.begin(), out.end(), [&]() { return distribution(generator); });
    if (sorted)
        std::sort(out.begin(), out.end());
    return out;
}

template<typename Dist, typename Gen>
std::vector<typename Dist::result_type> generate_unique(Dist &distribution, Gen &generator, size_t n, bool sorted) {
    using T = typename Dist::result_type;

    if constexpr (std::is_same<Dist, std::uniform_int_distribution<T>>::value) {
        std::vector<T> out(n);
        size_t i = 0;
        size_t u = distribution.max() - distribution.min();
        for (auto k = 0; k < u && i < n; ++k)
            if (generator() % (u - k) < n - i)
                out[i++] = k + distribution.min();

        if (!sorted)
            std::random_shuffle(out.begin(), out.end());

        return out;
    }

    std::unordered_set<T> set;
    set.reserve(n);

    while (set.size() < n)
        set.insert(distribution(generator));
    std::vector<T> out(set.begin(), set.end());

    if (sorted)
        std::sort(out.begin(), out.end());
    return out;
}

template<typename Dist, typename Gen>
std::vector<typename Dist::result_type> generate_gap_data(Dist &distribution, Gen &generator, size_t n) {
    using T = typename Dist::result_type;
    std::vector<T> out;
    out.reserve(n);
    out.push_back(0);
    Dist d;

    for (auto i = 0; i < n - 1; ++i)
        out.push_back(out.back() + distribution(generator));

    return out;
}


struct datasetStats {
    float averageLCP = 0;
    float averageLength = 0;
    int maxSize = 0;
    int minSize = INT_MAX;
    int size = 0;
    std::vector<char> chars;
};

bool checkIsPrintable(std::string s) {
    for (char c : s) {
        if (!isprint(c)) {
            return false;
        }
    }
    return true;
}

bool checkIsNotControl(std::string s) {
    for (char c : s) {
        if ((int) c < 31 || (int) c > 255) {
            return false;
        }
    }
    return true;
}

float lcp(std::string a, std::string b) {
    float result = 0;
    for (int i = 0; i < std::min(a.size(), b.size()); i++) {
        if (a[i] == b[i]) {
            result++;
        }
    }
    return result;
}

// Merging all the 5gram of the files in filenames
template<bool duplicate>
datasetStats *
load_data_from_Google_nGram_v1(std::vector<std::string> &data, std::string basePath,
                               std::vector<std::string> filenames) {
    auto ds = new datasetStats();
    int i = 0;
    for (const std::string &filename : filenames) {
        std::ifstream in(basePath + filename, std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Unable to open " << basePath + filename << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string delimiter = "\t";
        std::string token;
        for (std::string line; getline(in, line);) {
            token = line.substr(0, line.find(delimiter));
            int len = token.length();
            // deleting whitespaces in front and at the end of the ngram
            while (token[len - 1] == ' ') { token.erase(--len, 1); }
            while (token[0] == ' ') { token.erase(0, 1); }
            bool check = !token.empty() && checkIsNotControl(token);
            if constexpr (!duplicate){
                check = check && (data.empty() || data.back() != token);
            }
            if ( check ) {
                // to cheching for now if (data.empty() || data.back() != token) {
                if (!data.empty()) {
                    ds->averageLCP += lcp(data.back(), token);
                }
                for (auto c : token) {
                    auto foo = std::find(ds->chars.begin(), ds->chars.end(), c);
                    if (foo == std::end(ds->chars)) {
                        ds->chars.push_back(c);
                    }
                }
                data.push_back(token);
                ds->averageLength += token.size();
                if (token.size() > ds->maxSize) {
                    ds->maxSize = token.size();
                }
                if (token.size() < ds->minSize) {
                    ds->minSize = token.size();
                }
            }
        }
        i++;
        std::sort(data.begin(), data.end());
        if constexpr (!duplicate) {
            data.erase(unique(data.begin(), data.end()), data.end());
        }
        std::cout << i << "/" << filenames.size() << " " << filename << " done! Data size  = " << data.size()
                  << std::endl;
    }
    ds->size = data.size();
    ds->averageLength /= static_cast<float>(ds->size);
    ds->averageLCP /= static_cast<float>(ds->size);
    return ds;
}