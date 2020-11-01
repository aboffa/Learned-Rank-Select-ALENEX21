#pragma once

#include <cstdint>
#include <fstream>
#include <exception>

template<typename TypeIn, typename TypeOut>
std::vector<TypeOut> read_data_binary(const std::string &filename, bool first_is_size = true, size_t max_size = std::numeric_limits<TypeIn>::max()) {
    try {
        auto openmode = std::ios::in | std::ios::binary;
        if (!first_is_size)
            openmode |= std::ios::ate;

        std::fstream in(filename, openmode);
        in.exceptions(std::ios::failbit | std::ios::badbit);

        size_t size;
        if (first_is_size)
            in.read((char *) &size, sizeof(size_t));
        else {
            size = static_cast<size_t>(in.tellg() / sizeof(TypeIn));
            in.seekg(0);
        }
        size = std::min(max_size, size);

        std::vector<TypeIn> data(size);
        in.read((char *) data.data(), size * sizeof(TypeIn));
        if constexpr (std::is_same<TypeIn, TypeOut>::value)
            return data;

        return std::vector<TypeOut>(data.begin(), data.end());
    }
    catch (std::ios_base::failure &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << std::strerror(errno) << std::endl;
        exit(1);
    }
}

/**
 * Reads binary files containing integers lists coded as 32-bit little-endian unsigned integers. Specifically, each list
 * starts with an integer indicating its length followed by the corresponding number of integers in sorted order.
 *
 * Compatible both with the
 * (i) ds2i format (the first list has length 1 and the only integer is the number of documents in the collection);
 * (ii) Lemire/Boytsov format (there is not a singleton list at the beginning of the file).
 */
class IndexReader {
    std::ifstream in;
    size_t current_byte;
    size_t file_size;

    inline uint32_t read_uint32() {
        uint32_t n = 0;
        in.read(reinterpret_cast<char *>(&n), sizeof(n));
        return n;
    }

public:

    explicit IndexReader(std::string &file) : in(), current_byte(0) {
        in.exceptions(std::ios::failbit | std::ios::badbit);
        in.open(file, std::ios::binary);
        auto fs = in.tellg();
        in.seekg(0, std::ios::end);
        file_size = size_t(in.tellg() - fs);
        in.seekg(0, std::ios::beg);
        in.exceptions(std::ios::goodbit);

        uint32_t qty;
        std::ifstream::pos_type pos;
        read_next_pos_qty(pos, qty);
        if (qty != 1) // likely in Lemire/Boytsov format
            in.seekg(0, std::ios::beg);
    }

    bool read_next_pos_qty(std::ifstream::pos_type &pos, uint32_t &qty) {
        pos = in.tellg();
        qty = read_uint32();
        if (in.eof())
            return false;
        in.seekg(qty * sizeof(uint32_t), std::ifstream::cur);
        current_byte += (qty + 1) * sizeof(uint32_t);
        return true;
    }

    void seek(std::ifstream::pos_type pos) {
        in.clear();
        in.seekg(pos);
    }

    bool load_integers(std::vector<uint32_t> &buffer) {
        uint32_t qty = read_uint32();
        if (in.eof())
            return false;
        buffer.resize(qty);
        in.read(reinterpret_cast<char *>(buffer.data()), qty * sizeof(uint32_t));
        current_byte += (qty + 1) * sizeof(uint32_t);
        return true;
//        for (uint32_t i = 0; i < qty; ++i)
//            buffer[i] = read_uint32();
    }

    size_t get_progress() const {
        return current_byte;
    }

    size_t get_file_size() const {
        return file_size;
    }

};