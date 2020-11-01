//
// Created by boffa on 14/10/20.
//

/*! \file enc_vector.hpp
   \brief enc_vector.hpp contains the sdsl::enc_vector class.
   \author Simon Gog
*/

#include "sdsl/int_vector.hpp"
#include "sdsl/coder.hpp"
#include "sdsl/iterators.hpp"


//! A generic immutable space-saving vector class for unsigned integers.
/*! A vector v is stored more space-efficiently by self-delimiting coding
 *  the deltas v[i+1]-v[i] (v[-1]:=0). Space of the structure and random
 *  access time to it can be controlled by a sampling parameter t_dens.
 *
 *  \tparam t_coder  Self-delimiting coder.
 *  \tparam t_dens   Every t_dens-th element of v is sampled.
 *  \tparam t_width  Width of the int_vector used to store the samples and pointers.
 *  This class is a parameter of csa_sada.
 * @ingroup int_vector
 */

namespace sdsl {

    template<class t_coder=coder::elias_delta,
            uint32_t t_dens = 128, uint8_t t_width = 0>
    class enc_vector_rank {
    private:
        static_assert(t_dens > 1, "enc_vector: sample density must be larger than `1`");
    public:
        typedef uint64_t value_type;
        typedef random_access_const_iterator <enc_vector_rank> iterator;
        typedef iterator const_iterator;
        typedef const value_type reference;
        typedef const value_type const_reference;
        typedef const value_type *const_pointer;
        typedef ptrdiff_t difference_type;
        typedef int_vector<>::size_type size_type;
        typedef t_coder coder;
        typedef typename enc_vector_trait<t_width>::int_vector_type int_vector_type;
        typedef iv_tag index_category;
        static const uint32_t sample_dens = t_dens;
        typedef enc_vector_rank enc_vec_type;

        int_vector<0> m_z;                       // storage for encoded deltas
        int_vector_type m_sample_vals_and_pointer; // samples and pointers
        size_type m_size = 0;                // number of vector elements
    private:


        void clear() {
            m_z.resize(0);
            m_size = 0;
            m_sample_vals_and_pointer.resize(0);
        }

    public:
        enc_vector_rank() = default;

        enc_vector_rank(const enc_vector_rank &) = default;

        enc_vector_rank(enc_vector_rank &&) = default;

        enc_vector_rank &operator=(const enc_vector_rank &) = default;

        enc_vector_rank &operator=(enc_vector_rank &&) = default;

        //! Constructor for a Container of unsigned integers.
        /*! \param c A container of unsigned integers.
          */
        template<class Container>
        enc_vector_rank(const Container &c);

        //! Constructor for an int_vector_buffer of unsigned integers.
        /*
            \param v_buf A int_vector_buf.
        */
        template<uint8_t int_width>
        enc_vector_rank(sdsl::int_vector_buffer<int_width> &v_buf);

        //! Default Destructor
        ~enc_vector_rank() {}

        //! The number of elements in the enc_vector.
        size_type size() const {
            return m_size;
        }

        //! Return the largest size that this container can ever have.
        static size_type max_size() {
            return int_vector<>::max_size() / 2;
        }

        //!    Returns if the enc_vector is empty.
        bool empty() const {
            return 0 == m_size;
        }

        //! Swap method for enc_vector
        void swap(enc_vector_rank &v);

        //! Iterator that points to the first element of the enc_vector.
        const const_iterator begin() const {
            return const_iterator(this, 0);
        }

        //! Iterator that points to the position after the last element of the enc_vector.
        const const_iterator end() const {
            return const_iterator(this, this->m_size);
        }

        //! operator[]
        /*! \param i Index. \f$ i \in [0..size()-1]\f$.
         */
        value_type operator[](size_type i) const;

        //! Serialize the enc_vector to a stream.
        /*! \param out Out stream to write the data structure.
            \return The number of written bytes.
         */
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr, std::string name = "") const;

        //! Load the enc_vector from a stream.
        void load(std::istream &in);

        //! Returns the i-th sample of enc_vector
        /*! \param i The index of the sample. 0 <= i < size()/get_sample_dens()
         *  \return The value of the i-th sample.
         */
        value_type sample(const size_type i) const;

        uint32_t get_sample_dens() const {
            return t_dens;
        }

        /*!
         * \param i The index of the sample for which all values till the next sample should be decoded. 0 <= i < size()/get_sample_dens()
         * \param it A pointer to a uint64_t vector, whereto the values should be written
         */
        void get_inter_sampled_values(const size_type i, uint64_t *it) const {
            *(it++) = 0;
            if (i * t_dens + t_dens - 1 < size()) {
                t_coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i << 1) + 1], t_dens - 1,
                                                     it);
            } else {
                assert(i * t_dens < size());
                t_coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i << 1) + 1],
                                                     size() - i * t_dens - 1, it);
            }
        };
    };

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    inline typename enc_vector_rank<t_coder, t_dens, t_width>::value_type
    enc_vector_rank<t_coder, t_dens, t_width>::operator[](const size_type i) const {
        assert(i + 1 != 0);
        assert(i < m_size);
        size_type idx = i / get_sample_dens();
        return m_sample_vals_and_pointer[idx << 1] +
               t_coder::decode_prefix_sum(m_z.data(), m_sample_vals_and_pointer[(idx << 1) + 1], i - t_dens * idx);
    }

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    inline typename enc_vector_rank<t_coder, t_dens, t_width>::value_type
    enc_vector_rank<t_coder, t_dens, t_width>::sample(const size_type i) const {
        assert(i * get_sample_dens() + 1 != 0);
        assert(i * get_sample_dens() < m_size);
        return m_sample_vals_and_pointer[i << 1];
    }

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    void enc_vector_rank<t_coder, t_dens, t_width>::swap(enc_vector_rank<t_coder, t_dens, t_width> &v) {
        if (this != &v) { // if v and _this_ are not the same object
            m_z.swap(v.m_z);
            m_sample_vals_and_pointer.swap(v.m_sample_vals_and_pointer);
            std::swap(m_size, v.m_size);
        }
    }

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    template<class Container>
    enc_vector_rank<t_coder, t_dens, t_width>::enc_vector_rank(const Container &c) {
        // clear bit_vectors
        clear();

        if (c.empty())  // if c is empty there is nothing to do...
            return;
        typename Container::const_iterator it = c.begin(), end = c.end();
        typename Container::value_type v1 = *it, v2, max_sample_value = 0, x;
        size_type samples = 0;
        size_type z_size = 0;
//  (1) Calculate maximal value of samples and of deltas
        for (size_type i = 0, no_sample = 0; it != end; ++it, ++i, --no_sample) {
            v2 = *it;
            if (!no_sample) { // add a sample
                no_sample = get_sample_dens();
                if (max_sample_value < v2) max_sample_value = v2;
                ++samples;
            } else {
                z_size += t_coder::encoding_length(v2 - v1);
            }
            v1 = v2;
        }
//    (2) Write sample values and deltas
        {
            if (max_sample_value > z_size + 1)
                m_sample_vals_and_pointer.width(bits::hi(max_sample_value) + 1);
            else
                m_sample_vals_and_pointer.width(bits::hi(z_size + 1) + 1);
            m_sample_vals_and_pointer.resize(2 * samples + 2); // add 2 for last entry
            util::set_to_value(m_sample_vals_and_pointer, 0);
            typename int_vector_type::iterator sv_it = m_sample_vals_and_pointer.begin();
            z_size = 0;
            size_type no_sample = 0;
            for (it = c.begin(); it != end; ++it, --no_sample) {
                v2 = *it;
                if (!no_sample) { // add a sample
                    no_sample = get_sample_dens();
                    *sv_it = v2;
                    ++sv_it;
                    *sv_it = z_size;
                    ++sv_it;
                } else {
                    x = v2 - v1;
                    z_size += t_coder::encoding_length(x);
                }
                v1 = v2;
            }
            *sv_it = 0;
            ++sv_it;        // initialize
            *sv_it = z_size + 1;
            ++sv_it; // last entry

            m_z = int_vector<>(z_size, 0, 1);
            uint64_t *z_data = t_coder::raw_data(m_z);
            uint8_t offset = 0;
            no_sample = 0;
            for (it = c.begin(); it != end; ++it, --no_sample) {
                v2 = *it;
                if (!no_sample) { // add a sample
                    no_sample = get_sample_dens();
                } else {
                    t_coder::encode(v2 - v1, z_data, offset);
                }
                v1 = v2;
            }
        }
        m_size = c.size();
    }

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    template<uint8_t int_width>
    enc_vector_rank<t_coder, t_dens, t_width>::enc_vector_rank(int_vector_buffer <int_width> &v_buf) {
        // clear bit_vectors
        clear();
        size_type n = v_buf.size();
        if (n == 0)  // if c is empty there is nothing to do...
            return;
        value_type v1 = 0, v2 = 0, max_sample_value = 0;
        size_type samples = 0, z_size = 0;
        const size_type sd = get_sample_dens();
//  (1) Calculate maximal value of samples and of deltas
        for (size_type i = 0, no_sample = 0; i < n; ++i, --no_sample) {
            v2 = v_buf[i];
            if (!no_sample) { // is sample
                no_sample = sd;
                if (max_sample_value < v2) max_sample_value = v2;
                ++samples;
            } else {
                z_size += t_coder::encoding_length(v2 - v1);
            }
            v1 = v2;
        }

//    (2) Write sample values and deltas
//    (a) Initialize array for sample values and pointers
        if (max_sample_value > z_size + 1)
            m_sample_vals_and_pointer.width(bits::hi(max_sample_value) + 1);
        else
            m_sample_vals_and_pointer.width(bits::hi(z_size + 1) + 1);
        m_sample_vals_and_pointer.resize(2 * samples + 2); // add 2 for last entry
        util::set_to_value(m_sample_vals_and_pointer, 0);

//    (b) Initilize bit_vector for encoded data
        m_z = int_vector<>(z_size, 0, 1);
        uint64_t *z_data = t_coder::raw_data(m_z);
        uint8_t offset = 0;

//    (c) Write sample values and deltas
        z_size = 0;
        for (size_type i = 0, j = 0, no_sample = 0; i < n; ++i, --no_sample) {
            v2 = v_buf[i];
            if (!no_sample) { // is sample
                no_sample = sd;
                m_sample_vals_and_pointer[j++] = v2;    // write samples
                m_sample_vals_and_pointer[j++] = z_size;// write pointers
            } else {
                z_size += t_coder::encoding_length(v2 - v1);
                t_coder::encode(v2 - v1, z_data, offset);   // write encoded values
            }
            v1 = v2;
        }
        m_size = n;
    }

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    enc_vector_rank<>::size_type
    enc_vector_rank<t_coder, t_dens, t_width>::serialize(std::ostream &out, structure_tree_node *v,
                                                         std::string name) const {
        structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_size, out, child, "size");
        written_bytes += m_z.serialize(out, child, "encoded deltas");
        written_bytes += m_sample_vals_and_pointer.serialize(out, child, "samples_and_pointers");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    template<class t_coder, uint32_t t_dens, uint8_t t_width>
    void enc_vector_rank<t_coder, t_dens, t_width>::load(std::istream &in) {
        read_member(m_size, in);
        m_z.load(in);
        m_sample_vals_and_pointer.load(in);
    }

    namespace coder {

        class elias_delta_rank : public elias_delta {

        public:
            template<typename T>
            static size_t
            decode_until_succ(const uint64_t *d, const size_type start_idx, size_type n, const uint64_t sample,
                              const T val) {
                size_type i = 0;
                size_t to_return = 0;
                uint64_t value = sample;
                d += (start_idx >> 6);
                uint8_t offset = start_idx & 0x3F;
                size_type len_1_len, len;
                while (i++ < n) {
                    if (value >= val)
                        return to_return;
                    to_return++;
                    len_1_len = bits::read_unary_and_move(d, offset);
                    auto gap = 0;
                    if (!len_1_len) {
                        gap = 1;
                    } else {
                        len = bits::read_int_and_move(d, offset, len_1_len) + (1ULL << len_1_len);
                        gap += bits::read_int_and_move(d, offset, len - 1) + (len - 1 < 64) * (1ULL << (len - 1));
                    }
                    value += gap;
                }
                if (value >= val)
                    return to_return;
                return ++to_return;
            }
        };

        class elias_gamma_rank : public elias_gamma {

        public:
            template<typename T>
            static size_t
            decode_until_succ(const uint64_t *d, const size_type start_idx, size_type n, const uint64_t sample,
                              const T val) {
                size_type i = 0;
                size_t to_return = 0;
                uint64_t value = sample;
                d += (start_idx >> 6);
                uint8_t offset = start_idx & 0x3F;
                while (i++ < n) {
                    if (value >= val)
                        return to_return;
                    to_return++;
                    uint16_t len_1 = bits::read_unary_and_move(d, offset); // read length of x-1
                    auto gap = 0;
                    if (!len_1) {
                        gap = 1;
                    } else {
                        gap = bits::read_int_and_move(d, offset, len_1) + (len_1 < 64) * (1ULL << len_1);
                    }
                    value += gap;
                }
                if (value >= val)
                    return to_return;
                return ++to_return;
            }

        };
    }
}

template<class ForwardIterator, class T>
size_t even_lower_bound(ForwardIterator first, ForwardIterator last, const T &val) {
    ForwardIterator base = first;
    size_t n = distance(first, last);
    while (n > 2) {
        size_t half = n / 2;
        half &= ~(1UL << 0);
        base = (base[half] < val) ? base + half : base;
        n -= half;
    }
    auto to_return = (*base < val) + base - first;
    to_return &= ~(1UL << 0);
    return to_return;
}

template<typename V, typename D>
struct wrapper_sdsl_enc_vector : public V {
    typedef typename D::value_type T;

    wrapper_sdsl_enc_vector(const D &data) : V(data) {}

    inline double bits_per_element() const {
        return sdsl::size_in_bytes(*this) * 8. / this->size();
    }

    inline T select(size_t i) const { return this->operator[](i-1); }

    inline size_t rank(T x) const {
        // check in which block to search. (successor block)
        size_t index_sample = even_lower_bound(this->m_sample_vals_and_pointer.begin(),
                                               this->m_sample_vals_and_pointer.end() - 2, x);
        return ((index_sample / 2) * this->get_sample_dens()) +
               V::coder::decode_until_succ(this->m_z.data(),
                                           this->m_sample_vals_and_pointer[index_sample + 1],
                                           this->get_sample_dens() - 1,
                                           this->m_sample_vals_and_pointer[index_sample],
                                           x);
    }

};

template<typename T, typename D>
class wrapper_sdsl {
protected:
    D ef;
    typename D::select_1_type ef_select;
    typename D::rank_1_type ef_rank;
    size_t n;

public:
    wrapper_sdsl(std::vector<T> &a) {
        n = a.size();
        sdsl::bit_vector bv(a.back() + 1, 0);
        for (auto i = 0; i < a.size(); ++i) {
            bv[a[i]] = 1;
        }
        ef = *new D(bv);
        sdsl::util::init_support(ef_select, &ef);
        sdsl::util::init_support(ef_rank, &ef);
        return;
    }

    inline T select(size_t i) const {
        assert(i > 0 && i <= n);
        return ef_select(i);
    }

    inline double bits_per_element() const {
        double size = (double) (sdsl::size_in_bytes(ef) * 8);
        size += (double) (sdsl::size_in_bytes(ef_select) * 8);
        size += (double) (sdsl::size_in_bytes(ef_rank) * 8);
        size /= static_cast<double>(n);
        return size;
    };

    inline size_t rank(T x) const {
        return ef_rank.rank(x);
    }

    inline void decode(T *out) const {
        for (size_t i = 0; i < n; ++i) {
            *(out++) = ef_select(i + 1);
        }
    }

    inline std::vector<T> decode() const {
        std::vector<T> out(n);
        decode(out.data());
        return out;
    }

};

