/* sdsl - succinct data structures library
    Copyright (C) 2015 Dorle Osterode

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file wt_dna.hpp
    \brief wt_dna.hpp contains a specialized class to support select, rank
			and access on inputs over a dna alphabet (A,C,G,T).
    \author Dorle Osterode
*/
#ifndef INCLUDED_SDSL_WT_DNA_IL
#define INCLUDED_SDSL_WT_DNA_IL

#include <iostream> // for stdout
#include <cstdlib>
#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

// TODO: write a proper intro
//! A wavelet tree class for integer sequences.
/*!
 *  \tparam t_rac         Type of the random access container used for storing the permutation.
 *  \tparam t_inv_support Type of the support structure for inverse permutation
 *  \tparam t_bitvector   Type of the bitvector used for storing B and X.
 *  \tparam t_select      Type of the support structure for select on pattern `1`.
 *  \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 *
 *   @ingroup wt
 */
template<int t_exponent = 7> // exponent for the interval to count occurences in.
class wt_dna_il
{
        static_assert(t_exponent > 4,
                      "First template argument has to be greater then 4.");

    public:
        // need for typedefs of size_type and value_type
        typedef uint64_t size_type; // is this correct?
        typedef char value_type; // is this correct?
        typedef wt_tag index_category;
        typedef random_access_const_iterator<wt_dna_il> const_iterator;
        typedef byte_alphabet_tag alphabet_category;

    private:
        // the private variables
        int_vector<64> m_encseq; // encoded bwt with interleaved absolute counts
        uint64_t m_size; // length of input sequence
        uint8_t m_exponent = t_exponent; // exponent of the interval size
        uint64_t m_interval = 1 << m_exponent; // interval for which the counts are stored
        uint64_t m_mod_interval_mask = m_interval - 1; // mask to calculate modulo the interval

        // this is a counting function from BWA-mem bwt.c:98
        uint64_t count_c(uint64_t y, int c) const
        {
            // reduce nucleotide counting to bits counting
            uint64_t ret = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
            // count the number of 1s in y
            ret = (ret & 0x3333333333333333ull) + (ret >> 2 & 0x3333333333333333ull);
            return ((ret + (ret >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
        }

        // this is a modified counting function from BWA-mem bwt.c:98. it counts the pattern c
        // upto index idx. the index is enumerated from left to right starting with 0.
        uint64_t count_c_upto(uint64_t y, int c, uint64_t idx) const
        {
            uint64_t ret = y & ~((1ull<<((~idx&31)<<1)) - 1);
            ret = ((c&2)? ret : ~ret) >> 1 & ((c&1)? ret : ~ret) & 0x5555555555555555ull;
            // count the number of 1s in y
            ret = (ret & 0x3333333333333333ull) + (ret >> 2 & 0x3333333333333333ull);
            return (((ret + (ret >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56) - ((c == 0)? (~idx&31): 0);

        }

    public:

        //! Default constructor
        wt_dna_il() {}

        //! Semi-external constructor
        /*! \param buf         File buffer of a vector containing the already encoded
         *                     sequence for which the wt_dna_il should be build.
         *  \param size        Size of the prefix of input, which should be indexed.
         */
        template<uint8_t int_width>
        wt_dna_il(int_vector_buffer<int_width>& input, size_type size) : m_size(size)
        {
            // if m_size % m_interval == 0: this equals m_size/m_interval + 1
            // else: this is m_size/m_interval + 2.
            // we want to store the counts for the whole sequence
            uint64_t n_occ = (m_size + m_interval - 1) / m_interval + 1;
            uint64_t enc_size = (m_size / 32) + ((m_size % 31 != 0)? 1: 0) + (n_occ * 4);

            // create a new int_vector<64> to store the encoded sequence and the absolute counts
            int_vector<64> tmp(enc_size);
            uint64_t occ[4] = {0, 0, 0, 0};

            // iterate over the input sequence and store s[i]
            size_type j = 0;
            value_type c;
            size_type k = 0;
            uint64_t encoded = 0;

            // TODO: change iterativ array copying to memcpy.
            // store the first counts {0, 0, 0, 0}
            tmp[k] = occ[0];
            tmp[k + 1] = occ[1];
            tmp[k + 2] = occ[2];
            tmp[k + 3] = occ[3];
            k += 4;

            for (size_type i = 0; i < m_size; i++) {
                c = input[i];
                if (c >= 4) {
                    // got a non-valid character!!! error!!!1elf
                    // TODO: handle this case
                    std::cout << "wrong char: " << input[i]  << " at position: " << (c & 3ULL) << "\n";
                    assert(!"reached");
                } else {
                    if ((i & 31) == 0 && i > 0) {
                        // store the prepared uint64 with 32 chars in it
                        tmp[k] = encoded;
                        k += 1;
                        encoded = 0;
                    }
                    if ((i & m_mod_interval_mask) == 0 && i > 0) {
                        // store the counts if i % m_interval == 0
                        tmp[k] = occ[0];
                        tmp[k + 1] = occ[1];
                        tmp[k + 2] = occ[2];
                        tmp[k + 3] = occ[3];
                        k += 4;
                    }

                    // count the character and prepare the uint64
                    // encoded |= ((c & 3ULL) << (62 - (2 * (i % 32))));
                    encoded |= ((uint64_t)c << (62 - ((i & 31) << 1)));
                    occ[(int)c] += 1;
                }
            }

            // store the last encoded and occs. this has to happen every time.
            // if m_size % mod_intervall_mask == 0, then the last encoded and occurence is not
            // stored yet. the for-loop iterates until i == m_size - 1.
            // if m_size - 1 % mod_intervall_mask == 0, then the last char and occurence is not
            // stored yet. the loop breaks before it is stored.
            // in any other case the last encoded and occs are not stored.
            tmp[k] = encoded;
            tmp[k + 1] = occ[0];
            tmp[k + 2] = occ[1];
            tmp[k + 3] = occ[2];
            tmp[k + 4] = occ[3];

            // swap the created encseq into the private variable
            m_encseq.swap(tmp);
        }

        // just needed for debugging. can be removed later
        void print_encseq()
        {
            for (unsigned int i = 0; i < m_size; i++) {
                if (i % 128 == 0) {
                    std::cout << "occ[a]: " << m_encseq[i] << "\n";
                    std::cout << "occ[c]: " << m_encseq[i + 1] << "\n";
                    std::cout << "occ[g]: " << m_encseq[i + 2] << "\n";
                    std::cout << "occ[t]: " << m_encseq[i + 3] << "\n";
                    i += 4;
                }
                size_type index = (i / 32) + ((i / 128) * 4) + 1;
                std::cout << ((m_encseq[index] >> (62 - (2 * (index % 32)))) & 3ULL);
            }
            std::cout << "\n";
        }

        //! Copy constructor
        wt_dna_il(const wt_dna_il& wt)
        {
            m_encseq = wt.m_encseq;
            m_size = wt.m_size;
        }

        //! Assignment operator
        wt_dna_il& operator=(const wt_dna_il& wt)
        {
            wt_dna_il tmp(wt);
            tmp.swap(*this);
            return *this;
        }

        //! Swap operator
        void swap(wt_dna_il& fs)
        {
            if (this != &fs) {
                m_encseq.swap(fs.m_encseq);
                std::swap(m_size, fs.m_size);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const
        {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const
        {
            return m_size == 0;
        }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const
        {
            assert(i < size());
            // size_type index = (i / 32) + (((i / 128) + 1) * 4);
            size_type index = (i >> 5) + (((i >> m_exponent) + 1) << 2);
            // return (m_encseq[index] >> (62 - (2 * (i % 32)))) & 3ULL;
            return (m_encseq[index] >> (62 - ((i & 31) << 1))) & 3ULL;
        }

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const
        {
            assert(c < 4);
            assert(i <= size());

            uint64_t occ_pos = (i >> m_exponent) * ((m_interval/32) + 4);
            // std::cout << "occ_pos: " << occ_pos << "\n";
            uint64_t sum = m_encseq[occ_pos + c];
            // std::cout << "first sum: " << sum << "\n";
            uint64_t todo = i & m_mod_interval_mask;
            // std::cout << "bases to count: " << todo << "\n";
            occ_pos += 4;
            while (todo > 32) {
                // std::cout << "counting in word: " << m_encseq[occ_pos] << "\n";
                sum += count_c(m_encseq[occ_pos], c);
                occ_pos += 1;
                todo -= 32;
            }
            // std::cout << "sum after counting complete words: " << sum << "\n";
            // std::cout << "bases to count: " << todo << "\n";
            if (todo > 0) {
                // std::cout << "counting in word: " << m_encseq[occ_pos] << "\n";
                sum += count_c_upto(m_encseq[occ_pos], c, todo - 1);
            }
            // std::cout << "actual return value: " << sum << "\n";
            return sum;
        }

        //! Calculates how many occurrences of symbol input[i] are in the prefix [0..i-1] of the original input.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(input[i],i), input[i])
         *  \par Time complexity
         *       \f$ \Order{1} + One access to the inverse permutation \f$
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type> inverse_select(size_type i)const
        {
            assert(i < size());

            size_type index = (i >> 5) + (((i >> m_exponent) + 1) << 2);
            uint8_t enc_char = (m_encseq[index] >> (62 - ((i & 31) << 1))) & 3ULL;

            return std::make_pair(rank(i, enc_char), enc_char);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{1} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const
        {
            assert(1 <= i and i <= rank(size(), c));
            // TODO: select is not implemented! for the actual use-case select is not needed
            assert(!"reached");
            return 0;
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const
        {
            return const_iterator(this, size());
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += m_encseq.serialize(out, child, "encseq");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            m_encseq.load(in);
        }
};
}

#endif
