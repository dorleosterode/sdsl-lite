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
                      and access on inputs over an alphabet of size four, e.g.
		      (A,C,G,T).
    \author Dorle Osterode
*/
#ifndef INCLUDED_SDSL_WT_DNA
#define INCLUDED_SDSL_WT_DNA

#include <iostream> // needed for std::out
#include <cstdlib>
#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A wavelet tree class specialized for alphabets of size four.
/*!
 *  \tparam t_bitvector         Type of the bitvector used for storing the encoded sequence.
 *
 *   @ingroup wt
 */
template<class t_bitvector = bit_vector>
class wt_dna
{
    private:
        // typedefs for the used support structures
        typedef rank_support_v5<SDSL_WT_DNA_TYPE_A, 2> rank_00_type;
        typedef rank_support_v5<SDSL_WT_DNA_TYPE_C, 2> rank_01_type;
        typedef rank_support_v5<SDSL_WT_DNA_TYPE_G, 2> rank_10_type;
        typedef rank_support_v5<SDSL_WT_DNA_TYPE_T, 2> rank_11_type;

        typedef select_support_scan<SDSL_WT_DNA_TYPE_A, 2> select_00_type;
        typedef select_support_scan<SDSL_WT_DNA_TYPE_C, 2> select_01_type;
        typedef select_support_scan<SDSL_WT_DNA_TYPE_G, 2> select_10_type;
        typedef select_support_scan<SDSL_WT_DNA_TYPE_T, 2> select_11_type;

    public:
        // need for typedefs of size_type and value_type
        typedef uint64_t size_type;
        typedef char value_type;
        typedef wt_tag index_category;
        typedef random_access_const_iterator<wt_dna> const_iterator;
        typedef byte_alphabet_tag alphabet_category;

    private:
        // the private variables
        t_bitvector m_encseq; // encoded bwt
        uint64_t m_size; // length of input sequence

        // rank data structures needed for patterns 00 01 10 11
        rank_00_type m_rank_00;
        rank_01_type m_rank_01;
        rank_10_type m_rank_10;
        rank_11_type m_rank_11;

        // select data structures needed for patterns 00 01 10 11
        select_00_type m_select_00;
        select_01_type m_select_01;
        select_10_type m_select_10;
        select_11_type m_select_11;

    public:

        //! Default constructor
        wt_dna() {}

        //! Semi-external constructor
        /*! \param buf         File buffer of a vector containing the already encoded
         *                     sequence for which the wt_dna should be build.
         *  \param size        Size of the prefix of input, which should be indexed.
         */
        template<uint8_t int_width>
        wt_dna(int_vector_buffer<int_width>& input, size_type size) : m_size(size)
        {
            // create a new bit_vector to store the encoded sequence
            uint64_t enc_size = m_size << 1;
            t_bitvector tmp(enc_size, 0);

            // iterate over the input sequence and store s[i] at b[2i]b[2i+1]
            size_type j = 0;
            value_type c;
            for (size_type i = 0; i < m_size; i++) {
                c = input[i];
                if (c >= 4) {
                    // got a non-valid character
                    std::cout << "wrong char: " << input[i]  << " at position: " << i << "\n";
                    assert(!"reached");
                } else {
                    j = i << 1;
                    tmp[j] = (c >> 1) & 1ULL;
                    tmp[j + 1] = c & 1ULL;
                }
            }

            // swap the created encseq into the private variable
            m_encseq.swap(tmp);

            // init support of the rank structures
            util::init_support(m_rank_00, &m_encseq);
            util::init_support(m_rank_01, &m_encseq);
            util::init_support(m_rank_10, &m_encseq);
            util::init_support(m_rank_11, &m_encseq);

            // init support of the select structures
            util::init_support(m_select_00, &m_encseq);
            util::init_support(m_select_01, &m_encseq);
            util::init_support(m_select_10, &m_encseq);
            util::init_support(m_select_11, &m_encseq);
        }

        //! Copy constructor
        wt_dna(const wt_dna& wt)
        {
            m_encseq = wt.m_encseq;
            m_size = wt.m_size;

            // copy rank support
            m_rank_00 = wt.m_rank_00;
            m_rank_00.set_vector(&m_encseq);
            m_rank_01 = wt.m_rank_01;
            m_rank_01.set_vector(&m_encseq);
            m_rank_10 = wt.m_rank_10;
            m_rank_10.set_vector(&m_encseq);
            m_rank_11 = wt.m_rank_11;
            m_rank_11.set_vector(&m_encseq);

            // copy select support
            m_select_00 = wt.m_select_00;
            m_select_00.set_vector(&m_encseq);
            m_select_01 = wt.m_select_01;
            m_select_01.set_vector(&m_encseq);
            m_select_10 = wt.m_select_10;
            m_select_10.set_vector(&m_encseq);
            m_select_11 = wt.m_select_11;
            m_select_11.set_vector(&m_encseq);
        }

        //! Assignment operator
        wt_dna& operator=(const wt_dna& wt)
        {
            wt_dna tmp(wt);
            tmp.swap(*this);
            return *this;
        }

        //! Swap operator
        void swap(wt_dna& fs)
        {
            if (this != &fs) {
                m_encseq.swap(fs.m_encseq);
                std::swap(m_size, fs.m_size);

                // swap rank support
                util::swap_support(m_rank_00, fs.m_rank_00, &m_encseq, &(fs.m_encseq));
                util::swap_support(m_rank_01, fs.m_rank_01, &m_encseq, &(fs.m_encseq));
                util::swap_support(m_rank_10, fs.m_rank_10, &m_encseq, &(fs.m_encseq));
                util::swap_support(m_rank_11, fs.m_rank_11, &m_encseq, &(fs.m_encseq));

                // swap select support
                util::swap_support(m_select_00, fs.m_select_00, &m_encseq, &(fs.m_encseq));
                util::swap_support(m_select_01, fs.m_select_01, &m_encseq, &(fs.m_encseq));
                util::swap_support(m_select_10, fs.m_select_10, &m_encseq, &(fs.m_encseq));
                util::swap_support(m_select_11, fs.m_select_11, &m_encseq, &(fs.m_encseq));
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
            size_type j = i << 1;
            return (m_encseq[j] << 1) | m_encseq[j + 1];
        }

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{1} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const
        {
            assert(i <= size());

            if (0 == i)  {
                return 0;
            }

            size_type j = (i << 1);

            switch (c) {
                case 0:
                    return m_rank_00.rank(j);
                case 1:
                    return m_rank_01.rank(j);
                case 2:
                    return m_rank_10.rank(j);
                case 3:
                    return m_rank_11.rank(j);
                default:
                    std::cout << "wrong char " << (unsigned int)c << "\n";
                    assert(!"reached");
            }
        }

        //! Calculates how many occurrences of symbol input[i] are in the prefix [0..i-1] of the original input.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(input[i],i), input[i])
         *  \par Time complexity
         *       \f$ \Order{1} + One access to the encoded sequence \f$
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type> inverse_select(size_type i)const
        {
            assert(i < size());

            size_type j = (i << 1);
            uint8_t enc_char = (m_encseq[j] << 1) | m_encseq[j + 1];
            return std::make_pair(rank(i, enc_char), enc_char);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{n} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const
        {
            assert(1 <= i and i <= rank(size(), c));

            switch (c) {
                case 0:
                    return m_select_00.select(i) >> 1;
                case 1:
                    return m_select_01.select(i) >> 1;
                case 2:
                    return m_select_10.select(i) >> 1;
                case 3:
                    return m_select_11.select(i) >> 1;
                default:
                    assert(!"reached");
            }
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
            written_bytes += m_rank_00.serialize(out, child, "rank_00");
            written_bytes += m_rank_01.serialize(out, child, "rank_01");
            written_bytes += m_rank_10.serialize(out, child, "rank_10");
            written_bytes += m_rank_11.serialize(out, child, "rank_11");
            written_bytes += m_select_00.serialize(out, child, "select_00");
            written_bytes += m_select_01.serialize(out, child, "select_01");
            written_bytes += m_select_10.serialize(out, child, "select_10");
            written_bytes += m_select_11.serialize(out, child, "select_11");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            m_encseq.load(in);
            m_rank_00.load(in, &m_encseq);
            m_rank_01.load(in, &m_encseq);
            m_rank_10.load(in, &m_encseq);
            m_rank_11.load(in, &m_encseq);
            m_select_00.load(in, &m_encseq);
            m_select_01.load(in, &m_encseq);
            m_select_10.load(in, &m_encseq);
            m_select_11.load(in, &m_encseq);
        }
};
}

#endif
