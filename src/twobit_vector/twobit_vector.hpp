#ifndef SWARM_TWOBIT_VECTOR_H_
#define SWARM_TWOBIT_VECTOR_H_

/*
    SWARM Twobit Vector
    Copyright (C) 2014-2018 Lucas Czech and HITS gGmbH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include <cstddef>
#include <cstdint>
#include <vector>

/**
 * @brief Provides a vector of two-bit values, mainly used for storing nucleotide (ACGT) data.
 */
class TwobitVector
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Constants
    // -------------------------------------------------------------------------

    /**
     * @brief Underlying word type for the bitvector.
     *
     * We use 64bit words to store the 2bit values (of type ValueType), so that we get best
     * speed on modern architectures.
     */
    using WordType = uint64_t;

    /**
     * @brief Value Type enumeration for the elements of a TwobitVector.
     *
     * The values 0-3 are named `A`, `C`, `G` and `T`, respectively, in order to ease the
     * use with nucleotide sequences.
     *
     * The underlying value of the enum is WordType, so that a cast does not need to
     * convert internally.
     */
    enum class ValueType : WordType {
        A = 0,
        C = 1,
        G = 2,
        T = 3
    };

    /**
     * @brief Constant that holds the number of values (of tyoe ValueType) that are stored
     * in a single word in the vector.
     *
     * As we use 64bit words, this value is 32.
     */
    static const size_t kValuesPerWord = sizeof( WordType ) * 8 / 2;

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    /**
     * @brief Default constructor that initializes an empty vector.
     */
    TwobitVector() = default;

    /**
     * @brief Constructor that initializes the vector with `size` many zero values.
     */
    TwobitVector( size_t size );

    ~TwobitVector() = default;

    TwobitVector( TwobitVector const& ) = default;
    TwobitVector( TwobitVector&& )      = default;

    TwobitVector& operator= ( TwobitVector const& ) = default;
    TwobitVector& operator= ( TwobitVector&& )      = default;

    // -------------------------------------------------------------------------
    //     Accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Return the size of the vector, that is, how many values (of type ValueType)
     * it currently holds.
     */
    size_t size() const;

    /**
     * @brief Return the number of words (of type WordType) that are used to store the values
     * in the vector.
     */
    size_t data_size() const;

    /**
     * @brief Get the value at a position in the vector.
     */
    ValueType get( size_t index ) const;

    /**
     * @brief Alias for get().
     */
    ValueType operator [] ( size_t index ) const;

    /**
     * @brief Return a single word of the vector.
     *
     * This is useful for external functions that want to directly work on the underlying bit
     * representation.
     */
    WordType const& data_at( size_t index ) const;

    /**
     * @brief Return a single word of the vector.
     *
     * This is useful for external functions that want to directly work on the underlying bit
     * representation.
     */
    WordType& data_at( size_t index );

    /**
     * @brief Calculate a hash value of the vector, based on its size() and the xor of
     * all its words.
     *
     * This is a simple function, but might just be enough for using it in a hashmap.
     */
    WordType hash() const;

    // -------------------------------------------------------------------------
    //     Operators
    // -------------------------------------------------------------------------

    /**
     * @brief Equality operator.
     */
    bool operator == ( TwobitVector const& other ) const;

    /**
     * @brief Inequality operator, opposite of operator==().
     */
    bool operator != ( TwobitVector const& other ) const;

    /**
     * @brief Validation function that checks some basic invariants.
     *
     * This is mainly useful in testing. The function checks whether the vector is correctly
     * sized and contains zero padding at its end.
     */
    bool validate() const;

    // -------------------------------------------------------------------------
    //     Modifiers
    // -------------------------------------------------------------------------

    /**
     * @brief Clear the vector, so that it contains no data.
     */
    void clear();

    /**
     * @brief Set a value at a position in the vector.
     */
    void set( size_t index, ValueType value );

    /**
     * @brief Insert a value at a position.
     *
     * The size() is increased by one.
     */
    void insert_at( size_t index, ValueType value );

    /**
     * @brief Remove the value at a position.
     *
     * The size() is decreased by one.
     */
    void remove_at( size_t index );

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

private:

    /**
     * @brief Internal constant that holds an all-zero word.
     */
    static const WordType all_0_;

    /**
     * @brief Internal constant that holds an all-one word.
     */
    static const WordType all_1_;

    /**
     * @brief Internal bitmask that has two bits set to one for each value position in the word.
     *
     * The values are
     *
     *     bit_mask_[0]  == 00 00 ... 00 11
     *     bit_mask_[1]  == 00 00 ... 11 00
     *     ...
     *     bit_mask_[31] == 11 00 ... 00 00
     *
     * This is useful for setting or unsetting single values in a word.
     */
    static const WordType bit_mask_[  kValuesPerWord ];

    /**
     * @brief Internal mask that holds as many consecutive all-one values as the position in the
     * array tells.
     *
     * The element at position `i` in this mask contains `i` many all-one values, starting from
     * the right. (An all-one value for two bit values is 11.)
     *
     *     ones_mask_[0]  == 00 00 ... 00 00
     *     ones_mask_[1]  == 00 00 ... 00 11
     *     ones_mask_[2]  == 00 00 ... 11 11
     *     ...
     *     ones_mask_[31] == 00 11 ... 11 11
     *
     * This mask is used for extracting remainders of words (all values left or right of a
     * certain position).
     */
    static const WordType ones_mask_[ kValuesPerWord ];

    size_t                size_ = 0;
    std::vector<WordType> data_;

};

#endif // include guard
