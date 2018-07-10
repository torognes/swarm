#ifndef SWARM_TWOBIT_VECTOR_ITERATOR_SUBSTITUTIONS_H_
#define SWARM_TWOBIT_VECTOR_ITERATOR_SUBSTITUTIONS_H_

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

#include "range.hpp"
#include "twobit_vector.hpp"
#include "functions.hpp"

#include <cassert>
#include <iterator>

// =================================================================================================
//     Iterator Substitutions
// =================================================================================================

/**
 * @brief Take a TwobitVector sequence and iterate over all possible substitutions of single nucleotides.
 */
class IteratorSubstitutions
{

public:

    // -----------------------------------------------------
    //     Typedefs
    // -----------------------------------------------------

    using iterator_category = std::forward_iterator_tag;
    using self_type         = IteratorSubstitutions;
    using value_type        = TwobitVector;

    // -----------------------------------------------------
    //     Constructors and Rule of Five
    // -----------------------------------------------------

    IteratorSubstitutions ()
        : origin_( nullptr )
        , vec_()
        , pos_( 0 )
        , cnt_( 0 )
        , hash_( 0 )
    {}

    explicit IteratorSubstitutions( TwobitVector const& vector )
        : origin_( &vector )
        , vec_( vector )
        , pos_( 0 )
        , cnt_( 0 )
    {
        // Iterate to the first substitution, and do a first hash calculation.
        // Later iterations will just update all of this.
        operator++();
        hash_ = vec_.hash();
    }

    ~IteratorSubstitutions() = default;

    IteratorSubstitutions( IteratorSubstitutions const& ) = default;
    IteratorSubstitutions( IteratorSubstitutions&& )      = default;

    IteratorSubstitutions& operator= ( IteratorSubstitutions const& ) = default;
    IteratorSubstitutions& operator= ( IteratorSubstitutions&& )      = default;

    // -----------------------------------------------------
    //     Operators
    // -----------------------------------------------------

    value_type const& operator * ()
    {
        return vec_;
    }

    value_type const* operator -> ()
    {
        return &vec_;
    }

    self_type& operator ++ ()
    {
        // We use four xor's at the current position to cycle through the variants:
        // The first thee are the substitutions, the last one then restores the original value.
        // For this, we use the xor order 01 11 01 11.
        //
        // The table shows that this works for all four possible values.
        //
        //            | 00 01 10 11
        //     ---------------------
        //     0 | 01 | 01 00 11 10
        //     1 | 11 | 10 11 00 01
        //     2 | 01 | 11 10 01 00
        //     3 | 11 | 00 01 10 11

        // Helper function that cycles the value at a position.
        auto cycle = [&] ( size_t pos, size_t& cnt ) {

            // Shorthand.
            size_t shift = ( 2 * ( pos % TwobitVector::kValuesPerWord ));

            // Move the needed xor value to the position in the word and apply it.
            TwobitVector::WordType xor_val = ( cnt % 2 == 0 ? 0x1 : 0x3 );
            vec_.data_at( pos / TwobitVector::kValuesPerWord ) ^= ( xor_val << shift );

            // Update the hash: Remove the current value, store the new one.
            // (We can simply reuse the xor value, as a^b=c <=> b^c=a)
            hash_ ^= ( xor_val << shift );

            ++cnt;
        };

        // Do at least one cycle at the current position.
        cycle( pos_, cnt_ );

        // If we used all three possible substitution values at the current position:
        if( cnt_ == 4 ) {

            // If this is not the last possible position, move to the next one.
            if( pos_ < vec_.size() - 1 ) {
                // Restore original value.
                // cycle( pos_, cnt_ );

                // Move to next position and do a cycle.
                ++pos_;
                cnt_ = 0;
                cycle( pos_, cnt_ );

            // If we are done, set everything to zero, so that the iterator
            // is equal to the default-constructed end-iterator.
            } else {
                origin_ = nullptr;
                vec_.clear();
                pos_  = 0;
                cnt_  = 0;
                hash_ = 0;
            }
        }

        return *this;
    }

    self_type operator ++ (int)
    {
        self_type tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator == ( self_type const& other ) const
    {
        return ( origin_ == other.origin_ ) && ( pos_ == other.pos_ ) && ( cnt_ == other.cnt_ );
    }

    bool operator != ( self_type const& other ) const
    {
        return !( other == *this );
    }

    // -----------------------------------------------------
    //     Members
    // -----------------------------------------------------

    /**
     * @brief Get the position that is currently deleted.
     */
    size_t position() const
    {
        return pos_;
    }

    /**
     * @brief Get the hash value of the current vector.
     */
    TwobitVector::WordType hash() const
    {
        return hash_;
    }

    /**
     * @brief Get the current vector.
     */
    TwobitVector const& vector() const
    {
        return vec_;
    }

private:

    // Store the location of the original vector. This is mainly used for quickly checking
    // whether two iterators refer to the same underlying vector.
    // (We do not want to do a full vector equality check at each iteration.)
    TwobitVector const* origin_;

    // The current vector, which always has an additional value (compared to the original vector).
    TwobitVector vec_;

    // The position where currently a value is inserted.
    size_t pos_;

    // A counter for the possible substitutions possibilities.
    size_t cnt_;

    // The hash value of the current vector.
    TwobitVector::WordType hash_;

};

// =================================================================================================
//     Range Wrapper
// =================================================================================================

/**
 * @brief Helper function that allows range-based for loops.
 */
inline Range< IteratorSubstitutions > iterate_substitutions( TwobitVector const& vector )
{
    return {
        IteratorSubstitutions( vector ),
        IteratorSubstitutions()
    };
}

#endif // include guard
