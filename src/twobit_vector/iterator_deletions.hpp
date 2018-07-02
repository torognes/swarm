#ifndef SWARM_TWOBIT_VECTOR_ITERATOR_DELETIONS_H_
#define SWARM_TWOBIT_VECTOR_ITERATOR_DELETIONS_H_

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
//     Iterator Deletions
// =================================================================================================

/**
 * @brief Take a TwobitVector sequence and iterate over all possible deletions of single nucleotides.
 */
class IteratorDeletions
{
public:

    // -----------------------------------------------------
    //     Typedefs
    // -----------------------------------------------------

    using iterator_category = std::forward_iterator_tag;
    using self_type         = IteratorDeletions;
    using value_type        = TwobitVector;

    // -----------------------------------------------------
    //     Constructors and Rule of Five
    // -----------------------------------------------------

    IteratorDeletions ()
        : origin_( nullptr )
        , vec_()
        , pos_( 0 )
        , cur_( TwobitVector::ValueType::A )
        , hash_( 0 )
    {}

    explicit IteratorDeletions( TwobitVector const& vector )
        : origin_( &vector )
        , vec_( vector )
        , pos_( 0 )
        , cur_( vector[ 0 ] )
    {
        // Remove the first value (we stored it in the initializer),
        // and do a first hash calculation. Later iterations will just update
        // all of this.
        vec_.remove_at( 0 );
        hash_ = vec_.hash();
    }

    ~IteratorDeletions() = default;

    IteratorDeletions( IteratorDeletions const& ) = default;
    IteratorDeletions( IteratorDeletions&& )      = default;

    IteratorDeletions& operator= ( IteratorDeletions const& ) = default;
    IteratorDeletions& operator= ( IteratorDeletions&& )      = default;

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

    self_type operator ++ ()
    {
        // Toy example:
        // Original: ACGT
        // CGT --> pos = 0, cur = A
        // AGT --> pos = 1, cur = C
        // ACT --> pos = 2, cur = G
        // ACG --> pos = 3, cur = T

        // If we are not done yet:
        if( pos_ < vec_.size() ) {
            // We will swap the value at the current position. So first, get it.
            auto tmp = vec_.get( pos_ );

            // Update the hash: Remove the current value, and store the new one.
            hash_ ^= ((
                static_cast< TwobitVector::WordType >( tmp ) ^
                static_cast< TwobitVector::WordType >( cur_ )
                ) << ( 2 * ( pos_ % TwobitVector::kValuesPerWord ))
            );

            // Now do the swap and move to the next position, so that next time,
            // we will swap the next value.
            vec_.set( pos_, cur_ );
            cur_ = tmp;
            ++pos_;

        // If we are done, set everything to zero, so that the iterator
        // is equal to the default-constructed end-iterator.
        } else {
            origin_ = nullptr;
            vec_.clear();
            pos_  = 0;
            cur_  = TwobitVector::ValueType::A;
            hash_ = 0;
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
        return ( origin_ == other.origin_ ) && ( pos_ == other.pos_ );
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

    // The current vector, which always has a missing value (compared to the original vector).
    TwobitVector vec_;

    // The position that is currently deleted.
    size_t                  pos_;

    // The value that is currently deleted (need it to restore it later).
    TwobitVector::ValueType cur_;

    // The hash value of the current vector.
    TwobitVector::WordType  hash_;

};

// =================================================================================================
//     Range Wrapper
// =================================================================================================

/**
 * @brief Helper function that allows range-based for loops.
 */
inline Range< IteratorDeletions > iterate_deletions( TwobitVector const& vector )
{
    return {
        IteratorDeletions( vector ),
        IteratorDeletions()
    };
}

#endif // include guard
