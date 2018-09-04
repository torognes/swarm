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

#include "twobit_vector.hpp"

#include <cassert>
#include <iostream>
#include <stdexcept>

// ================================================================================================
//     Typedefs and Constants
// ================================================================================================

const TwobitVector::WordType TwobitVector::all_0_ = 0ul;

const TwobitVector::WordType TwobitVector::all_1_ = (((1ul << 32) - 1) << 32)  +  ((1ul << 32) - 1);

const TwobitVector::WordType TwobitVector::bit_mask_[ TwobitVector::kValuesPerWord ] =
{
    3ul << 0,   3ul << 2,   3ul << 4,   3ul << 6,   3ul << 8,   3ul << 10,  3ul << 12,  3ul << 14,
    3ul << 16,  3ul << 18,  3ul << 20,  3ul << 22,  3ul << 24,  3ul << 26,  3ul << 28,  3ul << 30,
    3ul << 32,  3ul << 34,  3ul << 36,  3ul << 38,  3ul << 40,  3ul << 42,  3ul << 44,  3ul << 46,
    3ul << 48,  3ul << 50,  3ul << 52,  3ul << 54,  3ul << 56,  3ul << 58,  3ul << 60,  3ul << 62
};

const TwobitVector::WordType TwobitVector::ones_mask_[ TwobitVector::kValuesPerWord ] =
{
    TwobitVector::all_0_,       TwobitVector::all_1_ >> 62,
    TwobitVector::all_1_ >> 60, TwobitVector::all_1_ >> 58,
    TwobitVector::all_1_ >> 56, TwobitVector::all_1_ >> 54,
    TwobitVector::all_1_ >> 52, TwobitVector::all_1_ >> 50,
    TwobitVector::all_1_ >> 48, TwobitVector::all_1_ >> 46,
    TwobitVector::all_1_ >> 44, TwobitVector::all_1_ >> 42,
    TwobitVector::all_1_ >> 40, TwobitVector::all_1_ >> 38,
    TwobitVector::all_1_ >> 36, TwobitVector::all_1_ >> 34,
    TwobitVector::all_1_ >> 32, TwobitVector::all_1_ >> 30,
    TwobitVector::all_1_ >> 28, TwobitVector::all_1_ >> 26,
    TwobitVector::all_1_ >> 24, TwobitVector::all_1_ >> 22,
    TwobitVector::all_1_ >> 20, TwobitVector::all_1_ >> 18,
    TwobitVector::all_1_ >> 16, TwobitVector::all_1_ >> 14,
    TwobitVector::all_1_ >> 12, TwobitVector::all_1_ >> 10,
    TwobitVector::all_1_ >> 8,  TwobitVector::all_1_ >> 6,
    TwobitVector::all_1_ >> 4,  TwobitVector::all_1_ >> 2
};

// ================================================================================================
//     Constructors and Rule of Five
// ================================================================================================

TwobitVector::TwobitVector( size_t size )
    : size_( size )
    , data_( ( size / kValuesPerWord ) + ( size % kValuesPerWord == 0 ? 0 : 1 ), 0 )
{}

// ================================================================================================
//     Accessors
// ================================================================================================

size_t TwobitVector::size() const
{
    return size_;
}

size_t TwobitVector::data_size() const
{
    assert( ( size_ / kValuesPerWord ) + ( size_ % kValuesPerWord == 0 ? 0 : 1 ) == data_.size() );
    return data_.size();
}

TwobitVector::ValueType TwobitVector::operator [] ( size_t index ) const
{
    return get( index );
}

TwobitVector::WordType const& TwobitVector::data_at( size_t index ) const
{
    return data_.at( index );
}

TwobitVector::WordType& TwobitVector::data_at( size_t index )
{
    return data_.at( index );
}

TwobitVector::WordType TwobitVector::hash() const
{
#ifdef ZOBRIST
    WordType result = zobrist_hash((unsigned char*)(data_.data()), size());
#else
    WordType result = static_cast< WordType >( size() );
    for( auto s : data_ ) {
        result ^= s;
    }
#endif
    return result;
}

// ================================================================================================
//     Operators
// ================================================================================================

bool TwobitVector::operator == ( TwobitVector const& other ) const
{
    if( size_ != other.size_ ) {
        return false;
    }
    for( size_t i = 0; i < data_.size(); ++i ) {
        if( data_[i] != other.data_[i] ) {
            return false;
        }
    }
    return true;
}

bool TwobitVector::operator != ( TwobitVector const& other ) const
{
    return !(*this == other);
}

bool TwobitVector::validate() const
{
    // Check if the size is correct.
    if( ( size_ / kValuesPerWord ) + ( size_ % kValuesPerWord == 0 ? 0 : 1 ) != data_.size() ) {
        std::cout << "Sizes do not match.\n";
        return false;
    }

    // Check if the zero padding at the end is correct
    // (only if we do have padding though).
    if(
        ( size_ % kValuesPerWord != 0 )
        && (( data_.back() & ~ ones_mask_[ size_ % kValuesPerWord ] ) != 0ul )
    ) {
        std::cout << "Invalid padding bits.\n";
        return false;
    }

    return true;
}

// ================================================================================================
//     Modifiers
// ================================================================================================

void TwobitVector::clear()
{
    size_ = 0;
    data_.clear();
}

// ================================================================================================
//     Insert
// ================================================================================================

void TwobitVector::insert_at( size_t index, TwobitVector::ValueType value )
{
    if( index > size_ ) {
        throw std::out_of_range( "TwobitVector::insert_at: Invalid index." );
    }

    // Shorthands.
    size_t const word_id = index / kValuesPerWord;
    size_t const segm_id = index % kValuesPerWord;

    // If the last word is fully used, we need to add a new one.
    if( size_ % kValuesPerWord == 0 ) {
        data_.push_back( 0ul );
    }

    // Shift all the data right of the insertion word by one.
    for( size_t i = data_.size() - 1; i > word_id; --i ) {

        // Shift the data by one value. We do not loose anything, because the value that is
        // shifted out of the word was already processed in a previous iteration of this
        // loop (or is zero anyway in the first iteration).
        data_[ i ] <<= 2;

        // Get the value of the next word (the one that will be shifted away in the next
        // iteration of this loop), and move it to the other end of the word.
        auto bleed = data_[ i - 1 ] >> ( sizeof( WordType ) * 8 - 2 );

        // Now add this to the current word, so that the bits that are zero (because of the
        // previous shifting) are filled again with the correct value.
        data_[ i ] |= bleed;
    }

    // Get the values in the insertion word that are right of the insertion position.
    auto const remainder = data_[ word_id ] & ( ~ ones_mask_[ segm_id ] );

    // Delete those values in the word.
    data_[ word_id ] &= ones_mask_[ segm_id ];

    // Restore them, shifted by one position. Now we have room for the actual insertion.
    data_[ word_id ] |= ( remainder << 2 );

    // Shift the insertion value to its position, store it in the word, and adjust the size.
    auto const val_shifted = static_cast< WordType >( value ) << ( 2 * segm_id );
    data_[ word_id ] |= val_shifted;
    ++size_;
}

// ================================================================================================
//     Remove
// ================================================================================================

void TwobitVector::remove_at( size_t index )
{
    if( index >= size_ ) {
        throw std::out_of_range( "TwobitVector::remove_at: Invalid index." );
    }

    // Shorthands.
    size_t const word_id = index / kValuesPerWord;
    size_t const segm_id = index % kValuesPerWord;

    // If the position in the word is not the last:
    if( segm_id < kValuesPerWord - 1 ) {

        // Get the part of the word that needs to be shifted.
        auto remainder = data_[ word_id ] & ( ~ ones_mask_[ segm_id + 1 ] );

        // Delete this part.
        data_[ word_id ] &= ones_mask_[ segm_id ];

        // Reset it with the shifted rest values.
        data_[ word_id ] |= ( remainder >> 2 );

    // If it is the last position in the word, we do not need to shift a rest,
    // but can just unset the last value.
    } else {
        data_[ word_id ] &= ones_mask_[ segm_id ];
    }

    // If the word of the deletion is not the last, we need to shift values.
    if( word_id < data_.size() - 1 ) {

        // Get the first value of the next word and store it as the last value in the
        // word where we just deleted a value.
        auto bleed = data_[ word_id + 1 ] << ( sizeof( WordType ) * 8 - 2 );
        data_[ word_id ] |= bleed;

        // Move all values in the remaining words (except the last one) by one.
        for( size_t i = word_id + 1; i < data_.size() - 1; ++i ) {
            bleed = data_[ i + 1 ] << ( sizeof( WordType ) * 8 - 2 );
            data_[ i ] >>= 2;
            data_[ i ] |= bleed;
        }

        // The last word does not need to store the last value of following words, so we
        // can just shift it.
        data_.back() >>= 2;
    }

    // Adjust the size. If we now have a useless word at the end of the vector, remove it.
    --size_;
    if( size_ % kValuesPerWord == 0 ) {
        data_.pop_back();
    }

    // Assert that the size of the vector is correct.
    assert(
        ( size_ / kValuesPerWord ) + ( size_ % kValuesPerWord == 0 ? 0 : 1 ) == data_.size()
    );
}
