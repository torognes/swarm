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

#include "functions.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

// =================================================================================================
//     Strings
// =================================================================================================

TwobitVector::ValueType translate_from_nucleic_acid( char site )
{
    switch( site ) {
        case 'a':
        case 'A':
            return TwobitVector::ValueType::A;

        case 'c':
        case 'C':
            return TwobitVector::ValueType::C;

        case 'g':
        case 'G':
            return TwobitVector::ValueType::G;

        case 't':
        case 'T':
            return TwobitVector::ValueType::T;

        default:
            throw std::runtime_error( "Invalid nucleic acid." );
    }
}

char translate_to_nucleic_acid( TwobitVector::ValueType value )
{
    switch( value ) {
        case TwobitVector::ValueType::A: return 'A';
        case TwobitVector::ValueType::C: return 'C';
        case TwobitVector::ValueType::G: return 'G';
        case TwobitVector::ValueType::T: return 'T';
        default:
            throw std::runtime_error( "Invalid twobit value." );
    }
}

TwobitVector from_nucleic_acids( std::string const& sequence )
{
    // We set each value individually.
    // Shifting them should be faster - future optimization!
    auto result = TwobitVector( sequence.size() );
    for( size_t i = 0; i < sequence.size(); ++i ) {
        result.set( i, translate_from_nucleic_acid( sequence[i] ));
    }
    return result;
}

std::string to_nucleic_acids( TwobitVector const& vec )
{
    std::string result;
    result.reserve( vec.size() );

    for( size_t i = 0; i < vec.size(); ++i ) {
        result += translate_to_nucleic_acid( vec[i] );
    }
    return result;
}

std::string bitstring( TwobitVector const& vec )
{
    // Put each word on a new line.
    std::string res = "";
    for( size_t i = 0; i < vec.data_size(); ++i ) {
        res += bitstring( vec.data_at( i )) + "\n";
    }
    return res;
}

std::string bitstring( TwobitVector::WordType const& vec )
{
    // This is an ugly quick hack function. Would be nicer to use bitmasks, but they are hidden
    // inside the TwobitVector class and exposing them just for this purpuse is also not nice.

    // Make a copy, so that we can shift away the processed parts.
    auto cpy = vec;
    std::string res;
    res.reserve( 96 );

    // Go through the word two bits at a time, and store the right-most bits
    // in the string. This way, we obtain a reverse order string, so we need to reverse
    // it later again. Also, note that the bit strings for 1 and 2 are reverses in the switch
    // statement because of this.
    for( size_t i = 0; i < TwobitVector::kValuesPerWord; ++i ) {
        auto tmp = cpy & 0x3;
        switch( tmp ) {
            case 0x0:
                res += "00";
                break;

            case 0x1:
                res += "10";
                break;

            case 0x2:
                res += "01";
                break;

            case 0x3:
                res += "11";
                break;

            default:
                assert( false );
        }
        if( i < TwobitVector::kValuesPerWord - 1 ) {
            res += " ";
        }
        cpy >>= 2;
    }

    std::reverse( res.begin(), res.end() );
    return res;
}
