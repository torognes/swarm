#ifndef SWARM_TWOBIT_VECTOR_FUNCTIONS_H_
#define SWARM_TWOBIT_VECTOR_FUNCTIONS_H_

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

#include <string>

// =================================================================================================
//     Strings
// =================================================================================================

/**
 * @brief Translate a char into TwobitVector::ValueType.
 *
 * Valid chars are `A`, `C`, `G`, `T`, and their lower case variants.
 */
TwobitVector::ValueType translate_from_nucleic_acid( char site );

/**
 * @brief Translate a TwobitVector::ValueType into its char representation.
 *
 * This gives one of the values `A`, `C`, `G` and `T`.
 */
char translate_to_nucleic_acid( TwobitVector::ValueType value );

/**
 * @brief Turn a string of nucleic acids into a TwobitVector.
 */
TwobitVector from_nucleic_acids( std::string const& sequence );

/**
 * @brief Turn a TwobitVector into its string representation of nucleic acids.
 */
std::string to_nucleic_acids( TwobitVector const& vec );

/**
 * @brief Return a string with a bit-representation of a TwobitVector.
 *
 * It returns the words of the vector with bits in the order of the underlying integer type.
 * This is mainly useful for debugging and testing.
 */
std::string bitstring( TwobitVector const& vec );

/**
 * @brief Return a string with a bit-representation of a TwobitVector::WordType.
 *
 * It returns the word with bits in the order of the underlying integer type.
 * This is mainly useful for debugging and testing.
 */
std::string bitstring( TwobitVector::WordType const& vec );

#endif // include guard
