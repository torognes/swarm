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
#include "iterator_deletions.hpp"
#include "iterator_insertions.hpp"
#include "iterator_substitutions.hpp"
#include "twobit_vector.hpp"

#include <iostream>
#include <string>

// =================================================================================================
//     Examples
// =================================================================================================

void example_basics( std::string const& nucleotides )
{
    // Given a string of nucleotides, turn it into a vector.
    auto vec = from_nucleic_acids( nucleotides );

    // Print the value at a certain position.
    if( vec.size() > 0 ) {
        std::cout << translate_to_nucleic_acid( vec[0] ) << "\n";
    }

    // Set a value at a position in the vector.
    if( vec.size() > 0 ) {
        vec.set( 0, TwobitVector::ValueType::G );
    }

    // Make another vector of a given size, initialized to 0 = A.
    auto vec2 = TwobitVector( 5 );

    // Compare for equality.
    if( vec == vec2 ) {
        std::cout << "Equal!\n";
    }

    // Print the hash.
    std::cout << vec.hash() << "\n";

    // Other member functions like insert_at() or remove_at() work accordingly.
}

void example_output( TwobitVector const& vec )
{
    // Output the vector as a nucleotide sequence.
    std::cout << to_nucleic_acids( vec ) << "\n";

    // Output its two-bit representation, mainly useful for debugging.
    std::cout << bitstring( vec ) << "\n";
}

void example_iterators( TwobitVector const& vec )
{
    // Example of how to use the iterators for insertions/deletions/substitions.
    // They all work the same, so we here just use deletions as an example.

    // Simple loop that yields all possible deletions of nucleotides.
    // As the range-based for loop automatically dereferences the iterator,
    // this version does not use the nice fast update of the hash offered by the iterator,
    // but instead re-calcualtes the hash from the vector.
    for( auto const& del_vec : iterate_deletions( vec ) ) {
        std::cout << to_nucleic_acids( del_vec ) << ": " << vec.hash() << "\n";
    }

    // Use the iterator explicitly. This way, we can access the fast update of the hash.
    // This should be used in actual production code.
    auto it_b = IteratorDeletions( vec );
    auto it_e = IteratorDeletions();
    while( it_b != it_e ) {

        // Get the current vector (with one nucleotide deleted), and its hash.
        // Instead of .vector(), one could also simply dereference the iterator.
        std::cout << to_nucleic_acids( it_b.vector() ) << ": " << it_b.hash() << "\n";

        // Iterate to next deletion
        ++it_b;
    }
}
