#ifndef SWARM_TWOBIT_VECTOR_RANGE_H_
#define SWARM_TWOBIT_VECTOR_RANGE_H_

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

// =================================================================================================
//     Range
// =================================================================================================

/**
 * @brief Simple helper wrapper class that allows to store a begin and and end iterator.
 * Mainly intended to be used for custom type range-based for loops.
 */
template <typename IteratorType>
class Range
{
public:

    // -------------------------------------------------------------------------
    //     Member Types
    // -------------------------------------------------------------------------

    using iterator = IteratorType;

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    Range () = default;

    template <class Container>
    Range ( Container& cont )
        : begin_( cont.begin() )
        , end_(   cont.end() )
    {}

    template <class Container>
    Range ( Container const& cont )
        : begin_( cont.begin() )
        , end_(   cont.end() )
    {}

    Range( iterator begin, iterator end )
        : begin_(begin)
        , end_(end)
    {}

    Range( Range const& ) = default;
    Range( Range&& )      = default;

    Range& operator= ( Range const& ) = default;
    Range& operator= ( Range&& )      = default;

    ~Range() = default;

    // -------------------------------------------------------------------------
    //     Iterators
    // -------------------------------------------------------------------------

    iterator begin()
    {
        return begin_;
    }

    iterator end()
    {
        return end_;
    }

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

private:

    iterator begin_;
    iterator end_;

};

#endif // include guard
