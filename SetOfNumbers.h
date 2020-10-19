#ifndef SETOFNUMBERS_H
#define SETOFNUMBERS_H

/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of my employers nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL CORNELIS JAN VAN DE STREEK BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************* */

#include <cstddef> // For definition of size_t
#include <vector>

/*
 * This class represents a set of numbers.
 *
 * Numbers are not necessarily unique, this is configurable.
 *
 *  @@ Add permutations, and "k-in-N"
 *
 */

class SetOfNumbers
{
public:
    
    
    // ALLOWED: duplicates are allowed
    // AUTO_REMOVE: duplicates are not allowed and are silently removed
    // THROW: Throw if a duplicate is encountered
    enum DuplicatesPolicy { ALLOWED, AUTO_REMOVE, THROW };

    // Default constructor, creates an empty set,
    explicit SetOfNumbers( const DuplicatesPolicy duplicates_policy = ALLOWED );

    // In keeping with C++ convention: zero-based.
    // Fills the set with the numbers 0, 1, ..., nvalues-1.
    explicit SetOfNumbers( const size_t nvalues, const DuplicatesPolicy duplicates_policy = ALLOWED );

    // Fills the set with the numbers [begin,end]
    SetOfNumbers( const size_t begin, const size_t end, const DuplicatesPolicy duplicates_policy = ALLOWED );

    explicit SetOfNumbers( const std::vector< size_t > & values, const DuplicatesPolicy duplicates_policy = ALLOWED );

    // Throws if set to THROW and set contains duplicates
    // If set to AUTO_REMOVE, duplicates are removed
    void set_duplicates_policy( const DuplicatesPolicy duplicates_policy );

    DuplicatesPolicy duplicates_policy() const { return duplicates_policy_; }

    // Throws if set to true and set currently empty
    void set_empty_is_allowed( const bool empty_is_allowed );

    bool empty_is_allowed() const { return empty_is_allowed_; }

    bool empty() const { return values_.empty(); }
    
    // Returns the number of numbers in the set.
    size_t size() const { return values_.size(); }

    void reserve( const size_t desired_size ) { values_.reserve( desired_size ); }

    size_t value( const size_t index ) const;

    std::vector< size_t > values() const;

    // DuplicatesPolicy is set to AUTO_REMOVE
    SetOfNumbers unique_values() const;

    // Returns how often value occurs in the set.
    size_t frequency( const size_t value ) const;

    // If value currently already in the set, adds it again
    void add( const size_t value );

    // Does nothing if value currently not in the set.
    // If multiple occurrences present, only removes one.
    void remove( const size_t value );

    bool contains( const size_t value ) const;

    bool contains_duplicates() const;

    // In a normal set, would be called union
    // The attributes duplicates_allowed_ and empty_is_allowed_ of the
    // argument values are ignored, the values of *this are kept.
    void add( const SetOfNumbers & set_of_numbers );

    // The attributes duplicates_allowed_ and empty_is_allowed_ of the
    // argument values are ignored, the values of *this are kept.
    void remove( const SetOfNumbers & set_of_numbers );

    // In a normal set, would be called intersection
    // The attribute duplicates_allowed_ of the
    // argument are ignored, the values of *this are kept.
    // The attribute empty_is_allowed_ is true
    // If *this contains three times the same value and the argument contains two times that same value, they have two in common.
    // If the argument contains three times the same value and *this contains two times that same value, they have two in common.
    SetOfNumbers in_common( const SetOfNumbers & set_of_numbers );

    bool operator==( const SetOfNumbers & rhs ) const;
    bool operator!=( const SetOfNumbers & rhs ) const;

    void show() const;

private:
    std::vector< size_t > values_;
    mutable std::vector< size_t > sorted_map_;
    DuplicatesPolicy duplicates_policy_;
    bool empty_is_allowed_;
    mutable bool is_sorted_;
    
    void remove_position( const size_t index );
    void sort() const;
    void check_for_duplicates();
    void check_if_empty() const;
};

// The attributes duplicates_allowed_ and empty_is_allowed_ of the
// arguments are ignored
// duplicates_allowed_ is set to true, empty_is_allowed_ is set to true
SetOfNumbers merge( const SetOfNumbers & lhs, const SetOfNumbers & rhs );

#endif // SETOFNUMBERS_H

