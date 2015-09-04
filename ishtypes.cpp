/* Copyright (c) 2015, UC San Diego CURE Program
  All rights reserved.
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  
  1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  
  2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  
  3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
  This program counts Ish types.
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <chrono>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <sstream>

using Ui64vec = std::vector<uint64_t>;
using Ui64vecvec = std::vector<Ui64vec*>;

// Stores Ish type data for a fixed Ish type
struct TypeData {
  int c;
  uint64_t total;
  Ui64vec primes;
};

// Stores the result of an enumeration of Ish types
struct EnumerationData {
  std::vector<TypeData*> type_data;
  uint64_t total;
  Ui64vec primes;

  void Destroy() {
    for (std::vector<TypeData*>::iterator itr = type_data.begin();
      itr != type_data.end(); ++itr) {
      delete *itr;
    }
    type_data.clear();
    type_data.shrink_to_fit();
  }
};

/*
  ThousandSeparate: adds thousand separators to a number

  Parameters:
    number    - The number to be separated

  Output:     The thousand separated number as a string
*/
std::string ThousandSeparate(uint64_t number) {
  std::ostringstream ss;
  std::string ret_val = "";
  int pos = 0;
  ss << number;
  std::string str = ss.str();
  ss.clear();

  for (std::string::iterator itr = str.end()-1; itr >= str.begin(); --itr)
  {
    pos++;

    ret_val = (*itr) + ret_val;

    if (pos % 3 == 0 && itr != str.begin())
      ret_val = "," + ret_val;
  }

  return ret_val;
}

/*
  PrintVec: Prints an integer vector on screen

  Parameters:
    vec   - The vector to print
*/
void PrintVec(Ui64vec vec) {
  for (Ui64vec::iterator itr = vec.begin(); itr != vec.end(); ++itr) {
    std::cout << ThousandSeparate(*itr);

    if (itr != vec.end() - 1)
      std::cout << " * ";
  }
}

/*
  RecurseSeqs: Recursive function to count all strictly increasing integer
          sequences bounded below by 1 and bounded above by an
          integer sequence

  Parameters:
    total       - reference to an integer holding the total number of
            sequences found
    upper_seq   - the sequence bounding from above
    pos      - current position in the sequence
    a_prev      - the value of the position prior to the current position
*/
void RecurseSeqs(uint64_t *total, Ui64vec* upper_seq, uint64_t pos,
  uint64_t a_prev) {
  // Count sequence when we have reached the end
  if (pos == upper_seq->size()) {
    (*total)++;
    return;
  }

  // Get lower and upper bounds for the next step
  uint64_t lower = a_prev + 1;
  uint64_t upper = *min_element(begin((*upper_seq)) + pos, end((*upper_seq)));

  // Go through all the possible next steps
  for (uint64_t i=lower; i <= upper; ++i) {
    RecurseSeqs(total, upper_seq, pos+1, i);
  }
}

/*
  CountSequences: Helper function to find the total number of strictly
            increasing integer sequences bounded below by 1 and
            above by another sequence

  Parameters:
    upper_seq      - The upper-bound sequence

  Output:        The number of sequences

*/
uint64_t CountSequences(Ui64vec* upper_seq) {
  uint64_t total = 0;

  RecurseSeqs(&total, upper_seq, 0, 0);

  return total;
}

/*
  GetUpperSeq: Constructs an upper-bound sequence from the permutation from
          which the Ish types are to be constructed

  Parameters:
    type_seq    - The permutation
    m        - the value of m in the corresponding Shi^m(n)
    n        - the value of n in the corresponding Shi^m(n)

  Output:     The upper-bound sequence
*/
Ui64vec* GetUpperSeq(Ui64vec *type_seq, uint64_t m, uint64_t n) {
  Ui64vec* ret_vec = new Ui64vec();

  for (Ui64vec::iterator itr = type_seq->begin(); itr != type_seq->end();
    ++itr) {
    // The upper bound is given by (m-1)n+i-1, where i is the number in
    //   the permutation
    uint64_t bound = (m-1) * n + *itr - 1;
    ret_vec->push_back(bound);
  }

  return ret_vec;
}

/*
  FindPermutations: Recursively finds all possible permutations of length k
            chosen from a set N

  Parameters:
    N        - Set to pick from
    k        - Length of permutation

  Output:     A list of permutations

  Algorithm:

    f(N,k)
    var list

    for each n in N
      if k=1
        add n to list
      else
        N' = N - {n}
        S = f(N', k-1)
        for each s in S
          add {S, n} to list

    ret list
*/
Ui64vecvec* FindPermutations(Ui64vec* N, uint64_t k) {
  Ui64vecvec* list = new Ui64vecvec();

  for (Ui64vec::iterator cur_N = N->begin(); cur_N != N->end(); ++cur_N) {
    if (k == 1) {
      list->push_back(new Ui64vec());
      list->back()->push_back(*cur_N);
    } else {
      // Create a new N with the same contents as N but without the
      //   current element of N
      Ui64vec* new_N = new Ui64vec();

      for (Ui64vec::iterator copy_itr = N->begin(); copy_itr != cur_N;
        ++copy_itr) {
        new_N->push_back(*copy_itr);
      }

      for (Ui64vec::iterator copy_itr = cur_N+1; copy_itr != N->end();
        ++copy_itr) {
        new_N->push_back(*copy_itr);
      }

      // Recursively find permutations on the smaller set
      Ui64vecvec* new_list = FindPermutations(new_N, k-1);

      // Add all the permutations found, with cur_N appended to each one
      for (Ui64vecvec::iterator rec_itr = new_list->begin();
        rec_itr != new_list->end(); ++rec_itr) {
        list->push_back(new Ui64vec());

        for (Ui64vec::iterator rec_inner_itr = (*rec_itr)->begin();
          rec_inner_itr != (*rec_itr)->end(); ++rec_inner_itr) {
          list->back()->push_back(*rec_inner_itr);
        }

        list->back()->push_back(*cur_N);
      }

      // Free memory
      new_list->clear();
      new_list->shrink_to_fit();
    }
  }

  return list;
}

/*
  PrimesBelow: Finds all prime numbers below an upper bound using the sieve
          of Eratosthenes

  Parameters:
    upper_bound    - The upper bound on the primes
*/
Ui64vec PrimesBelow(uint64_t upper_bound) {
  assert(upper_bound >= 2);

  Ui64vec prime_candidates(upper_bound, 0);
  Ui64vec primes;

  for (uint64_t i=2; i < floor(sqrt(upper_bound)); ++i) {
    for (uint64_t j=i*i; j < upper_bound; j+=i) {
      prime_candidates[j-1] = 1;
    }
  }

  for (uint64_t i=2; i < upper_bound; ++i) {
    if (prime_candidates[i - 1] == 0)
      primes.push_back(i);
  }

  return primes;
}

/*
  Factor: Gives the prime factorization of a number

  Parameters:
    number      - The number to factor
    primes      - A list of possible prime factors

  Output:        The list prime composition of the number
*/
Ui64vec Factor(uint64_t number, Ui64vec* primes) {
  Ui64vec factors;

  if (number < 1)
    number = -number;

  if (number == 0) {
    factors.push_back(0);
    return factors;
  }

  for (Ui64vec::iterator prime = primes->begin();
    prime != primes->end() && number != 1; ++prime) {
    while (number % (*prime) == 0) {
      number /= (*prime);
      factors.push_back((*prime));
    }
  }

  return factors;
}

/*
  IsInt: Checks if a string is an integer

  Parameters:
    str      - The string to check

  Output:     1 if str is an integer, 0 otherwise
*/
uint64_t IsInt(char* str) {
  uint64_t retval = 1;

  for (uint64_t i=0; i < strlen(str) && retval == 1; ++i) {
    if (!isdigit(str[i]))
      retval = 0;
  }

  return retval;
}

/*
  Enumerate: Gets the total number of Ish types for a given m and n together
              with additional statistics

  Parameters:
    m        - the value of m in the corresponding Shi^m(n)
    n        - the value of n in the corresponding Shi^m(n)

  Output:    enumeration data
*/
EnumerationData Enumerate(int m, int n)
{
  uint64_t total_sum = 0;
  EnumerationData ret_struct;
  Ui64vec primes;

  // N=[n]-{1}
  Ui64vec* N = new Ui64vec();
  for (uint64_t i=2; i <= n; ++i) {
    N->push_back(i);
  }

  // Loop over Ish types of a fixed length k
  for (uint64_t k=1; k < n; ++k) {
    Ui64vecvec* lower_seqs = FindPermutations(N, k);

    uint64_t running_total = 0;

    // Go through all the lower permutations
    for (Ui64vecvec::iterator lower_seq = lower_seqs->begin();
      lower_seq != lower_seqs->end(); ++lower_seq) {
      Ui64vec* upper_seq = GetUpperSeq((*lower_seq), m, n);

      running_total += CountSequences(upper_seq);

      // Free memory
      upper_seq->clear();
      upper_seq->shrink_to_fit();

      delete *lower_seq;
    }

    // Free memory
    lower_seqs->clear();
    lower_seqs->shrink_to_fit();

    ret_struct.type_data.push_back(new TypeData);
    ret_struct.type_data.back()->c = k;
    ret_struct.type_data.back()->total = running_total;

    total_sum += running_total;
  }

  delete N;

  ret_struct.total = total_sum;

  // Prime factorization

  primes = PrimesBelow(total_sum + 1);

  ret_struct.primes = Factor(ret_struct.total, &primes);

  for (std::vector<TypeData*>::iterator itr = ret_struct.type_data.begin();
    itr != ret_struct.type_data.end(); ++itr) {
    (*itr)->primes = Factor((*itr)->total, &primes);
  }

  primes.clear();
  primes.shrink_to_fit();

  return ret_struct;
}

/*
  PrintStatistics : prints the total number of Ish types and the number of
                    types of fixed length, together with prime factorizations

  Parameters:
    enum_data      - Data from Enumerate(m, n)
*/
void PrintStatistics(EnumerationData enum_data) {
  std::cout << "Total Ish types:     "
    << ThousandSeparate(enum_data.total);
  std::cout << "\nPrime Factorization: ";
  PrintVec(enum_data.primes);

  std::cout << std::endl << std::endl;

  std::cout << " c | " << std::setw(12) << "Total" << " | Factorization"
    << std::endl;

  for (std::vector<TypeData*>::iterator itr = enum_data.type_data.begin();
    itr != enum_data.type_data.end(); ++itr) {
    std::cout << " " << (*itr)->c << " | " << std::setw(12)
      << ThousandSeparate((*itr)->total) << " | ";
    PrintVec((*itr)->primes);
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

/*
  Entry point
*/
int main(int argc, char **argv) {
  uint64_t *total;
  uint64_t total_sum = 0;
  uint64_t m = 0, n = 0;
  EnumerationData enum_data;

  std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
  std::chrono::duration<double> total_time;

  if (argc != 3) {
    std::cout << "Usage: ishtypes m n" << std::endl;
    return EXIT_FAILURE;
  } else {
    if (!IsInt(argv[1]) || !IsInt(argv[2])) {
      std::cout << "Usage: ishtypes m n" << std::endl;
      return EXIT_FAILURE;
    }

    m = atoi(argv[1]);
    n = atoi(argv[2]);
  }

  start_time = std::chrono::system_clock::now();

  enum_data = Enumerate(m, n);
  PrintStatistics(enum_data);
  enum_data.Destroy();

  end_time = std::chrono::system_clock::now();
  total_time = end_time - start_time;

  std::cout << "Time elapsed: " << total_time.count() << "s" << std::endl;

  return EXIT_SUCCESS;
}
