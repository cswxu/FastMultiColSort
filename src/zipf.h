#ifndef     ZIPF_H
#define     ZIPF_H

#include    <cmath>
#include    <cassert>
#include    <iterator>
#include    <algorithm>

namespace multiAttrSort {
    
/**
 * From Byteslice project by Ziqiang Feng
 * 
 * Populate the range with Zipfian random data.
 * 
 * The algorithm starts from the highest ranked (i.e., most frequent) element,
 * and generates consecutive occurrences of the element.
 * The number of occurrences are calculated according to its expected frequency and rounded up.
 * It keeps generating until the end of the range 
 * and shuffles the whole sequence afterwards.
 * Least frequent elements may never be generated if the range is not large enough.
 * \c zipf=0.0 is supposed to generate uniform distributed data.
 * But again, items returned lately by \c next() may not appear.
 * 
 * @param first Beginning of the range
 * @param last  Past-the-end of the range
 * @param N     Cardinality (number of distinct values).
 * @param zipf  The zipf factor. Typically > 1.
 * @param next  A function object return ranked elements from most frequent to less.
 *              On each call to \c next(), it should return the next frequent element.
 */
//@{
template <class RandomAccessIterator, class Next>
void fill_zipf(RandomAccessIterator first, RandomAccessIterator last, 
               uint64_t N, double zipf,
               Next next){
    assert(zipf >= 0);
    // Calculate the denominator
    double denom = 0.0;
    for(uint64_t rank = 1; rank <= N; rank++){
        denom += 1.0 / std::pow(double(rank), zipf);
    }
    
    // Populate the range
    RandomAccessIterator it = first;
    for(uint64_t rank=1; rank <= N; rank++){
        auto elem = next();
        double prob = (1.0 / std::pow(double(rank), zipf)) / denom;
        uint64_t num_occurrences = std::max(1.0, std::floor(std::distance(first, last) * prob));
        it = std::fill_n(it,
                         std::min(num_occurrences, uint64_t(last - it)),
                         elem);
        if(last == it)
            break;
    }
    if(last != it){
        fill_zipf(it, last, N, zipf, next);
    }
    
}

template <class RandomAccessIterator, class Next>
void fill_zipf_random(RandomAccessIterator first, RandomAccessIterator last, 
                      uint64_t N, double zipf,
                      Next next){
    fill_zipf(first, last, N, zipf, next);
    //Shuffle it
    std::random_shuffle(first, last);
}
//@}

}   // namespace

#endif
