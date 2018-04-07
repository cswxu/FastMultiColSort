/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#ifndef COLUMN_H
#define COLUMN_H

#include <cstdint>
#include <iostream>

namespace multiAttrSort{

template <class T>
class SpecialiedColumn;

class Column {
public:
	virtual ~Column() {}
	virtual void* GetColumn() const = 0;
	virtual uint32_t GetWidth() const = 0;
	virtual void SetWidth(uint32_t w) = 0;
	virtual uint64_t GetCardinality() const = 0;
	virtual void SetCardinality(uint64_t card) = 0;

	virtual void SetValueAt(uint64_t val, uint64_t rowIdx) = 0;
	virtual void ShiftValueLeft(uint32_t nbits, uint64_t rowIdx) = 0;
	virtual void LogicalOr(uint64_t orValue, uint64_t rowIdx) = 0;
	virtual uint64_t GetValueAt(uint64_t rowIdx) const = 0;
	virtual uint64_t GetValPortion(uint64_t rowIdx, uint32_t start, uint32_t end) const = 0;
	virtual void ShiftLeftBatch(uint64_t num_rows, uint32_t nbits) = 0;
	virtual void AssembleFromColumn(SpecialiedColumn<uint8_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) = 0;
	virtual void AssembleFromColumn(SpecialiedColumn<uint16_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) = 0;
	virtual void AssembleFromColumn(SpecialiedColumn<uint32_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) = 0;
	virtual void AssembleFromColumn(SpecialiedColumn<uint64_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) = 0;
};

template <class T>
class SpecialiedColumn: public Column {
public:
	virtual ~SpecialiedColumn() {

	}
	void*	GetColumn() const override;
	void	SetColumn(T* values) {values_ = values;}
	uint32_t GetWidth() const override;
	void SetWidth(uint32_t w) override;
	uint64_t GetCardinality() const override;
	void SetCardinality(uint64_t card) override;

	void SetValueAt(uint64_t val, uint64_t rowIdx) override;
	void ShiftValueLeft(uint32_t nbits, uint64_t rowIdx) override;
	void LogicalOr(uint64_t orValue, uint64_t rowIdx) override;
	uint64_t GetValueAt(uint64_t rowIdx) const override;
	uint64_t GetValPortion(uint64_t rowIdx, uint32_t start, uint32_t end) const override;
	void ShiftLeftBatch(uint64_t num_rows, uint32_t nbits) override;

	void AssembleFromColumn(SpecialiedColumn<uint8_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) override;
	void AssembleFromColumn(SpecialiedColumn<uint16_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) override;
	void AssembleFromColumn(SpecialiedColumn<uint32_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) override;
	void AssembleFromColumn(SpecialiedColumn<uint64_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) override;

private:
	T *values_;
	uint32_t width_;	//# of bytes for T
	uint32_t cardinality_;	//# of distinct values
	const uint32_t kPrefetchDistance = 2 * CACHE_LINE_SIZE;
};

template <class T>
inline void* SpecialiedColumn<T>::GetColumn() const {
	return (void *)values_;
}

template <class T>
inline uint32_t SpecialiedColumn<T>::GetWidth() const {
	return width_;
}

template <class T>
inline void SpecialiedColumn<T>::SetWidth(uint32_t w) {
	width_ = w;
}

template <class T>
inline uint64_t SpecialiedColumn<T>::GetCardinality() const {
	return cardinality_;
}

template <class T>
inline void SpecialiedColumn<T>::SetCardinality(uint64_t card) {
	cardinality_ = card;
}


template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::SetValueAt(uint64_t val, uint64_t rowIdx) {
	values_[rowIdx] = (T)(val);

	//__builtin_prefetch(values_ + rowIdx + kPrefetchDistance/sizeof(T));
}

template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::ShiftValueLeft(uint32_t nbits, uint64_t rowIdx) {
	values_[rowIdx] <<= nbits;
}

template <>
inline void __attribute__((always_inline))
SpecialiedColumn<uint8_t>::ShiftLeftBatch(uint64_t num_rows, uint32_t nbits) {

	//unused
}

template <>
inline void __attribute__((always_inline))
SpecialiedColumn<uint16_t>::ShiftLeftBatch(uint64_t num_rows, uint32_t nbits) {
	assert(0 == ((uint64_t)values_ & 31));
	const uint32_t banksize = 16;
	uint64_t curIdx = 0;
	__m256i reg;
	while ((curIdx + banksize - 1) < num_rows) {
		reg = _mm256_load_si256((__m256i const *) (values_ + curIdx));
		reg = _mm256_slli_epi16(reg, nbits);

		_mm256_stream_si256((__m256i *)(values_ + curIdx), reg);
		//_mm256_store_si256((__m256i *)(values_ + curIdx), reg);	//difference is not that much
		curIdx += banksize;
	}

	//shift remaining elements
	while (curIdx < num_rows) {
		values_[curIdx] <<= nbits;

		curIdx++;
	}
}

template <>
inline void __attribute__((always_inline))
SpecialiedColumn<uint32_t>::ShiftLeftBatch(uint64_t num_rows, uint32_t nbits) {
	assert(0 == ((uint64_t)values_ & 31));
	const uint32_t banksize = 8;
	uint64_t curIdx = 0;
	__m256i reg;
	while ((curIdx + banksize - 1) < num_rows) {
		reg = _mm256_load_si256((__m256i const *) (values_ + curIdx));
		reg = _mm256_slli_epi32(reg, nbits);

		_mm256_stream_si256((__m256i *)(values_ + curIdx), reg);
		curIdx += banksize;
	}

	//shift remaining elements
	while (curIdx < num_rows) {
		values_[curIdx] <<= nbits;

		curIdx++;
	}
}

template <>
inline void __attribute__((always_inline))
SpecialiedColumn<uint64_t>::ShiftLeftBatch(uint64_t num_rows, uint32_t nbits) {
	assert(0 == ((uint64_t)values_ & 31));
	const uint32_t banksize = 4;
	uint64_t curIdx = 0;
	__m256i reg;
	while ((curIdx + banksize - 1) < num_rows) {
		reg = _mm256_load_si256((__m256i const *) (values_ + curIdx));
		reg = _mm256_slli_epi64(reg, nbits);

		_mm256_stream_si256((__m256i *)(values_ + curIdx), reg);
		curIdx += banksize;
	}

	//shift remaining elements
	while (curIdx < num_rows) {
		values_[curIdx] <<= nbits;

		curIdx++;
	}
}

template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::LogicalOr(uint64_t orValue, uint64_t rowIdx) {
	values_[rowIdx] |= orValue;

	//__builtin_prefetch(values_ + rowIdx + kPrefetchDistance/sizeof(T));
}

template <class T>
inline uint64_t __attribute__((always_inline))
SpecialiedColumn<T>::GetValueAt(uint64_t rowIdx) const {
	return (uint64_t)values_[rowIdx];
}

template <class T>
inline uint64_t __attribute__((always_inline))
SpecialiedColumn<T>::GetValPortion(uint64_t rowIdx, uint32_t start, uint32_t end) const {
	T val = values_[rowIdx];
	val >>= (start - 1);
	val &= ((1ULL << (end - start + 1)) - 1);

	//__builtin_prefetch(values_ + rowIdx + kPrefetchDistance/sizeof(T));

	return (uint64_t)(val);
}

template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::AssembleFromColumn(SpecialiedColumn<uint8_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) {

	switch(sizeof(T)) {
		case 1:	//from uint8_t to uint8_t
		{
			//std::cout<< "uint8_t to uint8_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 16;	//NOTE: it has to transfer to uint16_t, so with degree=16
			uint64_t curIdx = 0;

			__m128i loadSrcReg;
			__m128i loadDstReg;
			__m256i fromReg;
			__m256i toReg;
			__m128i storeReg;
			uint8_t *fromValues = (uint8_t *)fromColumn->GetColumn();
			uint8_t fromSingleVal;

#if 1
			while ((curIdx + banksize - 1) < num_rows) {
				//load 16 elements from <fromColumn>
				loadSrcReg = _mm_load_si128((__m128i const *) (fromValues + curIdx));

				//convert to uint16_t
				fromReg = _mm256_cvtepu8_epi16(loadSrcReg);

				//do the right shifting
				fromReg = _mm256_srli_epi16(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi16(static_cast<int16_t>(mask)));

				//load 16 elements from thisColumn
				loadDstReg =  _mm_load_si128((__m128i const *) (values_ + curIdx));

				//convert to uint16_t
				toReg = _mm256_cvtepu8_epi16(loadDstReg);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				storeReg = _mm_packus_epi16(_mm256_extracti128_si256(toReg, 0),
						_mm256_extracti128_si256(toReg, 1));

				//write it back to the memory
				_mm_stream_si128((__m128i *)(values_ + curIdx), storeReg);
				curIdx += banksize;
			}
#endif

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 2://from uint8_t to uint16_t
		{
			//std::cout<< "uint8_t to uint16_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 16;
			uint64_t curIdx = 0;

			__m128i loadSrcReg;
			//__m128i loadDstReg;
			__m256i fromReg;
			__m256i toReg;
			//__m128i storeReg;
			uint8_t *fromValues = (uint8_t *)fromColumn->GetColumn();
			uint8_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 16 elements from <fromColumn>
				loadSrcReg = _mm_load_si128((__m128i const *) (fromValues + curIdx));

				//convert to uint16_t
				fromReg = _mm256_cvtepu8_epi16(loadSrcReg);

				//do the right shifting
				fromReg = _mm256_srli_epi16(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi16(static_cast<int16_t>(mask)));

				//load 16 elements from thisColumn
				toReg =  _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 4://from uint8_t to uint32_t
		{
			//std::cout<< "uint8_t to uint32_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 8;
			uint64_t curIdx = 0;

			__m128i loadSrcReg;
			__m256i fromReg;
			__m256i toReg;

			uint8_t *fromValues = (uint8_t *)fromColumn->GetColumn();
			uint8_t fromSingleVal;

			//how to read 8 continues uint8_t elements
			//first _mm_loadl_epi64, then _mm256_cvtepu8_epi32
			//VERY IMPORTANT: big/little endianess may affect the correctness
			//important: <num_rows - 8> because of _mm_loadl_epi64
			while (((curIdx + banksize - 1) + 8) < num_rows) {
				//load 8 elements from <fromColumn>
				loadSrcReg = _mm_loadl_epi64((__m128i const *) (fromValues + curIdx));

				//convert to uint32_t
				fromReg = _mm256_cvtepu8_epi32(loadSrcReg);

				//do the right shifting
				fromReg = _mm256_srli_epi32(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi32(static_cast<int32_t>(mask)));

				//load 8 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 8://from uint8_t to uint64_t
		{
			//std::cout<< "uint8_t to uint64_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 8;
			uint64_t curIdx = 0;

			__m128i loadSrcReg;
			__m256i fromReg;
			__m256i toReg;

			uint8_t *fromValues = (uint8_t *)fromColumn->GetColumn();
			uint8_t fromSingleVal;

			//how to read 4 continues uint8_t elements?
			//still _mm_loadl_epi64 (only use lowest 32b), then _mm256_cvtepu8_epi64
			//VERY IMPORTANT: big/little endianess may affect the correctness
			//important: <num_rows - 12> because of _mm_loadl_epi64
			while (((curIdx + banksize - 1) + 12) < num_rows) {
				//load 4 elements from <fromColumn>
				loadSrcReg = _mm_loadl_epi64((__m128i const *) (fromValues + curIdx));

				//convert to uint64_t
				fromReg = _mm256_cvtepu8_epi64(loadSrcReg);

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
	}

}

template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::AssembleFromColumn(SpecialiedColumn<uint16_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) {

	switch(sizeof(T)) {
		case 1:	//from uint16_t to uint8_t
		{
			//std::cout<< "uint16_t to uint8_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 16;
			uint64_t curIdx = 0;

			__m128i loadDstReg;
			__m256i fromReg;
			__m256i toReg;
			__m128i storeReg;
			uint16_t *fromValues = (uint16_t *)fromColumn->GetColumn();
			uint16_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 16 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi16(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi16(static_cast<int16_t>(mask)));

				//load 16 elements from thisColumn
				loadDstReg =  _mm_load_si128((__m128i const *) (values_ + curIdx));

				//convert to uint16_t
				toReg = _mm256_cvtepu8_epi16(loadDstReg);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				storeReg = _mm_packus_epi16(_mm256_extracti128_si256(toReg, 0),
						_mm256_extracti128_si256(toReg, 1));

				//write it back to the memory
				_mm_stream_si128((__m128i *)(values_ + curIdx), storeReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 2:	//from uint16_t to uint16_t
		{
			//std::cout<< "uint16_t to uint16_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 16;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint16_t *fromValues = (uint16_t *)fromColumn->GetColumn();
			uint16_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 16 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi16(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi16(static_cast<int16_t>(mask)));

				//load 16 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 4://from uint16_t to uint32_t
		{
			//std::cout<< "uint16_t to uint32_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 8;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint16_t *fromValues = (uint16_t *)fromColumn->GetColumn();
			uint16_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 8 elements from <fromColumn>
				__m128i fromReg_short = _mm_load_si128((__m128i const *) (fromValues + curIdx));
				fromReg = _mm256_cvtepu16_epi32 (fromReg_short);

				//do the right shifting
				fromReg = _mm256_srli_epi32(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi32(static_cast<int32_t>(mask)));

				//load 8 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 8://from uint16_t to uint64_t
		{
			//std::cout<< "uint16_t to uint64_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 4;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint16_t *fromValues = (uint16_t *)fromColumn->GetColumn();
			uint16_t fromSingleVal;

			//how to read 4 continues uint16_t elements and store in __m256i???
			//first _mm_loadl_epi64, then _mm256_cvtepu16_epi64
			//use _mm_set_epi16 in the worst case
			//important: <num_rows - 4> because of _mm_loadl_epi64
			while (((curIdx + banksize - 1) + 4) < num_rows) {
				//load 4 elements from <fromColumn>
				__m128i fromReg_short = _mm_loadl_epi64((__m128i const *) (fromValues + curIdx));
				fromReg = _mm256_cvtepu16_epi64 (fromReg_short);

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
	}

}

template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::AssembleFromColumn(SpecialiedColumn<uint32_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) {

	switch(sizeof(T)) {
		case 1:	//from uint32_t to uint8_t
		{
			//std::cout<< "from uint32_t to uint8_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 8;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint32_t *fromValues = (uint32_t *)fromColumn->GetColumn();
			uint32_t fromSingleVal;

			//scalar store back
			//important: <num_rows - 8> because of _mm_loadl_epi64
			while (((curIdx + banksize - 1) + 8) < num_rows) {
				//load 8 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi32(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi32(static_cast<int32_t>(mask)));

				//load 8 elements from thisColumn
				__m128i loadDstReg = _mm_loadl_epi64((__m128i const *) (values_ + curIdx));
				toReg = _mm256_cvtepu8_epi32(loadDstReg);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//scalar store back
				values_[curIdx] = (uint8_t)_mm256_extract_epi32 (toReg, 0);
				values_[curIdx + 1] = (uint8_t)_mm256_extract_epi32 (toReg, 1);
				values_[curIdx + 2] = (uint8_t)_mm256_extract_epi32 (toReg, 2);
				values_[curIdx + 3] = (uint8_t)_mm256_extract_epi32 (toReg, 3);
				values_[curIdx + 4] = (uint8_t)_mm256_extract_epi32 (toReg, 4);
				values_[curIdx + 5] = (uint8_t)_mm256_extract_epi32 (toReg, 5);
				values_[curIdx + 6] = (uint8_t)_mm256_extract_epi32 (toReg, 6);
				values_[curIdx + 7] = (uint8_t)_mm256_extract_epi32 (toReg, 7);

				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 2:	//from uint32_t to uint16_t
		{
			//std::cout<< "from uint32_t to uint16_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 8;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint32_t *fromValues = (uint32_t *)fromColumn->GetColumn();
			uint32_t fromSingleVal;

			//how to shink __m256i to __m128i
			//_mm_packus_epi32
			while ((curIdx + banksize - 1) < num_rows) {
				//load 8 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi32(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi32(static_cast<int32_t>(mask)));

				//load 8 elements from thisColumn
				__m128i toReg_short = _mm_load_si128((__m128i const *) (values_ + curIdx));
				toReg = _mm256_cvtepu16_epi32(toReg_short);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//convert back to __m128i
				__m128i toReg_write = _mm_packus_epi32(_mm256_extracti128_si256(toReg, 0), _mm256_extracti128_si256(toReg, 1));

				//write it back to the memory
				_mm_stream_si128((__m128i *)(values_ + curIdx), toReg_write);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 4://from uint32_t to uint32_t
		{
			//std::cout<< "from uint32_t to uint32_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 8;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint32_t *fromValues = (uint32_t *)fromColumn->GetColumn();
			uint32_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 8 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi32(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi32(static_cast<int32_t>(mask)));

				//load 8 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 8://from uint32_t to uint64_t
		{
			//std::cout<< "from uint32_t to uint64_t" << std::endl;
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 4;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint32_t *fromValues = (uint32_t *)fromColumn->GetColumn();
			uint32_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 4 elements from <fromColumn>
				__m128i fromReg_short = _mm_load_si128((__m128i const *) (fromValues + curIdx));
				fromReg = _mm256_cvtepu32_epi64(fromReg_short);

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//write it back to the memory
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);
				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
	}

}

template <class T>
inline void __attribute__((always_inline))
SpecialiedColumn<T>::AssembleFromColumn(SpecialiedColumn<uint64_t> *fromColumn, uint32_t start, uint32_t end, uint64_t num_rows, uint32_t originWidth) {

	switch(sizeof(T)) {
		case 1:	//from uint64_t to uint8_t
		{
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 4;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint64_t *fromValues = (uint64_t *)fromColumn->GetColumn();
			uint64_t fromSingleVal;

			while (((curIdx + banksize - 1) + 12) < num_rows) {
				//load 4 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				__m128i loadDstReg = _mm_loadl_epi64((__m128i const *) (values_ + curIdx));
				toReg = _mm256_cvtepu8_epi64(loadDstReg);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//scalar store back;
				values_[curIdx] = (uint8_t)_mm256_extract_epi64 (toReg, 0);
				values_[curIdx + 1] = (uint8_t)_mm256_extract_epi64 (toReg, 1);
				values_[curIdx + 2] = (uint8_t)_mm256_extract_epi64 (toReg, 2);
				values_[curIdx + 3] = (uint8_t)_mm256_extract_epi64 (toReg, 3);

				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 2:	//from uint64_t to uint16_t
		{
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 4;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint64_t *fromValues = (uint64_t *)fromColumn->GetColumn();
			uint64_t fromSingleVal;

			// no way to convert 64-bank size to 16-bank size??
			while (((curIdx + banksize - 1) + 4) < num_rows) {
				//load 4 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				__m128i toReg_short = _mm_loadl_epi64((__m128i const *) (values_ + curIdx));
				toReg = _mm256_cvtepu16_epi64(toReg_short);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//store back; currently use sequential store; any better solution???
				values_[curIdx] = (uint16_t)_mm256_extract_epi64 (toReg, 0);
				values_[curIdx + 1] = (uint16_t)_mm256_extract_epi64 (toReg, 1);
				values_[curIdx + 2] = (uint16_t)_mm256_extract_epi64 (toReg, 2);
				values_[curIdx + 3] = (uint16_t)_mm256_extract_epi64 (toReg, 3);

				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 4://from uint64_t to uint32_t
		{
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 4;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint64_t *fromValues = (uint64_t *)fromColumn->GetColumn();
			uint64_t fromSingleVal;

			// no way to convert 64-bank size to 32-bank size??
			while ((curIdx + banksize - 1) < num_rows) {
				//load 4 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				__m128i toReg_short = _mm_load_si128((__m128i const *) (values_ + curIdx));
				toReg = _mm256_cvtepu32_epi64(toReg_short);

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//store back; currently use sequential store; any better solution???
				values_[curIdx] = (uint32_t)_mm256_extract_epi64 (toReg, 0);
				values_[curIdx + 1] = (uint32_t)_mm256_extract_epi64 (toReg, 1);
				values_[curIdx + 2] = (uint32_t)_mm256_extract_epi64 (toReg, 2);
				values_[curIdx + 3] = (uint32_t)_mm256_extract_epi64 (toReg, 3);

				curIdx += banksize;
			}

			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);
				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
		case 8://from uint64_t to uint64_t
		{
			uint64_t mask = (1ULL << (end - start + 1)) - 1;
			const uint32_t banksize = 4;
			uint64_t curIdx = 0;
			__m256i fromReg;
			__m256i toReg;
			uint64_t *fromValues = (uint64_t *)fromColumn->GetColumn();
			uint64_t fromSingleVal;

			while ((curIdx + banksize - 1) < num_rows) {
				//load 4 elements from <fromColumn>
				fromReg = _mm256_load_si256((__m256i const *) (fromValues + curIdx));

				//do the right shifting
				fromReg = _mm256_srli_epi64(fromReg, (originWidth-end));

				//do the masking
				fromReg = _mm256_and_si256(fromReg, _mm256_set1_epi64x(static_cast<long long>(mask)));

				//load 4 elements from thisColumn
				toReg = _mm256_load_si256((__m256i const *) (values_ + curIdx));

				//do the logical OR
				toReg = _mm256_or_si256(toReg, fromReg);

				//store back;
				_mm256_stream_si256((__m256i *)(values_ + curIdx), toReg);

				curIdx += banksize;
			}


			//operate on remaining elements
			while (curIdx < num_rows) {
				fromSingleVal = fromValues[curIdx];

				fromSingleVal >>= (originWidth-end);

				fromSingleVal &= mask;

				values_[curIdx] |= fromSingleVal;

				curIdx++;
			}
		}
		break;
	}

}


}

#endif
