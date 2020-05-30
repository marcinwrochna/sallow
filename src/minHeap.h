#pragma once
#include <cassert>
#include <type_traits>
#include <vector>

#include "graph.h"

/** Min heap of n elements identified by VertexID in 1..n. */
template<typename Value = int>
class MinHeap {
    private:
        size_t nElements;
        std::vector<int> positions; // Map from element id to position in heap.
        struct Node { VertexID id; Value value; };
        std::vector<Node> heap;
        // Position 0 is dummy, parent of i is i/2, children are 2i, 2i+1.
        // 1 is root, leafs are [heap.size()/2, heap.size()).
        // Positions [initVector.size()+1, heap.size()) are dummy.
        const VertexID dummyId = 0;
        const Value dummyValue = std::numeric_limits<Value>::max();

    public:
        /** Initialize heap with element i having value initVector[i-1].*/
        template <typename RandGen = void>
        MinHeap(const std::vector<Value>& initVector, RandGen* randGen = nullptr)
            : nElements(initVector.size()),
              positions(nElements + 1),
              heap(heapAllocSize(nElements)) {
            positions[dummyId] = 0;
            heap[0] = Node{dummyId, dummyValue};
            for (VertexID i = 1; i <= nElements; ++i)
                positions[i] = i;
            if constexpr (!std::is_void_v<RandGen>)
                std::shuffle(positions.begin() + 1, positions.end(), *randGen);
            for (VertexID i = 1; i <= nElements; ++i)
                heap[positions[i]] = Node{i, initVector[i - 1]};
            for (gsl::index i = nElements + 1; i < heap.size(); ++i)
                heap[i] = Node{dummyId, dummyValue};
            for (gsl::index i = nElements / 2; i >= 1; --i)
                increaseValue(heap[i].id, heap[i].value);
        }


        std::pair<VertexID, Value> top() {
            assert(nElements);
            return std::make_pair(heap[1].id, heap[1].value);
        }


        /** Remove element with smallest value and return its id. */
        std::pair<VertexID, Value> pop() {
            assert(nElements);
            --nElements;
            std::pair<VertexID, Value> result = std::make_pair(heap[1].id, heap[1].value);
            assert(heap[1].id != dummyId);
            positions[heap[1].id] = 0;
            int pos = 1;

            while (2 * pos < heap.size()) {
                // Invariant: heap[pos] is to be overwritten by a child or deleted.
                assert(2 * pos + 1 < heap.size());
                if (heap[2 * pos].value <= heap[2 * pos + 1].value) {
                    heap[pos] = std::move(heap[2 * pos]);
                    positions[heap[pos].id] = pos;
                    pos = 2 * pos;
                } else {
                    heap[pos] = std::move(heap[2 * pos + 1]);
                    positions[heap[pos].id] = pos;
                    pos = 2 * pos + 1;
                }
            }

            // Make the last node in the chain a dummy.
            heap[pos].id = dummyId;
            heap[pos].value = dummyValue;

            return result;
        }

        void increaseValue(VertexID id, Value value) {
            int pos = positions[id];
            assert(pos >= 1 && pos < heap.size());
            assert(pos < positions.size());

            while (2 * pos < heap.size()) {
                assert(2 * pos + 1 < heap.size());
                // Invariant: heap[pos] is to be overwritten (by a child or by (id,key)).
                if (heap[2 * pos].value < value) {
                    if (heap[2 * pos].value < heap[2 * pos + 1].value) {
                        heap[pos] = heap[2 * pos];
                        positions[heap[pos].id] = pos;
                        pos = 2 * pos;
                    } else {
                        heap[pos] = heap[2 * pos + 1];
                        positions[heap[pos].id] = pos;
                        pos = 2 * pos + 1;
                    }
                } else {
                    if (heap[2 * pos + 1].value < value) {
                        heap[pos] = heap[2 * pos + 1];
                        positions[heap[pos].id] = pos;
                        pos = 2 * pos + 1;
                    } else
                        break;
                }
            }

            heap[pos].id = id;
            heap[pos].value = value;
            positions[id] = pos;
        }

        void decreaseValue(VertexID id, Value value) {
            int pos = positions[id];
            assert(pos >= 1 && pos < heap.size());

            while (pos / 2) {
                // Invariant: heap[pos] is to be overwritten (by a parent or by (id,key)).
                if (heap[pos / 2].value <= value)
                    break;
                heap[pos] = heap[pos / 2];
                positions[heap[pos].id] = pos;
                pos /= 2;
            }

            heap[pos].id = id;
            heap[pos].value = value;
            positions[id] = pos;
        }

        void changeValueBy(VertexID id, Value amount) {
            int pos = positions[id];
            if (amount > 0)
                increaseValue(id, heap[pos].value + amount);
            else if (amount < 0)
                decreaseValue(id, heap[pos].value + amount);
        }

        void tryDecrementValue(VertexID id) {
            int pos = positions[id];
            if (!pos)
                return;
            decreaseValue(id, heap[pos].value - 1);
        }

        void maximizeValue(VertexID id, Value value) {
            int pos = positions[id];
            if (heap[pos].value < value)
                increaseValue(id, value);
        }

        bool contains(VertexID id) const {
            return (bool)positions[id];
        }

        Value getValue(VertexID id) const {
            int pos = positions[id];
            assert(pos);
            return heap[pos].value;
        }

        size_t size() const {
            return nElements;
        }

        bool empty() const {
            return nElements == 0;
        }

    protected:
        /** Return first power of two strictly larger than n. */
        constexpr static size_t heapAllocSize(size_t n) {
            n |= (n >> 1);
            n |= (n >> 2);
            n |= (n >> 4);
            n |= (n >> 8);
            n |= (n >> 16);
            n++;
            assert(n);
            return n;
        }
};
