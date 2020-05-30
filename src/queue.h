#pragma once

/** A FIFO queue with a bound on the total number of pushes. */
template <typename Value>
class Queue
{
 public:
    Queue() : _maxSize(0), mem(nullptr), begin(nullptr), end(nullptr) {}
    Queue(size_t maxSize) : _maxSize(maxSize), mem(reinterpret_cast<Value*>(std::calloc(maxSize, sizeof(Value)))), begin(mem), end(mem) {}
    ~Queue() { if (mem) std::free(mem); }
    void reset(size_t maxSize) { if (mem) std::free(mem); _maxSize = maxSize; mem = reinterpret_cast<Value*>(std::calloc(maxSize, sizeof(Value))); begin = end = mem; }
    void push(Value v) { assert(mem && end < mem + _maxSize); *(end++) = v; }
    Value pop() { assert(mem && begin < end); return *(begin++); }
    bool empty() const { return begin == end; }
    void clear() { begin = end = mem; }
    size_t size() const { return end - begin; }
    size_t capacity() const { return _maxSize; }
    Value front() const { assert(mem && begin < end); return *begin; }
    Value back() const { assert(mem && begin < end); return *(end - 1); }
 private:
    size_t _maxSize;
    Value* mem;
    Value* begin;
    Value* end;
};