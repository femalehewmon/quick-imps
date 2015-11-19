import math
import numpy.random as nprnd

def verify_heap(heap):
    heap = [-1] + heap
    for idx in range(2, len(heap)):
        parent_idx = int(math.floor(idx/2))
        assert(heap[parent_idx] >= heap[idx])

def NlogNHeapify(unsorted):
    heap = [-1]
    for idx, val in enumerate(unsorted,1):
        heap.append(val)
        curr_idx = idx
        while curr_idx > 1:
            parent_idx = int(math.floor(curr_idx/2))
            if heap[curr_idx] > heap[parent_idx]:
                parent = heap[parent_idx]
                heap[parent_idx] = heap[curr_idx]
                heap[curr_idx] = parent
                curr_idx = parent_idx
            else:
                break
    return heap[1:]

def NHeapify(unsorted):
    heap = [-1] + unsorted
    for idx in range(len(heap)/2, 0, -1):
        heap = MaxHeapify(heap, idx)
    return heap[1:]

def MaxHeapify(heap, idx):
    left = idx * 2
    right = left + 1
    largest = idx
    if left < len(heap) and heap[left] > heap[largest]:
        largest = left
    if right < len(heap) and heap[right] > heap[largest]:
        largest = right
    if largest != idx:
        hold = heap[largest]
        heap[largest] = heap[idx]
        heap[idx] = hold
        heap = MaxHeapify(heap, largest)
    return heap

def main():
    unsorted = nprnd.randint(100, size = 5)
    unsorted = unsorted.tolist()
    print(unsorted)
    print "-----------------------------------"

    print "O(nlogn) heapify"
    nlogn_heap = NlogNHeapify(unsorted)
    print(str(nlogn_heap))
    verify_heap(nlogn_heap)
    print "-----------------------------------"

    print "O(n) heapify"
    n_heap = NHeapify(unsorted)
    print(str(n_heap))
    verify_heap(n_heap)
    print "-----------------------------------"

    print "Are the two methods equal?"
    if n_heap == nlogn_heap:
        print "Yes!"
    else:
        print "No."

if __name__ == "__main__":
    main()
