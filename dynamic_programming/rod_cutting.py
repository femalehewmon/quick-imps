import os
import sys

run_counts = 0

# CLRS pg. 369
# Return array of lengths to cut and maximum revenue
def extended_bottom_up_cut_rod(p, n):
    global run_counts
    r = ["will be overwritten"] * (n + 1)
    s = ["will be overwritten"] * (n + 1)
    r[0] = 0
    for j in range(1, n+1):
        q = float("-inf")
        for i in range(1, j+1):
            if q < p[i] + r[j - 1]:
                q = p[i] + r[j-1]
                s[j] = i
                run_counts += 1
        r[j] = q
    return r, s

# CLRS pg. 366
# THETA(n^2) due to doubly nested for loop
def bottom_up_cut_rod(p, n):
    global run_counts
    r = ["will be overwritten"] * (n + 1)
    r[0] = 0
    for j in range(1, n + 1):
        q = float("-inf")
        for i in range(1, j+1):
            q = max(q, p[i] + r[j - i])
            run_counts += 1
        r[j] = q
    return r[n]

# CLRS pg. 365
# THETA(n^2) running time due to for loop running n times all n calls
def memoized_cut_rod(p, n, r):
    global run_counts
    run_counts += 1
    if r[n] >= 0:
        return r[n]
    q = 0
    if n != 0:
        q = float("-inf")
        for i in range(1,n+1):
            q = max(q, p[i] + memoized_cut_rod(p, n-i, r))
    r[n] = q
    return q

# CLRS pg. 363
# reason to use dynamic programming
# this algorithm re-computes the same subproblem over and over again
def cut_rod(p, n):
    global run_counts
    run_counts += 1
    if n == 0:
        return 0
    q = float("-inf")
    for i in range(1,n+1):
        q = max(q, p[i] + cut_rod(p, n-i))
    return q

def main():
    global run_counts
    # length i  |  1  2  3  4  5   6   7   8   9   10
    # price_array per length cut, i
    price_array = [0, 1, 5, 8, 9, 10, 17, 17, 20, 24, 30]
    price_array += price_array
    # length of rod
    n = 10

    # The following example is a demonstration of why dynamic programming is so important
    print("No dynamic programming")
    run_counts = 0
    max_price = cut_rod(price_array, n)
    print max_price
    print "run counts " + str(run_counts)

    # Dynamic programming uses additional memory to save computation time: a time-memory trade-off
    # A dynamic-programming approach runs in polynomial time when the number of distinct subproblems involved is polynomial in the input size and we can solve each such subproblem in polynomial time
    # 2 different ways to do dynamic programming:
    #    - top-down with memoization
    #         o written recursively in a natural manner
    #         o the result of each subproblem is saved
    #         o depth-first search
    #    - bottom-up method
    #         o typically depends on the "size" of the subproblem
    #         o sort the subproglems by size, solving smallest first
    #         o can take linear space
    #         o breadth-first search
    # Both methods yield algorithms with the same asymptotic running time, except in unusual circumstances where top-down does not recurse on all subproblems
    # Bottom-up often has much better constant factors since it has less overhead for procedure calls

    # memoization added for rod cutting
    print("Memoized")
    run_counts = 0
    r = [float("-inf")] * (n + 1) # space to save computed results
    max_price = memoized_cut_rod(price_array, n, r)
    print max_price
    print "run counts " + str(run_counts)

    # bottom up rod cutting
    print("Bottom up")
    run_counts = 0
    max_price = bottom_up_cut_rod(price_array, n)
    print max_price
    print "run counts " + str(run_counts)

    # What if we want to know the actual answer of how long of cuts we should make and not just the total revenue? The following algorithm returns both pieces of information.
    # The below code reflects PRINT_CUT_ROD_SOLUTION in CLRS
    run_counts = 0
    max_price, s = extended_bottom_up_cut_rod(
        price_array, n)
    print max_price
    print s
    k = n
    while k > 0:
        print s[k]
        k -= s[n]
    print "run counts " + str(run_counts)

main()
