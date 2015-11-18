import os
import sys

run_count = 0

def LongestPalindromeRecursive(x, i, j, buff=""):
    global run_count
    run_count += 1
    if(j >= i):
        #print("{3} {0} {1} {2}".format(x[i], x[j], x[i:j+1], buff))
        if(x[i] == x[j]):
            point = 2
            if(i == j):
                point = 1
            return point + LongestPalindromeRecursive(x, i+1, j-1, "")
        else:
            return max(
                LongestPalindromeRecursive(x, i+1, j, "   "),
                LongestPalindromeRecursive(x, i, j-1, "      ")
            )
    return 0

table = {}
def LongestPalindromeMemoization(x, i, j, buff=""):
    global run_count
    run_count += 1
    if(j >= i):
        #print("{3} {0} {1} {2}".format(x[i], x[j], x[i:j+1], buff))
        if(x[i] == x[j]):
            point = 2
            if(i == j):
                point = 1
            return point + LongestPalindrome(x, i+1, j-1, "")
        else:
            return max(
                LongestPalindrome(x, i+1, j, "   "),
                LongestPalindrome(x, i, j-1, "      ")
            )
    return 0

def LongestPalindrome(x, i, j, buff=""):
    if((i, j) in table):
        return table[(i, j)]
    else:
        value = LongestPalindromeMemoization(x, i, j, buff)
        table[(i, j)] = value
        return value

# Running time = O(n^2) given the two for loops
# However, the 2nd for loop is decreasing by a constant number each time
# This works out to taking O(((n^2)/2) + n) time and space
# Therefore, it remains O(n^2), while taking O(n^2) space
DIAG = "diag"
SING = "sing"
def LongestPalindromeDynamicProgramming(x):
    global run_count

    # create matrix with dummy '0' row and column
    matrix = [[0]] * (len(x)+2)
    matrix[0] = [0] * (len(x)+1)
    backpointers = [[(0, None)]] * (len(x)+2)
    backpointers[0] = [[(0, None)]] * (len(x)+1)

    for i in range(1, len(x) + 2):
        leftIdx = len(x) - i
        matrix[i] = [0] * (len(x) + 2 - i)
        backpointers[i] = [(0, None)] * (len(x) + 2 - i)
        for j in range(1, len(x) + 2 - i):
            rightIdx = j - 1
            run_count += 1
            if leftIdx == rightIdx:
                matrix[i][j] = matrix[i-1][j-1] + 1
                backpointers[i][j] = (matrix[i][j], SING)
            elif x[leftIdx] == x[rightIdx]:
                matrix[i][j] = matrix[i-1][j-1] + 2
                backpointers[i][j] = (matrix[i][j], DIAG)
            else:
                matrix[i][j] = max(matrix[i-1][j], matrix[i][j-1])
                backpointers[i][j] = (matrix[i][j], None)

    # remove dummy values from matrix
    #matrix = matrix[1:] # remove starting row
    #matrix = [row[1:] for row in matrix] # remove starting column
    matrix = matrix[:-1] # remove ending row
    backpointers = backpointers[:-1] # remove ending row

    cols = map(len, matrix)
    total_size = sum(cols[1:]) - (len(cols) - 1)
    #total_size = sum(cols)
    print("total size of matrix {0}".format(total_size))
    max_values = map(max, matrix)
    return max(max_values), matrix, backpointers

# Reverse output and print, skipping first value if LP is odd
def PrintLongestPalindrome(string, matrix, row_idx, col_idx):
    global run_count
    if row_idx == 0 or col_idx == 0:
        return
    run_count += 1
    current_val = matrix[row_idx][col_idx][0]
    current_pos = matrix[row_idx][col_idx][1]
    above_left_val = matrix[row_idx-1][col_idx-1][0]
    #print(
    #   "row {0} col {1} value {2}".format(row_idx, col_idx, current_val))
    if(current_pos == DIAG):
        print("diag found " + string[row_idx - 1])
        # on a diagonal with matching letter
        return PrintLongestPalindrome(string, matrix, row_idx-1, col_idx-1)
        #print(string[col_idx])
        #print("to end " + string[col_idx])
    if(current_pos == SING):
        print("single found " + string[row_idx - 1])
        # on a diagonal of matching letter to itself
        return PrintLongestPalindrome(string, matrix, row_idx-1, col_idx-1)
        #print(string[col_idx])
    else:
        above_val = matrix[row_idx - 1][col_idx][0]
        left_val = matrix[row_idx][col_idx - 1][0]
        max_above_left = max(above_val, left_val)
        if above_val == max_above_left:
            return PrintLongestPalindrome(
                string, matrix, row_idx - 1, col_idx)
        else:
            return PrintLongestPalindrome(
                string, matrix, row_idx, col_idx - 1)


def main():
    global run_count
    value = "XABCBBACXA"
    print("Total length of value: {0}".format(len(value)))

    startIdx = 0
    endIdx = len(value) - 1

    print("\nLongestPalindromeRecursive")
    run_count = 0
    #print("{0}".format(
    #    LongestPalindromeRecursive(value, startIdx, endIdx)))
    print("run_count " + str(run_count))

    print("\nLongestPalindromeMemoization")
    run_count = 0
    print("{0}".format(
        LongestPalindromeMemoization(value, startIdx, endIdx)))
    print("run_count " + str(run_count))

    print("\nLongestPalindromeDynamicProgramming")
    run_count = 0
    LP_len, matrix, bckpntrs = LongestPalindromeDynamicProgramming(value)
    #print(matrix)
    max_values = map(max, matrix)
    max_value = max(max_values)
    max_row_index = max_values.index(max_value)
    max_col_index = matrix[max_row_index].index(max_value)
    print("{0}".format(LP_len))
    print("run_count " + str(run_count))

    print("\nPrintLongestPalindrome")
    run_count = 0
    PrintLongestPalindrome(value[::-1], bckpntrs, max_row_index, max_col_index)
    print("run_count " + str(run_count))

main()
