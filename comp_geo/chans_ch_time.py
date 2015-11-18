import sys
import pylab as pl

'''
Plots the time complexity overheard of different guess steps for
Chan's point set convex hull algorithm.
See: http://cgm.cs.mcgill.ca/~athens/cs507/ChanHull.ps
Guesses and time complexities:
    m = 2^t,     O(nlog^2m)
    m = 2^2^t,   O(nlogm)
    m = 2^2^2^t, O(nlogn)
'''

def sum_arr(arr):
    sum_arr = map(lambda(i,t): sum(arr[0:i]), enumerate(arr))
    return sum_arr

def plot_t(x, func, plot_sum=False):
    t_arr = map(func, x)
    plt1, = pl.plot(x, t_arr)
    plt2 = None
    if plot_sum:
        col = plt1.get_color()
        t_sum_arr = sum_arr(t_arr)
        plt2 = pl.bar(x, t_sum_arr, color=col, alpha=0.5, width=1)
    return plt1, plt2

def main(num_iters, plot_type=1):
    x = []
    x.extend(range(1, num_iters+1))

    small_guess = lambda t: t
    orig_guess = lambda t: pow(2, t)
    large_guess = lambda t: pow(2, pow(2, t))

    if plot_type == '1':
        splot,_ = plot_t(x, small_guess)
        oplot,_ = plot_t(x, orig_guess)
        lplot,_ = plot_t(x, large_guess)
        pl.legend([splot, oplot, lplot],
                  ("2^t", "2^2^t", "2^2^2^t"), 'best', numpoints=1)
    elif plot_type == '2':
        splot,_ = plot_t(x, small_guess)
        oplot,_ = plot_t(x, orig_guess)
        pl.legend([splot, oplot], ("2^t", "2^2^t"), 'best', numpoints=1)
    elif plot_type == '3':
        splot,splot_sum = plot_t(x, small_guess, True)
        oplot,oplot_sum = plot_t(x, orig_guess, True)
        pl.legend([splot, splot_sum, oplot, oplot_sum],
                  ("2^t", "2^t sum of prev", "2^2^t", "2^2^t sum of prev"),
                  'best', numpoints=1)
    elif plot_type == 'small':
        splot,splot_sum = plot_t(x, small_guess, True)
        pl.legend([splot, splot_sum], ("2^t", "2^t sum of prev"), loc=2)
    elif plot_type == 'orig':
        oplot,oplot_sum = plot_t(x, orig_guess, True)
        pl.legend([oplot, oplot_sum], ("2^2^t", "2^2^t sum of prev"), loc=2)
    elif plot_type == 'large':
        lplot,lplot_sum = plot_t(x, large_guess, True)
        pl.legend([lplot, lplot_sum], ("2^2^2^t", "2^2^2^t sum of prev"),
                  loc=2)

    pl.title("Chan's Algorithm, Overhead of Guesses")
    pl.xlabel("t")
    pl.ylabel("log(m)")
    pl.show()

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print "Missing args"
        sys.exit()
    num_iters = int(sys.argv[1])
    plot_type = sys.argv[2]
    main(num_iters, plot_type)
