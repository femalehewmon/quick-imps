import sys
import pylab as pl
import math


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
        z = [1.0*t_sum_arr[i]/t_arr[i] for i in range(1, len(x))]
        y = [1.0*t_arr[i]/t_arr[i-1] for i in range(1, len(x))]
        v = [1.0*math.sqrt(t_arr[i]) for i in range(1, len(x))]
        vv = [1.0*math.sqrt(v[i-1]) for i in range(1, len(x))]
        print("t: " + str(x))
        print("x: " + str(t_arr))
        print("sum_prev: " + str(t_sum_arr))
        print("sum_prev/x: " + str(z))
        print("x/x_prev: " + str(y))
        print("sqrt(x): " + str(v))
        print("sqrt(x) + sqrt(sqrt(x)): " + str(vv))
        #y = [0.5*pow(t_arr[i],2) for i in range(1, len(x))]
        #pl.plot(x[1:], y)
        #print(y)
    return plt1, plt2

#small = list(x) # 2^t
#orig = map(lambda t: pow(2, t), small) # 2^2^t
#larger = map(lambda t: pow(2, t), orig) # 2^2^2^t
#splot, = pl.plot(x, small)
#oplot, = pl.plot(x, orig)
#lplot, = pl.plot(x, larger)

#small_sum = sum_arr(small)
#orig_sum = sum_arr(orig)
#splot_sum = pl.bar(x, small_sum, alpha=0.5)
#oplot_sum = pl.bar(x, orig_sum, alpha=0.5)

def main(num_iters, plot_type=1):
    x = []
    x.extend(range(1, num_iters+1))

    if plot_type == '1':
        print "plottype 1"
        splot,_ = plot_t(x, lambda t: t)
        oplot,_ = plot_t(x, lambda t: pow(2, t))
        lplot,_ = plot_t(x, lambda t: pow(2, pow(2,t)))
        pl.legend([splot, oplot, lplot],
                  ("2^t", "2^2^t", "2^2^2^t"), 'best', numpoints=1)
    elif plot_type == '2':
        print "plottype 2"
        splot,_ = plot_t(x, lambda t: t)
        oplot,_ = plot_t(x, lambda t: pow(2, t))
        pl.legend([splot, oplot], ("2^t", "2^2^t"), 'best', numpoints=1)
    elif plot_type == '3':
        print "plottype 3"
        splot,splot_sum = plot_t(x, lambda t: t, True)
        oplot,oplot_sum = plot_t(x, lambda t: pow(2, t), True)
        pl.legend([splot, splot_sum, oplot, oplot_sum],
                  ("2^t", "2^t sum of prev", "2^2^t", "2^2^t sum of prev"),
                  'best', numpoints=1)
    elif plot_type == 'small':
        splot,splot_sum = plot_t(x, lambda t: t, True)
        pl.legend([splot, splot_sum], ("2^t", "2^t sum of prev"), loc=2)
    elif plot_type == 'orig':
        oplot,oplot_sum = plot_t(x, lambda t: pow(2, t), True)
        pl.legend([oplot, oplot_sum], ("2^2^t", "2^2^t sum of prev"), loc=2)
    elif plot_type == 'large':
        lplot,lplot_sum = plot_t(x, lambda t: pow(2, pow(2,t)), True)
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
