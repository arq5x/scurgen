from matplotlib import pyplot as plt
import matplotlib

def debug_plot(h):
    imshow_kwargs = dict(interpolation='nearest', 
            cmap=matplotlib.cm.Spectral_r)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mappable = ax.imshow(h.matrix, **imshow_kwargs)
    plt.colorbar(mappable)
    x, y, labels = h.curve()
    ax.plot(x, y, 'k')
    for i in range(10):
        ax.text(x[i], y[i], labels[i])
    ax.axis('tight')
    plt.show()
