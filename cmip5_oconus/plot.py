
import numpy as np
import cartopy
import cartopy.crs as ccrs

from matplotlib.colors import Normalize

# to standardize plots
PROJECTION = ccrs.Mercator()


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def custom_div_cmap(numcolors=11, name='custom_div_cmap',
                    mincol='blue', midcol='white', maxcol='red'):
    """ Create a custom diverging colormap with three colors

    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    https://xkcd.com/color/rgb/
    """

    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list(name=name,
                                             colors =[mincol, midcol, maxcol],
                                             N=numcolors)
    return cmap


def make_plot(da, ax=None, **plot_kwargs):

    if ax is None:
        ax = plt.gca()

    # make the actual plot
    da.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), **plot_kwargs)

    # geometries
    ax.coastlines();
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor='lightgray')

#     gl = ax.gridlines(draw_labels=True, linewidth=0)
#     gl.xlabels_top = False
#     gl.ylabels_left = False

    # boundaries
    states_provinces = cartopy.feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='gray')
    return ax


def add_ylabel(ax, text, fontsize=12):
    return ax.text(-0.07, 0.55, text, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes, fontsize=fontsize)
