from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.io import fits
from astropy.table import Table
from glob import glob
from icecream import ic
from matplotlib import pyplot as plt, ticker
import numpy as np
import os
from plottery.plotutils import savefig
from tqdm import tqdm
import warnings

from astro.clusters import ClusterCatalog

warnings.simplefilter('ignore')


def main():
    act = Table(fits.getdata(
        '~/Documents/actpol/catalogs/act-dr6_probable_confirmed.fits'))
    act = ClusterCatalog(
        'act-dr6', act, base_cols=('name','RADeg','decDeg','redshift'))
    ic(act['name','ra','dec'][act['name'] == 'ACT-CL J0841.3+1141'])

    # format_images(act, 2, 0.262, barsize=20*u.arcsec)
    # format_images(act, 3, 0.262, barsize=30*u.arcsec)
    # format_images(act, 5, 0.262, barsize=1*u.arcmin)
    # format_images(act, 10, 0.262, barsize=1*u.arcmin)
    format_images(act, 20, 0.262, barsize=2*u.arcmin)
    # format_images(act, 30, 0.524, barsize=5*u.arcmin, xyname=(0.04,0.93))
    return


def annotate_image(act, imgfile, size, scale, dz=0, barsize=1*u.arcmin,
                   xyname=(0.04,0.90), xybar=(0.85,0.06), barwidth=4,
                   barcolor='w', fontsize=16):
    """barsize is in arcmin"""
    output = imgfile.replace('images/', 'images/annotated/')
    if dz > 0:
        output = output.replace('.jpg', '_photoz.jpg')
    cluster = imgfile.split('/')[-1].split('__')[0].replace('_', ' ')
    cl = act[act['name'] == cluster][0]
    img = plt.imread(imgfile)
    imgsize = img.shape[0]
    scale = scale * u.arcsec / u.pix
    if barsize.unit in ('Mpc', 'kpc'):
        barname = f'{barsize.value} {barsize.unit}'
        barsize = barsize / cosmo.kpc_proper_per_arcmin(cl['redshift']).value
    elif barsize.unit == 'arcmin':
        barname = f"{barsize.value:.0f}'"
    elif barsize.unit == 'arcsec':
        barname = f'{barsize.value:.0f}"'
    barpix = (barsize / scale).to(u.pix).value
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, origin='lower')
    ax.annotate(f'{cluster} (z={cl["z"]:.3f})\nSNR = {cl["SNR"]:.1f}',
                xy=xyname, xycoords='axes fraction', ha='left',
                va='center', fontsize=fontsize, color=barcolor)
    xybar = (xybar[0] * imgsize, xybar[1] * imgsize)
    ax.plot([xybar[0]-barpix/2, xybar[0]+barpix/2], [xybar[1],xybar[1]],
            lw=barwidth, color=barcolor)
    ax.annotate(barname, (xybar[0], xybar[1]+0.03*imgsize), color=barcolor,
                fontsize=fontsize, ha='center', va='center')
    ax.axis('off')
    path = os.path.split(output)[0]
    os.makedirs(path, exist_ok=True)
    savefig(output, fig=fig, verbose=False)
    return output


def galaxy_contours(ax, cluster, dz):

    return


def format_images(act, size, scale, **kwargs):
    path = f'images/{size}_{scale}'
    imgfiles = sorted(glob(os.path.join(path, 'ACT-CL_J*.jpg')))
    for imgfile in tqdm(imgfiles):
        output = annotate_image(act, imgfile, size, scale, **kwargs)


if __name__ == "__main__":
    main()