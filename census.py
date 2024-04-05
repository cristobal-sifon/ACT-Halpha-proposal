from astLib.astCoords import decimal2hms, decimal2dms
from astroplan import FixedTarget, Observer
from astroplan import AltitudeConstraint, AirmassConstraint, AtNightConstraint
from astroplan import is_observable, is_always_observable, months_observable
from astroplan.plots import plot_airmass
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time
from icecream import ic
from matplotlib import pyplot as plt, ticker
import numpy as np
import os
from plottery.plotutils import savefig, update_rcParams

update_rcParams()

from astro.clusters import ClusterCatalog


def main():
    # act = ClusterCatalog('act-dr5')
    # there appear to be the same number in DR6, probably because the
    # spectroscopic redshifts have not been updated?
    act = Table(
        fits.getdata(
            "~/Documents/actpol/catalogs/act-dr6nemo_probable_confirmed_20240403.fits"
        )
    )
    act = ClusterCatalog(
        "act-dr6", act, base_cols=("name", "RADeg", "decDeg", "redshift")
    )
    act.catalog["hms"] = [decimal2hms(ra, ":") for ra in act["ra"]]
    act.catalog["dms"] = [decimal2dms(dec, ":") for dec in act["dec"]]
    if "M500cCal" not in act.colnames:
        if "SZMass_M500cCal" in act.colnames:
            act.catalog["M500cCal"] = act["SZMass_M500cCal"]
        else:
            # for newer catalogs this appears to be the calibrated mass
            act.catalog["M500cCal"] = act["SZMass_M500c"]

    # wave_min = cwl - fwhm / 2
    # wave_max = cwl + fwhm / 2
    # zmin = wave_min / Halpha - 1
    # zmax = wave_max / Halpha - 1
    # ic(zmin, zmax)
    # zmin = np.round(zmin, 3)
    # zmax = np.round(zmax, 3)
    nb964 = targets(act, 0.462, 0.477, "NB964")
    nb1066 = targets(act, 0.621, 0.627, "NB1066")
    nb118 = targets(act, 0.799, 0.817, "NB118")
    ic(act["M500cCal"][nb118])
    ic(nb964.sum(), nb1066.sum(), nb118.sum(), nb964.sum() + nb1066.sum() + nb118.sum())
    # a little beyond one FWHM
    # targets(act, 0.460, 0.479, "NB964")
    # targets(act, 0.619, 0.629, "NB1066")
    # targets(act, 0.797, 0.819, "NB118")

    fig, ax = plt.subplots()
    ax.plot(act["z"], act["M500cCal"], ".", color="0.4", ms=2)
    ax.plot(act["z"][nb964], act["M500cCal"][nb964], "C0o", ms=5, label="NB964")
    ax.plot(act["z"][nb1066], act["M500cCal"][nb1066], "C3o", ms=5, label="NB1066")
    ax.plot(act["z"][nb118], act["M500cCal"][nb118], "C1o", ms=5, label="NB118")
    # January science verification
    sv = np.in1d(act["name"], ["ACT-CL J1038.7+0236", "ACT-CL J1104.5+1110"])
    ax.plot(
        act["z"][sv], act["M500cCal"][sv], "C9x", ms=4, mew=1.5, label="Jan 2024 SV"
    )
    ax.legend(loc="upper right", fontsize=12)
    ax.set(
        xlim=(0.01, 2),
        ylim=(1, 20),
        yscale="log",
        xlabel="Redshift",
        ylabel="$M_{500}$ ($10^{14}\,\mathrm{M_\odot}$)",
    )
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    output = "plots/zmass.pdf"
    savefig(output, fig=fig)

    return


def targets(act, zmin, zmax, sample_name):
    print(f"*** {sample_name} ***")
    ic(np.unique(act["redshiftType"]))
    jz = (act["z"] >= zmin) & (act["z"] <= zmax)
    spec = act["redshiftType"] == "spec"
    ic(jz.sum(), spec.sum(), (jz & spec).sum())
    sample = jz & spec

    coords = SkyCoord(act["ra"], act["dec"], unit="deg")
    fixed_targets = [
        FixedTarget(c, name=name)
        for c, name in zip(coords[sample], act["name"][sample])
    ]
    ctio = Observer.at_site("ctio")
    constraints = [
        AltitudeConstraint(40 * u.deg, 80 * u.deg),
        AirmassConstraint(1.6),
        AtNightConstraint.twilight_astronomical(),
    ]
    # NEWFIRM science verification
    fig, axes = plt.subplots(
        1, 3, sharey=True, figsize=(14, 4), constrained_layout=True
    )
    sv_times = Time(["2024-05-10", "2024-05-17", "2024-05-24"])
    # ic(sv_times.mjd)
    for i, (ax, date) in enumerate(zip(axes, sv_times)):
        times = Time(
            np.linspace(
                date.mjd - (1 * u.hr).to(u.day).value,
                date.mjd + (10 * u.hr).to(u.day).value,
                50,
            ),
            format="mjd",
        )
        good = is_observable(constraints, ctio, fixed_targets, times)
        # ic(good)
        # ic(times.iso)
        plot_airmass(
            [targ for isgood, targ in zip(good, fixed_targets) if isgood],
            ctio,
            times.iso,
            ax=ax,
            brightness_shading=True,
            max_airmass=2,
            use_local_tz=True,
        )
    for ax in axes[1:]:
        ax.set_ylabel("")
        ax.legend(fontsize=12)
    output = f"plots/airmass_sv_{sample_name}.png"
    savefig(output, fig=fig, tight=False)
    # all clusters are observable for at least 3 months during 2024A
    # so let's forget about that for now
    time_range_24a = Time(["2024-02-01 12:00", "2024-07-31 12:00"])
    best_months_24a = months_observable(
        constraints, ctio, fixed_targets, time_range_24a
    )
    time_range_24b = Time(["2024-08-01 12:00", "2025-01-31 12:00"])
    best_months_24b = months_observable(
        constraints, ctio, fixed_targets, time_range_24b
    )
    ic(best_months_24a, best_months_24b)
    # priorities for semesters A and B
    targets = [[], []]
    months = np.zeros((12, 2), dtype=int)
    for i, (obs_a, obs_b) in enumerate(zip(best_months_24a, best_months_24b)):
        cl = act[sample][i]
        sem = np.argmax([len(obs_a), len(obs_b)])
        obs = [obs_a, obs_b][sem]
        # prioritized in this semester if it is hard to observe in the other one
        priority = 2 - int(len([obs_a, obs_b][1 - sem]) < 3)
        block = "B" if sample_name == "NB964" else "A"
        targets[sem].append(
            [block]
            + list(cl["name", "hms", "dms"])
            # + ["-".join([str(i) for i in list(obs)])]
            + [priority]
        )
        for i in range(12):
            if i + 1 in obs_a or i + 1 in obs_b:
                months[i, 0] += 1
                if priority == 1:
                    months[i, 1] += 1
    # this is necessary for the CNTAC csv files
    ic(targets[0], targets[1], len(targets[0]), len(targets[1]))
    ic(months.T, 1 + np.argsort(months[:, 0])[::-1], 1 + np.argsort(months[:, 1])[::-1])
    targets = [
        [[t[i] for t in targ] for i in range(len(targets[0][0]))] for targ in targets
    ]
    targets[0][1] = [name.replace(" ", "_") for name in targets[0][1]]
    targets[1][1] = [name.replace(" ", "_") for name in targets[1][1]]
    # names = ("block", "name", "ra", "dec", "months", "priority")
    names = ("block", "name", "ra", "dec", "priority")
    targets_yr = []
    for targlist, sem in zip(targets, "AB"):
        targets_sem = Table(targlist, names=names)
        output = f"act_Ha_targets_sem{sem}_{sample_name}.csv"
        targets_sem.write(output, format="ascii.csv", overwrite=True)
        print(f"Saved to {output}")
        targets_yr.append(targets_sem)
    targets = vstack(targets_yr)
    targets.write("act_Ha_targets.csv", format="ascii.csv", overwrite=True)
    targets["name", "ra", "dec"].write(
        "act_Ha_mast_search.csv", format="ascii.csv", overwrite=True
    )
    # now do the INT Object Visibility

    scales = (0.262, 0.262, 0.262, 0.262, 0.262, 0.524, 0.524)
    sizes = (2, 3, 5, 10, 20, 30, 60)  # arcmin
    scales = (0.262, 0.262, 0.262)
    sizes = (10, 20, 30)
    for size, scale in zip(sizes, scales):
        url_root = (
            f"https://www.legacysurvey.org/viewer/"
            "jpeg-cutout?ra=<ra>&dec=<dec>&layer=ls-dr10&"
            f"pixscale={scale}&size={int(60*size/scale)}"
        )
        with open(f"decals_download_jpeg__{size}_{scale}.sh", "w") as f:
            for name, ra, dec in act["name", "ra", "dec"][sample]:
                url = url_root.replace("<ra>", f"{ra:.5f}").replace(
                    "<dec>", f"{dec:.5f}"
                )
                path = f"images/{size}_{scale}"
                os.makedirs(path, exist_ok=True)
                filename = os.path.join(
                    path, f"{name}__{size}_{scale}.jpg".replace(" ", "_")
                )
                if os.path.isfile(filename):
                    os.remove(filename)
                print(f'wget -nc -nd -O {filename} "{url}"', file=f)
        url_root = (
            f"https://www.legacysurvey.org/viewer/"
            f"fits-cutout?ra=<ra>&dec=<dec>&layer=ls-dr10&"
            f"pixscale={scale}&size={int(60*size/scale)}&bands=griz"
        )
        with open(f"decals_download_fits__{size}_{scale}.sh", "w") as f:
            for name, ra, dec in act["name", "ra", "dec"][sample]:
                url = url_root.replace("<ra>", f"{ra:.5f}").replace(
                    "<dec>", f"{dec:.5f}"
                )
                path = f"images/fits_{size}_{scale}"
                os.makedirs(path, exist_ok=True)
                filename = os.path.join(
                    path, f"{name}__{size}_{scale}.fits".replace(" ", "_")
                )
                if os.path.isfile(filename):
                    os.remove(filename)
                print(f'wget -nc -nd -O {filename} "{url}"', file=f)

    # RA, DEC color- or size-coding by SZ mass

    return sample


if __name__ == "__main__":
    main()
