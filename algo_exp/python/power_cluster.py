import argparse
from sys import exit, stdout

import geopandas as gpd
import numpy as np
import pandas as pd
from libpysal.weights import Rook

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import ticker

from tqdm import tqdm

# Needed to keep the centroids in the boundary.
from shapely.geometry import Polygon, Point, LinearRing

import os, glob


def ens_dir(f, quiet=False):
    try:
        os.makedirs(f)
    except FileExistsError:
        if not quiet:
            print("File {} exists.  Cleaning it.".format(f))

    for f in glob.glob(f + "*"):
        os.remove(f)


def wavg(grp):
    return grp._get_numeric_data().multiply(grp["a"], axis=0).sum() / grp["a"].sum()


class cluster:
    def __init__(
        self,  geo_path: str, name: str, crs, seats: int, seed: int, attrs: list = [], quiet=False
    ):

        geo = gpd.read_file(geo_path).to_crs(crs)
        self.name = name.upper()
        self.seed = seed
        self.nloops = 0
        self.attrs = ['x', 'y', 'pct_black', 'pct_hispanic', 'pct_asian', 'pct_white']
        print(self.attrs)

        self.crs = crs
        self.seats = seats

        boundary = geo.dissolve().geometry[0]
        self.area = boundary.area
        self.poly = boundary.simplify(5)

        if type(self.poly) is Polygon:
            self.rings = [LinearRing(self.poly.exterior.coords)]
        else:
            self.rings = [LinearRing(p.exterior.coords) for p in self.poly]

        self.initialize_tract_frame(geo)
        self.initialize_districts()

        self.voronoi_classify()

        self.folder = f"power/{self.name}/{self.seed}/"
        ens_dir(self.folder, quiet)

        if not quiet:
            self.plot_map("init.pdf".format(name, seed))

        self.print_freq = 500

    def initialize_districts(self):

        # A random sample of tracts, as seeds.
        np.random.seed(self.seed)
        cdf = self.gdf.iloc[np.random.choice(self.ntracts, self.seats, replace=False)][self.attrs
        ].reset_index(drop=True)

        # Initalize to just the largest radius...
        distances = (
            (
                cdf[self.attrs].to_numpy()
                - cdf[self.attrs].to_numpy()[:, np.newaxis, :]
            )
            ** 2
        ).sum(axis=2)
        cdf["r2"] = np.amin(np.ma.masked_equal(distances, 0), axis=1).data

        self.cdf = cdf

    def initialize_tract_frame(self, geo):

        geo['wrook'] = pd.Series(
            Rook.from_dataframe(geo).neighbors
        )

        centroids = geo.geometry.centroid
        geo['x'] = centroids.apply(lambda x: x.x)
        geo['y'] = centroids.apply(lambda x: x.y)
        geo['points'] = [Point(xy) for xy in zip(geo.x, geo.y)]

        geo['a'] = geo.geometry.area
        self.target = geo["pop"].sum() / self.seats
        self.ntracts = geo.shape[0]

        self.gdf = geo

    def scale(self, rate=0.01, clip=0.01):

        groups = self.gdf.groupby("C")

        dist_pop = groups["pop"].sum()
        pop = pd.Series(np.zeros(self.seats))
        pop.iloc[dist_pop.index] = dist_pop

        dfactor = (self.target - pop) / self.target

        # Signed squared difference.
        clipped_factor = (
            np.clip(np.abs(dfactor) * dfactor, -clip, +clip) * groups["a"].sum()
        )
        clipped_factor.fillna(clip * self.area / self.seats, inplace=True)

        self.cdf["S"] = clipped_factor
        self.cdf["r2"] += clipped_factor

    def voronoi_classify(self):

        # scipy.spatial.distance is equally fast...
        d2 = (
            (
                self.gdf[self.attrs].to_numpy()
                - self.cdf[self.attrs].to_numpy()[:, np.newaxis, :]
            )
            ** 2
        ).sum(axis=2)

        self.gdf["C"] = np.argmin(d2 - self.cdf["r2"].to_numpy()[:, np.newaxis], axis=0)

    def plot_map(self, filename=None, annotate=True, figsize=4):

        """
        Somewhat fancy function, for nicely monitoring
        the evolution of district centroids in cdf and
        tract points in gdf.
        """

        dis = self.gdf.dissolve("C", aggfunc="sum")
        dis["frac"] = dis["pop"] / self.target

        bounds = self.poly.bounds
        xr = bounds[2] - bounds[0]
        yr = bounds[3] - bounds[1]

        col, alpha, trunc = "coolwarm", 0.7, ""

        if dis["frac"].max() > 2:
            norm = Normalize(vmin=0, vmax=2)
            trunc = " (Truncated)"
        elif dis["frac"].max() - 1 < 1e-5:
            norm = Normalize(vmin=0.9999, vmax=1.0001)
        else:  # regardless, keep it centered
            larger = max(1 - dis["frac"].min(), dis["frac"].max() - 1)
            norm = Normalize(vmin=1 - larger, vmax=1 + larger)

        cmap = plt.cm.ScalarMappable(norm=norm, cmap=col)

        ax = dis.plot(
            color="white",
            edgecolor="black",
            figsize=(figsize * np.sqrt(xr / yr), figsize * np.sqrt(yr / xr)),
        )
        for xi, row in dis.iterrows():
            dis[dis.index == xi].plot(
                ax=ax, alpha=alpha, facecolor=cmap.to_rgba(row["frac"])
            )

        if annotate:
            ax.scatter(
                self.cdf.x, self.cdf.y, c="k", s=5, edgecolor="w", marker="o", zorder=3
            )

            for ri, row in self.cdf.iterrows():
                ax.annotate(str(ri), (row.x, row.y), zorder=4)
                ax.add_artist(
                    plt.Circle(
                        (row.x, row.y), np.sqrt(max(row.r2, 1)), color="k", fill=False
                    )
                )

        fig = ax.get_figure()
        cax = fig.add_axes([0.16, 0.13, 0.70, 0.015 * np.sqrt(xr / yr)])
        sm = plt.cm.ScalarMappable(cmap=col, norm=norm)
        sm._A = []  # gross

        cb = fig.colorbar(
            sm,
            cax=cax,
            alpha=alpha,
            label="Population / Target" + trunc,
            orientation="horizontal",
        )
        cb.locator = ticker.MaxNLocator(nbins=5)
        cb.formatter.set_useOffset(False)
        cb.update_ticks()

        ax.set_xlim([bounds[0] - 0.1 * xr, bounds[2] + 0.1 * xr])
        ax.set_ylim([bounds[1] - 0.1 * yr, bounds[3] + 0.1 * yr])

        ax.set_axis_off()

        if not filename:
            return ax

        ax.figure.savefig(self.folder + filename, bbox_inches="tight", pad_inches=0.05)
        plt.close("all")

    def closest_point_on_boundary(self, x, y):

        point = Point(x, y)

        mini = np.inf
        for ring in self.rings:

            d = ring.project(point)
            if d < mini:
                mini = d
                pt_in_poly = ring.interpolate(d)

        # return pt_in_poly
        return list(pt_in_poly.coords)[0]

    def go_where_its_hot(self, centers):

        # get the directions and distance squared to other cells.
        vectors = (
            centers[["x", "y"]].to_numpy()
            - centers[["x", "y"]].to_numpy()[:, np.newaxis]
        )
        distance = np.sqrt(np.square(vectors).sum(axis=2))

        distance = np.ma.masked_equal(distance, 0)
        unit_vectors = vectors / distance[:, :, np.newaxis]

        # populations and strenths of the pulls.
        pop_diff = (
            centers["pop"].to_numpy() - centers["pop"].to_numpy()[:, np.newaxis]
        ) / self.target
        pop_diff[pop_diff < 0] = 0

        print("\n:::: POP :::", (pop_diff / distance)[:, :, np.newaxis])
        move = (pop_diff / distance).T[:, :, np.newaxis]

        step = np.sum(vectors * move, axis=0).data
        print("\n:::: STEP\n", step)
        return step

    def loop_voronoi(self, nloops, stop=2.5e-3):

        if self.seats == 1:
            return  # no point!

        nloops = int(nloops)

        ## for i in tqdm(range(self.nloops, self.nloops + nloops), ncols = 80,
        ##               desc = "{} ({})".format(self.name, self.seed)): # status bar!

        for i in range(self.nloops, self.nloops + nloops):

            # calculate a scaling from the populations
            self.voronoi_classify()

            self.cdf["pop"] = self.gdf.groupby("C")["pop"].sum()

            if np.abs(1 - self.cdf["pop"] / self.target).max() < stop:
                break

            self.scale(rate=0.01, clip=0.01)

            # watch out for the nearest (district) neighbors
            # you get "zapped" for bullying smaller districts...
            distance2 = (
                (
                    self.cdf[["x", "y"]].to_numpy()
                    - self.cdf[["x", "y"]].to_numpy()[  # broadcast to seats × seats
                        :, np.newaxis, :
                    ]
                )
                ** 2
            ).sum(axis=2)
            distance2 = np.ma.masked_equal(distance2, 0)  # robust...
            idx = np.nonzero(distance2 <= self.cdf.r2.to_numpy()[np.newaxis, :])
            idx = np.array([np.append(idx[0], idx[1]), np.append(idx[1], idx[0])])

            # The zapping is away from the "adversary"
            # with each one stepping in proportion to its population.
            self.cdf["pop"].fillna(1, inplace=True)
            pop_ratio = self.cdf["pop"].to_numpy()[:, np.newaxis] / (
                self.cdf["pop"].to_numpy()[np.newaxis, :] + self.cdf["pop"].to_numpy()[:, np.newaxis]
            )

            D = pop_ratio[idx[0], idx[1]][:, np.newaxis] * (
                self.cdf.loc[idx[1], self.attrs].values
                - self.cdf.loc[idx[0], self.attrs]
            )

            D = D.groupby(D.index).sum()
            self.cdf.loc[D.index, self.attrs] -= 0.05 * D

            # Make sure that it hasn't left the state!
            for idx, r in self.cdf.loc[D.index].iterrows():
                if not self.poly.contains(Point(r.x, r.y)):
                    self.cdf.loc[idx, ["x", "y"]] = self.closest_point_on_boundary(
                        r.x, r.y
                    )

            # at any rate, the power radius mustn't exceed
            # the distance to nearest neighbor's centroid.
            # recalculate distances post-zapping.
            distance2 = (
                (
                    self.cdf[["x", "y"]].to_numpy()
                    - self.cdf[["x", "y"]].to_numpy()[  # broadcast to seats × seats
                        :, np.newaxis, :
                    ]
                )
                ** 2
            ).sum(axis=2)
            distance2 = np.ma.masked_equal(
                distance2, 0
            )  # notes on masking method, below.
            closests = np.amin(distance2, axis=1).data
            self.cdf["r2"] = np.minimum(self.cdf.r2, closests)

            # "regularize" -- can't go below negative avg area.
            # cdf.loc[(cdf["pop"] < 1.1 * self.target) & (cdf["r2"] < 1), "r2"] = 1
            # cdf["r2"] = np.maximum(cdf.r2, - state_area/seats)

            ############
            # Also move *slowly* towards the geographic center.
            self.voronoi_classify()

            # calculate lines towards the new points
            baryctr = self.gdf.groupby("C").apply(wavg)[["x", "y"]].dropna()
            old_pts = self.cdf.loc[baryctr.index, ["x", "y"]]

            # how far to step.
            dist_pop = self.gdf.groupby("C")["pop"].sum().loc[baryctr.index]
            mv_rate = np.power(10.0, -0.5 - dist_pop / self.target)
            mv_rate = mv_rate.values.reshape(1, mv_rate.size)

            self.cdf.loc[baryctr.index, ["x", "y"]] = (
                1 - mv_rate.T
            ) * old_pts + mv_rate.T * baryctr

            ### annd.... move towards "hot" districts.
            # if i % 25 == 0: go_where_its_hot(self.cdf)
            # self.cdf[["x", "y"]] -= go_where_its_hot(self.cdf)

            # Monitoring and plotting.
            self.voronoi_classify()

            if i % self.print_freq == 0:
                self.plot_map("{:05d}.pdf".format(i))

        self.plot_map("final.pdf")
        self.plot_map("final_clear.pdf", annotate=False)

        self.nloops += nloops


def main(file, name, crs, seats, seed, nloop, stop):

    weights = {
        'pct_black' : 2,
        'pct_hispanic' : 1.5,
        'pct_asian': 1.5,
        'pct_white': 1
    }

    races = ['pct_black',
    'pct_hispanic',
    'pct_asian',
    'pct_white']

    clust = cluster(geo_path=file, name=name, crs=crs, seats=seats, seed=seed)
    clust.loop_voronoi(nloop, stop=stop)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file", help="The path to the file with geometry and population data."
    )
    parser.add_argument("--name", help="The name to give this run.")
    parser.add_argument("--seats", type=int, help="Number of seats to create.")
    parser.add_argument("--crs", help="Projection to enforce on the geometry.")
    parser.add_argument("--seed", default=1, type=int)
    parser.add_argument("--nloop", default=int(4e3), type=int)
    parser.add_argument("--stop", default=2.5e-3, type=float)
    args = parser.parse_args()

    main(args.file, args.name, args.crs, args.seats, args.seed, args.nloop, args.stop)
