from datetime import date
import os
import json
from collections import defaultdict
import statistics as stats

import requests
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import animation
from matplotlib import gridspec

os.makedirs("plots", exist_ok=True)


def get_data():
    """Retrieve the data we will be working with."""
    response = requests.get(
        "http://earthquake.usgs.gov/fdsnws/event/1/query.geojson",
        params={
            'starttime': "2000-01-01",
            "maxlatitude": "58.723",
            "minlatitude": "50.008",
            "maxlongitude": "1.67",
            "minlongitude": "-9.756",
            "minmagnitude": "1",
            "endtime": "2018-10-11",
            "orderby": "time-asc"}
    )
   
    text = response.text
    return json.loads(text)  


def get_year(earthquake):
    """Extract the year in which an earthquake happened."""
    timestamp = earthquake['properties']['time']
    year = date.fromtimestamp(timestamp/1000).year
    return year


def get_magnitude(earthquake):
    """Retrive the magnitude of an earthquake item."""
    return earthquake['properties'].get('mag', None)



def get_magnitudes_per_year(earthquakes):
    """Retrieve the magnitudes of all the earthquakes in a given year.
    
    Returns a dictionary with years as keys, and lists of magnitudes as values.
    """
    mags_by_year = defaultdict(list)
    for eq in earthquakes:
        mag = get_magnitude(eq)
        if mag is None:
            continue
        mags_by_year[get_year(eq)].append(float(mag))
    return mags_by_year


def plot_average_magnitude_per_year(earthquakes):
    mags_by_year = get_magnitudes_per_year(earthquakes)
    if not mags_by_year:
        print("No magnitudes available to plot.")
        return

    years = sorted(mags_by_year.keys())
    avgs = [stats.fmean(mags_by_year[y]) for y in years]

    plt.figure()
    plt.plot(years, avgs, marker='o', linewidth=1.8)
    plt.xlabel("Year")
    plt.ylabel("Average magnitude")
    plt.title("Average earthquake magnitude per year")
    plt.xticks(years, rotation=45)
    plt.grid(True, linestyle='--', alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/image1.png", dpi=150)
    print("Saved: plots/image1.png")
    plt.close()  




def plot_number_per_year(earthquakes):
    counts = defaultdict(int)
    for eq in earthquakes:
        counts[get_year(eq)] += 1

    if not counts:
        print("No earthquakes to plot.")
        return

    years = sorted(counts.keys())
    nums = [counts[y] for y in years]

    plt.figure()
    plt.bar(years, nums)
    plt.xlabel("Year")
    plt.ylabel("Number of earthquakes")
    plt.title("Earthquake count per year")
    plt.xticks(years, rotation=45)
    plt.tight_layout()
    plt.savefig("plots/image2.png", dpi=150)
    print("Saved: plots/image2.png")
    plt.close() 


def extract_lon_lat_mag(earthquakes, year=None):
    if year is not None:
        earthquakes = [eq for eq in earthquakes if get_year(eq) == year]

    lons, lats, mags = [], [], []
    for eq in earthquakes:
        mag = get_magnitude(eq)
        if mag is None:
            continue
        lon, lat = eq["geometry"]["coordinates"][:2]
        if lon is None or lat is None:
            continue
        lons.append(float(lon))
        lats.append(float(lat))
        mags.append(float(mag))

    return np.asarray(lons), np.asarray(lats), np.asarray(mags)






def plot_map_with_marginals(earthquakes, year=None, bin_deg=1.0, outfile="plots/map_marginals.png"):
    
    if year is not None:
        earthquakes = [eq for eq in earthquakes if get_year(eq) == year]

    lons, lats, mags = [], [], []
    for eq in earthquakes:
        mag = get_magnitude(eq)
        if mag is None: 
            continue
        lon, lat = eq["geometry"]["coordinates"][:2]
        if lon is None or lat is None:
            continue
        lons.append(float(lon)); lats.append(float(lat)); mags.append(float(mag))

    if not lons:
        print(f"No data to plot{f' for year {year}' if year else ''}.")
        return

    lons = np.asarray(lons); lats = np.asarray(lats); mags = np.asarray(mags)

    
    pad = 0.5
    lon_min, lon_max = lons.min() - pad, lons.max() + pad
    lat_min, lat_max = lats.min() - pad, lats.max() + pad
    lon_bins = np.arange(np.floor(lon_min), np.ceil(lon_max) + bin_deg, bin_deg)
    lat_bins = np.arange(np.floor(lat_min), np.ceil(lat_max) + bin_deg, bin_deg)

    
    fig = plt.figure(figsize=(9, 10))
    gs = gridspec.GridSpec(4, 4, figure=fig, wspace=0.06, hspace=0.06)
    ax_histx = fig.add_subplot(gs[0, 0:3])            
    ax_map   = fig.add_subplot(gs[1:4, 0:3], projection=ccrs.PlateCarree())  
    ax_histy = fig.add_subplot(gs[1:4, 3])            

    
    ax_map.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax_map.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.4)
    ax_map.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.2)
    ax_map.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax_map.add_feature(cfeature.BORDERS, linestyle=':')
    ax_map.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    
    size_scale = 18.0
    sc = ax_map.scatter(
        lons, lats,
        s=mags * size_scale,
        c=mags, cmap='viridis', alpha=0.65, edgecolor='k', linewidth=0.3,
        transform=ccrs.PlateCarree()
    )
    cbar = plt.colorbar(sc, ax=ax_map, orientation='vertical', shrink=0.85)
    cbar.set_label("Magnitude")

    
    ax_histx.hist(lons, bins=lon_bins)
    ax_histx.set_xlim(lon_min, lon_max)
    ax_histx.set_ylabel("Count")
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histx.grid(True, linestyle=':', alpha=0.25)
    ax_histx.set_title(f"Earthquake locations{f' — {year}' if year else ''}")

    
    ax_histy.hist(lats, bins=lat_bins, orientation='horizontal', color='tab:orange')
    ax_histy.set_ylim(lat_min, lat_max)
    ax_histy.set_xlabel("Count")
    ax_histy.tick_params(axis='y', labelleft=False)
    ax_histy.grid(True, linestyle=':', alpha=0.25)


    fig.suptitle("Earthquakes — locations and marginal histograms", y=0.995)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    print(f"Saved: {outfile}")
    plt.close(fig)

def animate_earthquakes_with_marginals(
    earthquakes, by='year', fps=2, bin_deg=1.0,
    outfile="plots/earthquakes_years_marginals.gif"
):


    if by == 'year':
        frames = sorted({ get_year(eq) for eq in earthquakes })
        label_fmt = "{}"
        def data_for_frame(frame_label):
            return extract_lon_lat_mag(earthquakes, year=int(frame_label))
    elif by == 'day':
        from datetime import datetime
        def daykey(eq):
            ts = eq['properties']['time'] / 1000.0
            return datetime.utcfromtimestamp(ts).date().isoformat()
        frames = sorted({ daykey(eq) for eq in earthquakes })
        label_fmt = "{}"
        def data_for_frame(frame_label):
            sel = []
            for eq in earthquakes:
                ts = eq['properties']['time'] / 1000.0
                if datetime.utcfromtimestamp(ts).date().isoformat() == frame_label:
                    sel.append(eq)
            return extract_lon_lat_mag(sel, year=None)
    else:
        raise ValueError("by must be 'year' or 'day'")


    all_lons, all_lats, _ = extract_lon_lat_mag(earthquakes, year=None)
    pad = 0.5
    lon_min, lon_max = all_lons.min() - pad, all_lons.max() + pad
    lat_min, lat_max = all_lats.min() - pad, all_lats.max() + pad
    lon_bins = np.arange(np.floor(lon_min), np.ceil(lon_max) + bin_deg, bin_deg)
    lat_bins = np.arange(np.floor(lat_min), np.ceil(lat_max) + bin_deg, bin_deg)


    max_lon_count = 0
    max_lat_count = 0
    for frame_label in frames:
        lons, lats, _ = data_for_frame(frame_label)
        lon_hist, _ = np.histogram(lons, bins=lon_bins)
        lat_hist, _ = np.histogram(lats, bins=lat_bins)
        if lon_hist.max() > max_lon_count:
            max_lon_count = lon_hist.max()
        if lat_hist.max() > max_lat_count:
            max_lat_count = lat_hist.max()

   
    fig = plt.figure(figsize=(9, 10))
    gs = gridspec.GridSpec(4, 4, figure=fig, wspace=0.06, hspace=0.06)
    ax_histx = fig.add_subplot(gs[0, 0:3])
    ax_map   = fig.add_subplot(gs[1:4, 0:3], projection=ccrs.PlateCarree())
    ax_histy = fig.add_subplot(gs[1:4, 3])

    ax_map.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax_map.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.4)
    ax_map.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.2)
    ax_map.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax_map.add_feature(cfeature.BORDERS, linestyle=':')
    ax_map.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    size_scale = 18.0
    sc = ax_map.scatter([], [], s=[], c=[], cmap='viridis', alpha=0.65,
                        edgecolor='k', linewidth=0.3, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(sc, ax=ax_map, orientation='vertical', shrink=0.85)
    cbar.set_label("Magnitude")

    def style_histx():
        ax_histx.set_xlim(lon_min, lon_max)
        ax_histx.set_ylim(0, max_lon_count * 1.1)   
        ax_histx.set_ylabel("Count")
        ax_histx.tick_params(axis='x', labelbottom=False)
        ax_histx.grid(True, linestyle=':', alpha=0.25)

    def style_histy():
        ax_histy.set_ylim(lat_min, lat_max)
        ax_histy.set_xlim(0, max_lat_count * 1.1)   
        ax_histy.set_xlabel("Count")
        ax_histy.tick_params(axis='y', labelleft=False)
        ax_histy.grid(True, linestyle=':', alpha=0.25)

    def init():
        sc.set_offsets(np.empty((0, 2)))
        sc.set_sizes(np.array([]))
        sc.set_array(np.array([]))
        ax_histx.cla(); ax_histy.cla()
        style_histx(); style_histy()
        ax_histx.set_title("Earthquake locations — animation")
        return (sc,)

    def update(frame_label):
        lons, lats, mags = data_for_frame(frame_label)

        if lons.size == 0:
            sc.set_offsets(np.empty((0, 2)))
            sc.set_sizes(np.array([]))
            sc.set_array(np.array([]))
        else:
            sc.set_offsets(np.column_stack([lons, lats]))
            sc.set_sizes(mags * size_scale)
            sc.set_array(mags)

        
        ax_histx.cla()
        ax_histx.hist(lons, bins=lon_bins)
        style_histx()

        ax_histy.cla()
        ax_histy.hist(lats, bins=lat_bins, orientation='horizontal', color='tab:orange')
        style_histy()

        ax_map.set_title(f"Earthquake locations — {label_fmt.format(frame_label)}")
        return (sc,)

    anim = animation.FuncAnimation(
        fig, update, frames=frames, init_func=init,
        interval=int(1000 / fps), blit=False
    )

    try:
        anim.save(outfile, writer='pillow', dpi=200)
    except Exception:
        mp4 = outfile.rsplit('.', 1)[0] + ".mp4"
        anim.save(mp4, writer='ffmpeg', dpi=200)
        outfile = mp4
    print(f"Saved animation: {outfile}")
    plt.close(fig)






# Get the data we will work with
quakes = get_data()['features']

# Plot the results - this is not perfect since the x axis is shown as real
# numbers rather than integers, which is what we would prefer!
plot_number_per_year(quakes)
plt.clf()  # This clears the figure, so that we don't overlay the two plots
plot_average_magnitude_per_year(quakes)






plot_map_with_marginals(quakes, year=None, bin_deg=1.0, outfile="plots/map_marginals_all.png")
plot_map_with_marginals(quakes, year=2008, bin_deg=1.0, outfile="plots/map_marginals_2008.png")




animate_earthquakes_with_marginals(
    quakes, by='year', fps=2, bin_deg=1.0,
    outfile="plots/earthquakes_years_marginals.gif"
)
