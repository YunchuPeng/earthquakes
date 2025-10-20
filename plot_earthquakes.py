from datetime import date

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import os
os.makedirs("plots", exist_ok=True)

import requests
import json
from collections import defaultdict
import statistics as stats



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
    # The time is given in a strange-looking but commonly-used format.
    # To understand it, we can look at the documentation of the source data:
    # https://earthquake.usgs.gov/data/comcat/index.php#time
    # Fortunately, Python provides a way of interpreting this timestamp:
    # (Question for discussion: Why do we divide by 1000?)
    year = date.fromtimestamp(timestamp/1000).year
    return year


def get_magnitude(earthquake):
    """Retrive the magnitude of an earthquake item."""
    return earthquake['properties'].get('mag', None)


# This is function you may want to create to break down the computations,
# although it is not necessary. You may also change it to something different.
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



# Get the data we will work with
quakes = get_data()['features']

# Plot the results - this is not perfect since the x axis is shown as real
# numbers rather than integers, which is what we would prefer!
plot_number_per_year(quakes)
plt.clf()  # This clears the figure, so that we don't overlay the two plots
plot_average_magnitude_per_year(quakes)