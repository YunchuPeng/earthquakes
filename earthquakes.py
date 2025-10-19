# The Python standard library includes some functionality for communicating
# over the Internet.
# However, we will use a more powerful and simpler library called requests.
# This is external library that you may need to install first.
import requests
import json



def get_data():
    # With requests, we can ask the web service for the data.
    # Can you understand the parameters we are passing here?
    response = requests.get(
        "http://earthquake.usgs.gov/fdsnws/event/1/query.geojson",
        params={
            "starttime": "2000-01-01",
            "maxlatitude": "58.723",
            "minlatitude": "50.008",
            "maxlongitude": "1.67",
            "minlongitude": "-9.756",
            "minmagnitude": "1",
            "endtime": "2018-10-11",
            "orderby": "time-asc",
        },
        timeout=30,
    )
    response.raise_for_status()

    # The response we get back is an object with several fields.
    # The actual contents we care about are in its text field:
    # 方式一：用 text + json.loads
    text = response.text
    # To understand the structure of this text, you may want to save it
    # to a file and open it in VS Code or a browser.
    # See the README file for more information.
    data = json.loads(text)

    # We need to interpret the text to get values that we can work with.
    # What format is the text in? How can we load the values?
    return data

def count_earthquakes(data):
    """Get the total number of earthquakes in the response."""
    # USGS GeoJSON 顶层有 'features' 列表；metadata.count 也给了数量，但以实际 features 为准更稳妥
    return len(data.get("features", []))


def get_magnitude(earthquake):
    """Retrive the magnitude of an earthquake item."""
    # 震级在 properties.mag；可能为 None，调用方需要处理
    return earthquake.get("properties", {}).get("mag")


def get_location(earthquake):
    """Retrieve the latitude and longitude of an earthquake item."""
    # There are three coordinates, but we don't care about the third (altitude)
    coords = earthquake.get("geometry", {}).get("coordinates", [])
    if len(coords) >= 2:
        lon, lat = coords[0], coords[1]
        return (lat, lon)
    return (None, None)


def get_maximum(data):
    """Get the magnitude and location of the strongest earthquake in the data."""
    features = data.get("features", [])
    max_mag = None
    max_loc = (None, None)

    for feat in features:
        mag = get_magnitude(feat)
        if mag is None:  # 跳过缺失震级的条目
            continue
        if (max_mag is None) or (mag > max_mag):
            max_mag = mag
            max_loc = get_location(feat)

    return max_mag, max_loc


# With all the above functions defined, we can now call them and get the result
data = get_data()
print(f"Loaded {count_earthquakes(data)}")
max_magnitude, max_location = get_maximum(data)
print(f"The strongest earthquake was at {max_location} with magnitude {max_magnitude}")