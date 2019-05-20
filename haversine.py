from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r


def get_lat2 (dy, lat1): #dy in km
    lat2 = dy/6371. +lat1
    return lat2


def get_lon2 (dx, lat1, lon1): #dy in km
    inner = sin((dx/(2.*6371))/(cos(lat1)))

    lon2 = 2*asin(inner)+lon1
    return lon2

