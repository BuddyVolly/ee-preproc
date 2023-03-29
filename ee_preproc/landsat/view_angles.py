import ee
from ..misc.util import *

MAX_DISTANCE = 1000000


def create(footprint):
    return azimuth(footprint), zenith(footprint)


def azimuth(footprint):
    upper_center = (
        line_from_coords(footprint, UPPER_LEFT, UPPER_RIGHT).centroid().coordinates()
    )
    lower_center = (
        line_from_coords(footprint, LOWER_LEFT, LOWER_RIGHT).centroid().coordinates()
    )
    slope = ((y(lower_center)).subtract(y(upper_center))).divide(
        (x(lower_center)).subtract(x(upper_center))
    )
    slope_perp = ee.Number(-1).divide(slope)
    azimuth_left = ee.Image(PI().divide(2).subtract(slope_perp.atan()))
    return azimuth_left.rename(["viewAz"])


def zenith(footprint):
    left_line = line_from_coords(footprint, UPPER_LEFT, LOWER_LEFT)
    right_line = line_from_coords(footprint, UPPER_RIGHT, LOWER_RIGHT)
    left_distance = ee.FeatureCollection(left_line).distance(MAX_DISTANCE)
    right_distance = ee.FeatureCollection(right_line).distance(MAX_DISTANCE)
    view_zenith = (
        right_distance.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2))
        .divide(right_distance.add(left_distance))
        .subtract(ee.Number(MAX_SATELLITE_ZENITH))
        .clip(ee.Geometry.Polygon(footprint))
        .rename(["viewZen"])
    )
    return deg_to_rad(view_zenith)
