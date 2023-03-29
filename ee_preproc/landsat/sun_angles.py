import ee
from ..misc.util import *


def create(date, footprint):
    jdp = date.getFraction("year")
    seconds_in_hour = 3600
    hour_gmt = ee.Number(date.getRelative("second", "day")).divide(seconds_in_hour)

    lat_rad = deg_to_rad(ee.Image.pixelLonLat().select("latitude"))
    long_deg = ee.Image.pixelLonLat().select("longitude")

    # Julian day proportion in radians
    jdpr = jdp.multiply(PI()).multiply(2)

    a = ee.List([0.000075, 0.001868, 0.032077, 0.014615, 0.040849])
    mean_solar_time = long_deg.divide(15.0).add(ee.Number(hour_gmt))
    local_solar_diff1 = (
        value(a, 0)
        .add(value(a, 1).multiply(jdpr.cos()))
        .subtract(value(a, 2).multiply(jdpr.sin()))
        .subtract(value(a, 3).multiply(jdpr.multiply(2).cos()))
        .subtract(value(a, 4).multiply(jdpr.multiply(2).sin()))
    )

    local_solar_diff2 = local_solar_diff1.multiply(12 * 60)

    local_solar_diff = local_solar_diff2.divide(PI())
    true_solar_time = mean_solar_time.add(local_solar_diff.divide(60)).subtract(12.0)

    # Hour as an angle
    ah = true_solar_time.multiply(deg_to_rad(ee.Number(MAX_SATELLITE_ZENITH * 2)))
    b = ee.List([0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480])
    delta = (
        value(b, 0)
        .subtract(value(b, 1).multiply(jdpr.cos()))
        .add(value(b, 2).multiply(jdpr.sin()))
        .subtract(value(b, 3).multiply(jdpr.multiply(2).cos()))
        .add(value(b, 4).multiply(jdpr.multiply(2).sin()))
        .subtract(value(b, 5).multiply(jdpr.multiply(3).cos()))
        .add(value(b, 6).multiply(jdpr.multiply(3).sin()))
    )
    cos_sun_zen = (
        lat_rad.sin()
        .multiply(delta.sin())
        .add(lat_rad.cos().multiply(ah.cos()).multiply(delta.cos()))
    )

    sun_zen = cos_sun_zen.acos()

    # sun azimuth from south, turning west
    sin_sun_az_sw = ah.sin().multiply(delta.cos()).divide(sun_zen.sin())
    sin_sun_az_sw = sin_sun_az_sw.clamp(-1.0, 1.0)

    cos_sun_az_sw = (
        lat_rad.cos()
        .multiply(-1)
        .multiply(delta.sin())
        .add(lat_rad.sin().multiply(delta.cos()).multiply(ah.cos()))
    ).divide(sun_zen.sin())
    sun_az_sw = sin_sun_az_sw.asin()

    sun_az_sw = where(cos_sun_az_sw.lte(0), sun_az_sw.multiply(-1).add(PI()), sun_az_sw)
    sun_az_sw = where(
        cos_sun_az_sw.gt(0).And(sin_sun_az_sw.lte(0)), sun_az_sw.add(PI().multiply(2)), sun_az_sw
    )

    sun_az = sun_az_sw.add(PI())
    # Keep within [0, 2pi] range
    sun_az = where(sun_az.gt(PI().multiply(2)), sun_az.subtract(PI().multiply(2)), sun_az)

    footprint_polygon = ee.Geometry.Polygon(footprint)
    sun_az = sun_az.clip(footprint_polygon)
    sun_az = sun_az.rename(["sunAz"])
    sun_zen = sun_zen.clip(footprint_polygon).rename(["sunZen"])

    return sun_az, sun_zen
