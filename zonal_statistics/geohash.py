import logging
import numpy as np
import geopandas as gpd
import geohash
import rasterio as rio
from polygon_geohasher.polygon_geohasher import geohash_to_polygon
from zonal_statistics.utils import get_vector_gdf, retrieve_overlapping, check_datafield_aggmethod
from zonal_statistics.zonal_functions import zonal_raster

logging.basicConfig(level=logging.INFO)


def vector_to_geohash(data_vector_src, data_field,
                      geohash_precision=4, agg_method='weight_avg', output_geohash_src=None):
    """
    Compute statistics values from a vector based on a geohash grid.
    Args:
        data_vector_src: a reference to a vector source. It can be a vector file as read by gpd.read_file (type str),
            a reference to a postgis layer (type dict), or a geodataframe (type gpd.GeoDataFrame).
        data_field: str or list, columns names of data vector variables to be considered in the analysis
        geohash_precision (int): geohash precision in the range 0-8
        agg_method: str or list, one or more among these possibilities: 'weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'.
            If it is a list each item correspond to a variable in data_field. If data_field is a list and agg_method is
            a str, then the same type is applied to each variable.
        output_geohash_src (str): vector output file source as written by gpd.to_file.

    Returns:
        gpd.GeoDataFrame: a geohash geodataframe with values from data_vector_src

    """
    # TODO validate results

    data_field, agg_method = check_datafield_aggmethod(data_field, agg_method)

    assert geohash_precision in range(0, 9)

    polygons, geom_col = get_vector_gdf(data_vector_src)

    # check data_field type: if it is not numerical (categorical) set agg_method to max_area
    for e, df in enumerate(data_field):
        if not np.issubdtype(polygons[df].dtype, np.number):
            agg_method[e] = 'max_area'
            logging.info('data_field[{}] is not numeric ({}: {}). '
                         'agg_method is set to max_area.'.format(e, df, polygons[df].dtype))

    # transform crs
    polygons = polygons.to_crs({'init': 'epsg:4326'})

    # select only geometry and field of interest
    polygons = polygons[data_field + [geom_col]].reset_index()

    # get bounds for geohashes identification
    bounds = polygons.bounds
    extent = [bounds.minx.min(), bounds.miny.min(), bounds.maxx.max(), bounds.maxy.max()]

    # estimate a list of geohashes based on polygons extent, transform into geometries and store in a geodataframe
    geohash_list = _compute_geohash_tiles(extent, geohash_precision)
    geohash_grid = geohashes_to_gpd(geohash_list)

    # overlap datasets and retrieve data_field values for base_gdf geometries based on agg_method
    output_gdf = retrieve_overlapping(geohash_grid, polygons, data_field, base_field='geohash', agg_methods=agg_method)

    # export if requested
    if output_geohash_src:
        if output_geohash_src[-3:] == 'csv':
            output_gdf[data_field + ['geohash']].to_csv(output_geohash_src, header=True)
        else:
            output_gdf.to_file(output_geohash_src)

    return output_gdf


def raster_to_geohash(raster_src, geohash_precision=4, agg_method='weight_avg',
                      categorical=False, output_geohash_src=None):
    """
    Compute statistics values from a raster based on a geohash grid.
    Args:
        raster_src (str): a reference to a raster source
        geohash_precision (int): geohash precision in the range 0-8
        agg_method: str or list, one or more among these possibilities: 'weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'.
            If it is a list each item correspond to a variable in data_field. If data_field is a list and agg_method is
            a str, then the same type is applied to each variable.
        categorical (bool): whether the raster variable should be considered as categorical (agg_method = 'max_area')
        output_geohash_src (str): vector output file source as written by gpd.to_file.

    Returns:
        gpd.GeoDataFrame: a geohash geodataframe with values from raster_src

    """

    # get raster bounds and related geohash grid
    with rio.open(raster_src, 'r') as src:
        if src.crs.to_epsg() != 4326:
            raise ValueError('input raster must have epsg code 4326')
        extent = list(src.bounds)
    geohash_list = _compute_geohash_tiles(extent, geohash_precision)
    geohash_grid = geohashes_to_gpd(geohash_list)

    # execute zonal_raster
    output_gdf = zonal_raster(
        base_vector_src=geohash_grid,
        raster_src=raster_src,
        agg_method=agg_method,
        categorical=categorical,
        output_src=None
    )

    # export if requested
    if output_geohash_src:
        if output_geohash_src[-3:] == 'csv':
            output_gdf[['raster_val', 'geohash']].to_csv(output_geohash_src, header=True)
        else:
            output_gdf.to_file(output_geohash_src)

    return output_gdf


def geohashes_to_gpd(geohashes):
    GEOM_COLUMN = 'geometry'
    geohash_grid = gpd.GeoDataFrame.from_records(
        zip(
            geohashes,
            map(geohash_to_polygon, geohashes)
        ),
        columns=['geohash', GEOM_COLUMN]
    )
    geohash_grid.set_geometry(GEOM_COLUMN, inplace=True)
    geohash_grid.crs = {'init': 'epsg:4326'}
    return geohash_grid


def _is_geohash_in_bounding_box(current_geohash, bbox_coordinates):
    """Checks if the box of a geohash is inside the bounding box

    :param current_geohash: a geohash
    :param bbox_coordinates: bounding box coordinates
    :return: true if the the center of the geohash is in the bounding box
    """

    coordinates = geohash.decode(current_geohash)
    geohash_in_bounding_box = (bbox_coordinates[1] < coordinates[0] < bbox_coordinates[3]) and (
            bbox_coordinates[0] < coordinates[1] < bbox_coordinates[2])
    return geohash_in_bounding_box


def _compute_geohash_tiles(bbox_coordinates, geohash_precision):
    """
    Computes all geohash tile in the given bounding box
    (modified from https://blog.tafkas.net/2018/09/28/creating-a-grid-based-on-geohashes/)

    Args:
        bbox_coordinates: the bounding box coordinates of the geohashes (minx, miny, maxx, maxy)
        geohash_precision (int): geohash precision (level) of output geohashes

    Returns:
        list: geohashes

    """

    checked_geohashes = set()
    geohash_stack = set()
    geohashes = []
    # get center of bounding box, assuming the earth is flat ;)
    center_latitude = (bbox_coordinates[1] + bbox_coordinates[3]) / 2
    center_longitude = (bbox_coordinates[0] + bbox_coordinates[2]) / 2

    center_geohash = geohash.encode(center_latitude, center_longitude, precision=geohash_precision)
    geohashes.append(center_geohash)
    geohash_stack.add(center_geohash)
    checked_geohashes.add(center_geohash)
    while len(geohash_stack) > 0:
        current_geohash = geohash_stack.pop()
        neighbors = geohash.neighbors(current_geohash)
        for neighbor in neighbors:
            if neighbor not in checked_geohashes and _is_geohash_in_bounding_box(neighbor, bbox_coordinates):
                geohashes.append(neighbor)
                geohash_stack.add(neighbor)
                checked_geohashes.add(neighbor)
    return geohashes
