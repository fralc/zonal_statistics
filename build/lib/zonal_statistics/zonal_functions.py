import logging
import operator
import numpy as np
from rasterstats import zonal_stats
import geopandas as gpd
from zonal_statistics.utils import get_vector_gdf, retrieve_overlapping, weight_average, check_datafield_aggmethod

logging.basicConfig(level=logging.INFO)


def zonal_raster(base_vector_src, raster_src, agg_method='max', categorical=False, output_src=None):
    """
    Compute statistics values from a raster based on a vector.
    Args:
        base_vector_src: a reference to a vector source. It can be a vector file as read by gpd.read_file (type str),
            a reference to a postgis layer (type dict), or a geodataframe (type gpd.GeoDataFrame).
        raster_src (str): a reference to a raster source
        agg_method: str or list, one or more among these possibilities: 'weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'.
            If it is a list each item correspond to a variable in data_field. If data_field is a list and agg_method is
            a str, then the same type is applied to each variable.
        categorical (bool): whether the raster variable should be considered as categorical (agg_method = 'max_area')
        output_src (str): vector output file source as written by gpd.to_file.

    Returns:
        gpd.GeoDataFrame: a base_vector geodataframe with values from raster_src

    """
    # TODO read raster from postgis
    # TODO check CRSs
    # TODO validate results
    # TODO check whether it make sense to consider more raster at once
    if categorical:
        agg_method = 'max_area'
        logging.info('categorical==True, then agg_method has been set to max_area.')

    if agg_method in ['max_area', 'weight_avg']:
        stats = None
        categorical = True
    else:
        stats = agg_method

    base_gdf, geom_col_base = get_vector_gdf(base_vector_src)
    zs = zonal_stats(base_gdf,
                     raster_src,
                     all_touched=True,
                     stats=stats,
                     categorical=categorical)

    if agg_method == 'max_area':
        raster_values = [max(z.items(), key=operator.itemgetter(1))[0] for z in zs]
    elif agg_method == 'weight_avg':
        raster_values = [weight_average(z) for z in zs]
    elif agg_method in ['mean', 'median', 'sum', 'max', 'min', 'std']:
        raster_values = [z[agg_method] for z in zs]
    else:
        raise ValueError("agg_method must either be 'weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'")
    base_gdf['raster_val'] = raster_values

    if output_src:
        base_gdf.to_file(output_src)

    return base_gdf


def zonal_vector(base_vector_src, data_vector_src, data_field, output_src=None, agg_method='weight_avg'):
    """
    Compute statistics values from a vector based on another vector.
    Args:
        base_vector_src: a reference to a vector source. It can be a vector file as read by gpd.read_file (type str),
            a reference to a postgis layer (type dict), or a geodataframe (type gpd.GeoDataFrame).
        data_vector_src: a reference to a vector source. It can be a vector file as read by gpd.read_file (type str),
            a reference to a postgis layer (type dict), or a geodataframe (type gpd.GeoDataFrame).
        data_field: str or list, columns names of data vector variables to be considered in the analysis
        output_src (str): vector output file source as written by gpd.to_file.
        agg_method: str or list, one or more among these possibilities: 'weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'.
            If it is a list each item correspond to a variable in data_field. If data_field is a list and agg_method is
            a str, then the same type is applied to each variable.

    Returns:
        gpd.GeoDataFrame: a base_vector geodataframe with values from data_vector_src

    """
    # TODO validate results

    data_field, agg_method = check_datafield_aggmethod(data_field, agg_method)

    # read data
    base_gdf, geom_col_base = get_vector_gdf(base_vector_src)
    data_gdf, geom_col_data = get_vector_gdf(data_vector_src)

    # check data_field type: if it is not numerical (categorical) set agg_method to max_area
    for e, df in enumerate(data_field):
        if not np.issubdtype(data_gdf[df].dtype, np.number):
            agg_method[e] = 'max_area'
            logging.info('data_field[{}] is not numeric ({}: {}). agg_method is set to max_area.'.format(e, df, data_gdf[df].dtype))

    # check crs
    if base_gdf.crs != data_gdf.crs:
        logging.info('CRS of layers are different. Data layer is going to be reprojected.')
        data_gdf = data_gdf.to_crs(base_gdf.crs)

    # overlap datasets and retrieve data_field values for base_gdf geometries based on agg_method
    output_gdf = retrieve_overlapping(base_gdf, data_gdf, data_field, base_field=None, agg_methods=agg_method)

    # export if requested
    if output_src:
        output_gdf.to_file(output_src)

    return output_gdf

