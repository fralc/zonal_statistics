import numpy as np
import geopandas as gpd
from sqlalchemy import create_engine


def retrieve_overlapping(base_gdf, data_gdf, data_fields, agg_methods, base_field=None):

    # check whether a base_field is provided; if not, a temporary one is created
    adhoc_base_field = False
    if not base_field:
        base_field = 'this_index'
        base_gdf[base_field] = base_gdf.reset_index().index
        adhoc_base_field = True

    # select only field of interest and geometry
    base_gpd_only_geom = base_gdf[[base_field, base_gdf.geometry.name]]
    data_gpd_only_geom = data_gdf[data_fields + [data_gdf.geometry.name]]

    # intersect polygons and geohashes
    union_gdf = gpd.overlay(data_gpd_only_geom, base_gpd_only_geom, how='union')

    # compute the intersection area
    if any(am in agg_methods for am in ['weight_avg', 'max_area']):
        union_gdf['area'] = union_gdf.area

    data_values_to_base_list = []
    idx_max_area = None
    for df, am in zip(data_fields, agg_methods):
        if am == 'weight_avg':
            # compute products for weighted average
            union_gdf['product'] = union_gdf[df] * union_gdf['area']
            # groupby base field (geometry); retain just products and areas
            data_values_to_base = union_gdf.groupby(base_field)['product', 'area'].agg(np.nansum)
            # compute weighted averages as ratios between products and areas
            data_values_to_base[df] = data_values_to_base['product'] / data_values_to_base['area']
            # data_values_to_base = data_values_to_base[[base_field, df]]
            data_values_to_base = data_values_to_base[df]
        elif am == 'max_area':
            if not idx_max_area:
                # groupby base field (geometry); retain just max area indexes
                idx_max_area = union_gdf.groupby(base_field)['area'].agg(lambda x: x.argmax())
            # get base_field (geometry indexes) and data_field where area is max
            data_values_to_base = union_gdf[[base_field, df]].iloc[idx_max_area.values].set_index(base_field)
        elif am in ['mean', 'median', 'sum', 'max', 'min', 'std']:
            data_values_to_base = union_gdf.groupby(base_field)[df].agg(am)
        else:
            raise ValueError("agg_method must either be 'weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'")
        data_values_to_base_list.append(data_values_to_base)
    # join values to base_gdf
    for dvtb in data_values_to_base_list:
        base_gdf = base_gdf.join(dvtb, on=base_field)

    if adhoc_base_field:
        base_gdf = base_gdf.drop(base_field, axis=1)

    return base_gdf


def get_vector_gdf(vector_src):
    """
    Read a vector data source and returns a geodataframe with data.
    Args:
        vector_src: a reference to a vector source. It can be a vector file as read by gpd.read_file (type str),
            a reference to a postgis layer (type dict), or a geodataframe (type gpd.GeoDataFrame).

    Returns:
        gpd.GeoDataFrame: a geodataframe relative to vector_src
    """
    GEOM_COL = 'geometry'
    if isinstance(vector_src, str):
        polygons = gpd.read_file(vector_src)
    elif isinstance(vector_src, dict):
        polygons, GEOM_COL = gdf_from_postgis(**vector_src)
    elif isinstance(vector_src, gpd.GeoDataFrame):
        polygons = vector_src

    return polygons, GEOM_COL


def gdf_from_postgis(**conn_dict):
    if 'geom_col' not in conn_dict.keys():
        conn_dict['geom_col'] = 'geom'
    conn_str = 'postgresql://{}:{}@{}:{}/{}'.format(
        conn_dict['user'],
        conn_dict['password'],
        conn_dict['host'],
        conn_dict['port'],
        conn_dict['database']
    )
    engine = create_engine(conn_str).raw_connection()
    sql = 'select * from "{}"'.format(conn_dict['table'])
    return gpd.GeoDataFrame.from_postgis(sql, engine, geom_col=conn_dict['geom_col']), conn_dict['geom_col']


def check_datafield_aggmethod(data_field, agg_method):
    if isinstance(data_field, str):
        data_field = [data_field]
    elif not isinstance(data_field, list):
        raise TypeError('data_field must either a str or list')

    if isinstance(agg_method, str):
        agg_method = [agg_method] * len(data_field)
    elif isinstance(agg_method, list):
        if len(agg_method) != len(data_field):
            raise ValueError('len(agg_method) must be equal of len(data_field).')
    else:
        raise TypeError('agg_method must be either str or list')

    for am in agg_method:
        if am not in ['max_area', 'weight_avg', 'max', 'min', 'std']:
            raise ValueError('agg_method is {}. It should be in [max_area, weight_avg, max, min, std]'.format(am))

    return data_field, agg_method


def weight_average(input_dict):
    """
    Compute weighted average from a dict like {value_1: weight_1, ...}
    Args:
        input_dict (dict):

    Returns:

    """
    values = [(v * w, w) for v, w in list(input_dict.items())]
    num_div = [np.sum(v) for v in zip(*values)]
    return (num_div[0] / num_div[1]) if num_div[1] != 0 else None

