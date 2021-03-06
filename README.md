# Zonal Statistics

Tools for extracting zonal statistics from geographic layers.

It provides four functions:
 * `zonal_vector`: extracts data from a vector data source based on another vector source.
 * `zonal_raster`: extracts data from a raster data source based on a vector source.
 * `vector_to_geohash`: extracts data from a vector data source based on a geohash grid.
 * `raster_to_geohash`: extracts data from a raster data source based on a geohash grid.

Input data can be files or layers retrieved from a PostGIS server. 

Geohashes are built by functions based on the data source extent and desired geohash precision.

Each function estimates values from a data source based on a given aggregation method.
Aggregation methods are managed by *agg_method* argument.
If the data source is a vector it is possible to retrieve information about several fields at once providing a list 
of fields to *data_field* argument.

*data_field* and *agg_method* arguments have following features:

* *data_field* can either be a string that identifies a data source field or a list of them;
* *agg_method* can either be a string among possible aggregation method identifiers or a list of them. Aggregation
    method identifiers are ['weight_mean', 'max_area', 'mean', 'median', 'sum', 'max', 'min', 'std'].
    
    'max_area' provides the data field value corresponding to the data source feature with max intersection area.
    
    'weight_avg' provides a weighted average data field value considering intersection areas as weights.
* if *data_field* is a list, *agg_method* must either be a single aggregation method (to be applied to all *data_field*
    items) or a list of corresponding aggregation methods;
* either the data source variable is categorical (not numerical) or it must be considered as such, the aggregation method
    to use is 'max_area'.
    
    `zonal_raster` and `raster_to_geohash` have a *categorical* argument that must be set to True for categorical 
    variables. In this case *agg_method* value is ignored;
    
    `zonal_vector` and `vector_to_geohash` manage categorical variables by means of *agg_method*. 
    If *agg_method* (or one of its items) is not 'max_area' but corresponding data field dtype is not numerical, 
    it is set to 'max_area'.

## Usage examples:

### 1. Get information from a **Raster** layer based on a **Vector** layer:
```
from zonal_statistics import zonal_raster

input_vector_file = 'comuni.shp'
input_raster_file = 'temperature.tif'
output_vector_file = 'avg_temp_on_comuni.shp'

zonal_raster(
    raster_src=input_raster_file, 
    vector_src=input_vector_file, 
    output_src=output_vector_file, 
    agg_method='mean'
)
```   

### 2. Get information from a **Vector** layer based on a **Vector layer**:
```
from zonal_statistics import zonal_vector

input_vector_base_file = 'comuni.shp'
input_vector_data_file = 'temperature_on_zones.shp'
output_vector_file = 'avg_temp_on_comuni.shp'

output_gdf = zonal_vector(
    base_vector_src=input_vector_base_file, 
    data_vector_src=input_vector_data_file, 
    data_field='temperature', 
    output_src=output_vector_file, 
    agg_method='mean'
)
```

### 3. Get information from a **Raster** layer rearranged into **gehashes**:
```
from zonal_statistics import zonal_raster

raster_fn = r"myvar.tif"
geohash_fn = r"myvar_on_geohashes.csv"

output_gdf = zonal_raster(
    raster_src=raster_fn, 
    geohash_precision=4, 
    agg_method='weight_avg',
    categorical=False,
    output_geohash_src=geohash_fn
)
```

### 4. Get information from a **Vector** layer rearranged into **gehashes**:
```
from zonal_statistics import vector_to_geohash

vector_fn = r"myvar_on_zones.shp"
geohash_fn = r"myvar_on_geohashes.csv"

output_gdf = vector_to_geohash(
    data_vector_src=vector_fn, 
    data_field='myvar', 
    geohash_precision=4, 
    agg_method='weight_avg',
    categorical=False,
    output_geohash_src=geohash_fn
)
```
    
### 5. Get information from a **PostGIS Vector** layer based on a **Vector layer**
```
from zonal_statistics import zonal_vector

conn_dict = {
    'user': 'myusername',
    'password': 'mypassword',
    'host': 'myhost.it',
    'port': 5432,
    'database': 'dbname',
    'table': 'temperature_on_zones'
}

input_vector_base_file = 'comuni.shp'
output_vector_file = 'avg_temp_on_comuni.shp'

output_gdf = zonal_vector(
    base_vector_src=input_vector_base_file, 
    data_vector_src=conn_dict, 
    data_field='temperature', 
    output_src=output_vector_file, 
    agg_method='mean'
)
```

### 6. Get information from a **PostGIS Vector** layer rearranged into **gehashes**
```
from zonal_statistics import vector_to_geohash

conn_dict = {
    'user': 'myusername',
    'password': 'mypassword',
    'host': 'myhost.it',
    'port': 5432,
    'database': 'dbname',
    'table': 'myvar_on_zones'
}

input_vector_base_file = 'comuni.shp'
geohash_fn = r"myvar_on_geohashes.csv"

output_gdf = vector_to_geohash(
    data_vector_src=conn_dict, 
    data_field='myvar', 
    geohash_precision=4, 
    agg_method='weight_avg',
    output_geohash_src=geohash_fn
)
```
