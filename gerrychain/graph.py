import enum
import json
import warnings
from collections import Counter

import geopandas as gp
import networkx
import pandas as pd
import pysal
from networkx.readwrite import json_graph
from shapely.ops import cascaded_union

from gerrychain.utm import from_latlon


def utm_of_point(point):
    return from_latlon(point.y, point.x)[2]


def identify_utm_zone(df):
    wgs_df = df.to_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    utm_counts = Counter(utm_of_point(point) for point in wgs_df["geometry"].centroid)
    # most_common returns a list of tuples, and we want the 0,0th entry
    most_common = utm_counts.most_common(1)[0][0]
    return most_common


def reprojected(df):
    """Returns a copy of `df`, projected into the coordinate reference system of a suitable
        `Universal Transverse Mercator`_ zone.
    :param df: :class:`geopandas.GeoDataFrame`
    :rtype: :class:`geopandas.GeoDataFrame`

    .. _`Universal Transverse Mercator`: https://en.wikipedia.org/wiki/UTM_coordinate_system
    """
    utm = identify_utm_zone(df)
    return df.to_crs(
        f"+proj=utm +zone={utm} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    )


class Adjacency(enum.Enum):
    """There are two concepts of "adjacency" that one can apply when constructing
    an adjacency graph.

    - **Rook adjacency**: Two polygons are adjacent if they share one or more *edges*.
    - **Queen adjacency**: Two polygons are adjacent if they share one or more *vertices*.

    All Rook edges are Queen edges, but not the other way around. Many congressional
    districts are only Queen-contiguous (i.e., they consist of two polygons connected
    at a single corner).
    """
    Rook = "rook"
    Queen = "queen"


adjacency_functions = {
    Adjacency.Rook: pysal.weights.Rook.from_dataframe,
    Adjacency.Queen: pysal.weights.Queen.from_dataframe,
}


def get_neighbors(df, adjacency):
    return adjacency_functions[adjacency](df, geom_col=df.geometry.name).neighbors


class Graph(networkx.Graph):
    @classmethod
    def from_json(cls, json_file):
        """Load a graph from a JSON file in the NetworkX json_graph format.
        :json_file: Path to JSON file.
        :returns: Graph
        """
        with open(json_file) as f:
            data = json.load(f)
        g = json_graph.adjacency_graph(data)
        return cls(g)

    @classmethod
    def from_file(cls, filename, cols_to_add=None, reproject=True):
        df = gp.read_file(filename)
        return cls.from_geodataframe(df, cols_to_add, reproject)

    @classmethod
    def from_geodataframe(cls, dataframe, adj=Adjacency.Rook, reproject=True):
        """Creates the adjacency :class:`Graph` of geometries described by `dataframe`.
        The areas of the polygons are included as node attributes (with key `area`).
        The shared perimeter of neighboring polygons are included as edge attributes
        (with key `shared_perim`).
        Nodes corresponding to polygons on the boundary of the union of all the geometries
        (e.g., the state, if your dataframe describes VTDs) have a `boundary_node` attribute
        (set to `True`) and a `boundary_perim` attribute with the length of this "exterior"
        boundary.

        By default, areas and lengths are computed in a UTM projection suitable for the
        geometries. This prevents the bizarro area and perimeter values that show up when
        you accidentally do computations in Longitude-Latitude coordinates. If the user
        specifies `reproject=False`, then the areas and lengths will be computed in the
        GeoDataFrame's current coordinate reference system. This option is for users who
        have a preferred CRS they would like to use.

        :param dataframe: :class:`geopandas.GeoDataFrame`
        :param adj: (optional) The adjacency type to use. Default is `Adjacency.Rook`.
            Other options are `Adjacency.Queen`, "rook" or "queen".
        :return: The adjacency graph of the geometries from `dataframe`.
        :rtype: :class:`Graph`
        """
        # Project the dataframe to an appropriate UTM projection unless
        # explicitly told not to.
        if reproject:
            df = reprojected(dataframe)
        else:
            df = dataframe

        # Generate rook neighbor lists from dataframe.
        neighbors = get_neighbors(df, adj)

        # Add shared ("interior") perimeters to edges between nodes
        adjacency_dict = neighbors_with_shared_perimeters(neighbors, df)
        graph = cls.from_dict_of_dicts(adjacency_dict)

        # Add "exterior" perimeters to the boundary nodes
        add_boundary_perimeters(graph, neighbors, df)

        # Add area data to the nodes
        areas = df["geometry"].area.to_dict()
        networkx.set_node_attributes(graph, name="area", values=areas)

        return graph

    def add_data(self, df, columns=None):
        """Add columns of a DataFrame to a graph as node attributes using
        by matching the DataFrame's index to node ids.

        :param df: Dataframe containing given columns.
        :param columns: (optional) List of dataframe column names to add.
        :return: Nothing.
        """

        if columns is None:
            columns = df.columns

        check_dataframe(df[columns])

        column_dictionaries = df.to_dict("index")
        networkx.set_node_attributes(self, column_dictionaries)

    def assignment(self, node_attribute_key):
        """Create an assignment dictionary using an attribute of the nodes
        of the graph. For example, if you created your graph from Census data
        and each node has a `CD` attribute that gives the congressional district
        the node belongs to, then `graph.assignment("CD")` would return the
        desired assignment of nodes to CDs.

        :param graph: NetworkX graph.
        :param node_attribute_key: Attribute available on all nodes.
        :return: Dictionary of {node_id: attribute} pairs.
        """
        return networkx.get_node_attributes(self, node_attribute_key)

    def join(self, dataframe, columns=None, left_index=None, right_index=None):
        """Add data from a dataframe to the graph, matching nodes to rows when
        the node's `left_index` attribute equals the row's `right_index` value.

        :param dataframe: DataFrame.
        :columns: (optional) The columns whose data you wish to add to the graph.
            If not provided, all columns are added.
        :left_index: (optional) The node attribute used to match nodes to rows.
            If not provided, node IDs are used.
        :right_index: (optional) The DataFrame column name to use to match rows
            to nodes. If not provided, the DataFrame's index is used.
        """
        if right_index is not None:
            df = dataframe.set_index(right_index)
        else:
            df = dataframe

        if columns is not None:
            df = df[columns]

        check_dataframe(df)

        column_dictionaries = df.to_dict()

        if left_index is not None:
            ids_to_index = networkx.get_node_attributes(self, left_index)
        else:
            # When the left_index is node ID, the matching is just
            # a redundant {node: node} dictionary
            ids_to_index = dict(zip(self.nodes, self.nodes))

        node_attributes = {
            node_id: {
                column: values[index] for column, values in column_dictionaries.items()
            }
            for node_id, index in ids_to_index.items()
        }

        networkx.set_node_attributes(self, node_attributes)


def add_boundary_perimeters(graph, neighbors, df):
    """Add shared perimeter between nodes and the total geometry boundary.

    :graph: NetworkX graph.
    :neighbors: Adjacency information generated from pysal.
    :df: Geodataframe containing geometry information.
    :returns: The updated graph.

    """
    # creates one shape of the entire state to compare outer boundaries against
    boundary = cascaded_union(df.geometry).boundary

    intersections = df.intersection(boundary)
    is_boundary = intersections.apply(bool)

    # Add boundary node information to the graph.
    intersection_df = gp.GeoDataFrame(intersections)
    intersection_df["boundary_node"] = is_boundary

    # List-indexing here to get the correct dictionary format for NetworkX.
    attr_dict = intersection_df[["boundary_node"]].to_dict("index")
    networkx.set_node_attributes(graph, attr_dict)

    # For the boundary nodes, set the boundary perimeter.
    boundary_perims = intersections[is_boundary].length
    boundary_perims = gp.GeoDataFrame(boundary_perims)
    boundary_perims.columns = ["boundary_perim"]

    attribute_dict = boundary_perims.to_dict("index")
    networkx.set_node_attributes(graph, attribute_dict)


def neighbors_with_shared_perimeters(neighbors, df):
    """Construct a graph with shared perimeter between neighbors on the edges.

    :neighbors: Adjacency information generated from pysal.
    :df: Geodataframe containing geometry information.
    :returns: A dict of dicts of the following form::

        { node: { neighbor: { shared_perim: <value> }}}

    """
    adjacencies = {}

    geom = df.geometry
    for shape in neighbors:
        shared_perim = geom[neighbors[shape]].intersection(geom[shape]).length
        shared_perim.name = "shared_perim"
        adjacencies[shape] = pd.DataFrame(shared_perim).to_dict("index")

    return adjacencies


def check_dataframe(df):
    for column in df.columns:
        if sum(df[column]) > 0:
            warnings.warn("NA values found in column {}!".format(column))
