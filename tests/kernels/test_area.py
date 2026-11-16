"""The area kernel and ``toolkit.area``: SPARTN v2.0.3 8.10, Table 6.24, the
polynomial of 8.11 and 8.12 (SF041, SF054), and the grid interpolation of
8.11 and 8.12 with the weights of 8.13.

Synthetic areas only. The figures of 8.10 (area layout, node order,
centroid; inclusion examples) and 8.13 (rectangular cell) are read from the
PDF of the specification; the tests below restate what they show.
"""

import math

import pytest

from ssrcore._kernels import area as k
from ssrcore.toolkit import area
from ssrcore.toolkit.area import Area

# The layout of the 8.10 figure: 5 latitude nodes, 6 longitude nodes.
FIGURE = Area(50.0, 10.0, 5, 6, 0.5, 0.75)


def node_position(a: Area, index: int) -> tuple[float, float]:
    """Position of a node from the 8.10 figure's order, written independently of the kernel."""
    row, col = divmod(index, a.longitude_nodes)
    return a.reference_latitude - row * a.latitude_spacing, a.reference_longitude + col * a.longitude_spacing


# ---- extent and inclusion (8.10) ---------------------------------------------

def test_extent_runs_south_and_east_from_the_reference_point():
    # SF032 is the northern-most latitude, SF033 the western-most longitude;
    # SF034/SF035 count nodes from them (Table 6.2).
    assert area.extent(FIGURE) == (50.0, 48.0, 10.0, 13.75)


def test_inside_edges_and_outside():
    assert area.contains(FIGURE, 49.0, 12.0)
    for lat, lon in ((50.0, 10.0), (48.0, 13.75), (50.0, 13.75), (48.0, 10.0), (49.0, 10.0)):
        assert area.contains(FIGURE, lat, lon), (lat, lon)
    for lat, lon in ((50.01, 12.0), (47.99, 12.0), (49.0, 9.99), (49.0, 13.76), (-49.0, 12.0), (49.0, -170.0)):
        assert not area.contains(FIGURE, lat, lon), (lat, lon)


def test_second_figure_of_8_10_adjacent_and_separated_areas():
    # Left of the figure: two areas side by side share an edge; positions in
    # either are served, a position outside both is not. Right: two areas one
    # above the other with a gap; a position in the gap is in neither.
    west = Area(50.0, 0.0, 4, 6, 1.0, 1.0)   # 47..50 N, 0..5 E
    east = Area(49.0, 5.0, 3, 5, 1.0, 1.0)   # 47..49 N, 5..9 E
    assert area.contains(west, 48.5, 2.0) and not area.contains(east, 48.5, 2.0)
    assert area.contains(east, 48.0, 7.0) and not area.contains(west, 48.0, 7.0)
    assert area.contains(west, 48.0, 5.0) and area.contains(east, 48.0, 5.0)  # the shared edge
    assert not area.contains(west, 49.5, 6.0) and not area.contains(east, 49.5, 6.0)
    assert not area.contains(west, 46.5, 2.0) and not area.contains(east, 46.5, 2.0)

    upper = Area(40.0, 0.0, 4, 6, 1.0, 1.0)  # 37..40 N
    lower = Area(36.0, 0.0, 5, 6, 1.0, 1.0)  # 32..36 N
    assert area.contains(upper, 37.5, 2.0) and area.contains(lower, 33.0, 4.0)
    assert not area.contains(upper, 36.5, 2.0) and not area.contains(lower, 36.5, 2.0)


def test_not_above_85_degrees():
    # 8.10: "should not be used above +/- 85 degrees latitude".
    north = Area(88.0, 0.0, 5, 3, 1.0, 1.0)   # 84..88 N
    assert area.contains(north, 84.5, 1.0)
    assert area.contains(north, 85.0, 1.0)
    assert not area.contains(north, 85.5, 1.0)
    assert not area.contains(north, 88.0, 1.0)
    south = Area(-84.0, 0.0, 5, 3, 1.0, 1.0)  # -88..-84
    assert area.contains(south, -85.0, 1.0)
    assert not area.contains(south, -85.01, 1.0)
    # The grid itself is still defined there; only the use is refused.
    assert area.grid(north, 86.0, 1.0)[1] is not None


def test_longitude_wraps_at_the_antimeridian():
    a = Area(10.0, 178.0, 3, 5, 1.0, 1.0)  # 178 E .. 182 E = 178 W
    assert area.extent(a) == (10.0, 8.0, 178.0, 182.0)
    assert area.contains(a, 9.0, 179.5)
    assert area.contains(a, 9.0, -179.5)
    assert area.contains(a, 9.0, -178.0)
    assert area.contains(a, 9.0, 181.0)
    assert not area.contains(a, 9.0, -177.9)
    assert not area.contains(a, 9.0, 177.9)
    # A reference longitude written as west longitude is the same area.
    b = a._replace(reference_longitude=-182.0)
    assert area.contains(b, 9.0, -179.5) and area.contains(b, 9.0, 179.0)


@pytest.mark.parametrize("bad", [
    Area(50.0, 10.0, 0, 3, 1.0, 1.0),
    Area(50.0, 10.0, 3, 9, 1.0, 1.0),
    Area(50.0, 10.0, 3, 3, 0.0, 1.0),
    Area(50.0, 10.0, 3, 3, 1.0, -1.0),
    Area(90.5, 10.0, 3, 3, 1.0, 1.0),
    Area(-89.0, 10.0, 3, 3, 1.0, 1.0),  # extends below -90
])
def test_definitions_that_are_not_areas(bad):
    for f in (area.extent, area.centroid):
        with pytest.raises(ValueError):
            f(bad)
    with pytest.raises(ValueError):
        area.contains(bad, 0.0, 0.0)
    assert k.extent(*bad)[0] == k.INVALID


def test_node_counts_of_table_6_2_are_accepted():
    for n in range(1, 9):
        assert k.extent(0.0, 0.0, n, n, 3.2, 0.1)[0] == k.OK


# ---- centroid and offsets (8.10 figure, 8.11, 8.12; Q8) -----------------------

def test_centroid_is_the_centre_of_the_node_grid_of_the_figure():
    # In the 8.10 figure (5 x 6 nodes) the polynomial expansion point lies on
    # the middle (third) row, half-way between the third and fourth columns.
    lat, lon = area.centroid(FIGURE)
    assert lat == pytest.approx(50.0 - 2 * 0.5)
    assert lon == pytest.approx(10.0 + 2.5 * 0.75)
    # The mean of the four corner nodes, an independent route.
    corners = [node_position(FIGURE, i) for i in (0, 5, 24, 29)]
    assert lat == pytest.approx(sum(p[0] for p in corners) / 4)
    assert lon == pytest.approx(sum(p[1] for p in corners) / 4)


def test_centroid_across_the_antimeridian_is_wrapped():
    lat, lon = area.centroid(Area(10.0, 178.0, 3, 6, 1.0, 1.0))  # 178 .. 183 E
    assert (lat, lon) == pytest.approx((9.0, -179.5))


def test_offsets_from_the_centroid():
    dlat, dlon = area.offsets(FIGURE, 49.5, 11.0)
    assert dlat == pytest.approx(49.5 - 49.0)
    assert dlon == pytest.approx(11.0 - 11.875)
    a = Area(10.0, 178.0, 3, 4, 1.0, 1.0)  # centroid 9 N, 179.5 E
    assert area.offsets(a, 8.5, -179.5) == pytest.approx((-0.5, 1.0))
    assert area.offsets(a, 8.5, 178.5) == pytest.approx((-0.5, -1.0))


# ---- grid nodes and weights (8.11, 8.12, 8.13; Q9) ----------------------------

def test_node_order_of_the_figure():
    # The 8.10 figure orders nodes from the reference point west to east,
    # then row after row southwards. At a node, that node has weight 1.
    for index in range(FIGURE.latitude_nodes * FIGURE.longitude_nodes):
        lat, lon = node_position(FIGURE, index)
        nodes, weights = area.grid(FIGURE, lat, lon)
        got = {n: 0.0 for n in nodes}
        for n, w in zip(nodes, weights):
            got[n] += w
        assert got[index] == pytest.approx(1.0, abs=1e-12), index
        assert sum(w for n, w in got.items() if n != index) == pytest.approx(0.0, abs=1e-12)


def test_cell_corners_follow_8_13():
    # 8.13 rectangular: W1 (lambda2, phi2) = x y, W2 (lambda1, phi2) = (1-x) y,
    # W3 (lambda1, phi1) = (1-x)(1-y), W4 (lambda2, phi1) = x (1-y).
    nodes, weights = area.grid(FIGURE, 49.4, 11.2)
    # Row 1 (49.5 N) and row 2 (49.0 N); columns 1 (10.75 E) and 2 (11.5 E).
    assert nodes == (1 * 6 + 2, 1 * 6 + 1, 2 * 6 + 1, 2 * 6 + 2)
    x = (11.2 - 10.75) / 0.75
    y = (49.4 - 49.0) / 0.5
    assert weights == pytest.approx((x * y, (1 - x) * y, (1 - x) * (1 - y), x * (1 - y)))


def test_weights_sum_to_one_and_reproduce_bilinear_fields():
    # A bilinear field in latitude and longitude is reproduced exactly by
    # bilinear interpolation of its node values: an independent check of
    # nodes and weights together.
    def field(lat, lon):
        return 3.0 - 2.0 * lat + 0.5 * lon + 0.25 * lat * lon

    a = Area(-20.0, 175.0, 6, 8, 0.7, 1.1)  # crosses the antimeridian
    north, south, west, east = area.extent(a)
    n = 37
    for i in range(n):
        for j in range(n):
            lat = south + (north - south) * i / (n - 1)
            lon_east = west + (east - west) * j / (n - 1)
            lon = lon_east - 360.0 if lon_east >= 180.0 else lon_east
            nodes, weights = area.grid(a, lat, lon)
            assert math.fsum(weights) == pytest.approx(1.0, abs=1e-14)
            assert all(-1e-12 <= w <= 1.0 + 1e-12 for w in weights)
            value = sum(w * field(*node_position(a, nd)) for nd, w in zip(nodes, weights))
            assert value == pytest.approx(field(lat, lon_east), abs=1e-9), (lat, lon)


def test_one_node_in_a_direction():
    row = Area(45.0, 0.0, 1, 4, 1.0, 1.0)  # a line of latitude
    nodes, weights = area.grid(row, 45.0, 1.25)
    assert sum(weights) == pytest.approx(1.0)
    assert nodes[:2] == (2, 1) and weights[0] == pytest.approx(0.25) and weights[1] == pytest.approx(0.75)
    column = Area(45.0, 0.0, 3, 1, 1.0, 1.0)
    nodes, weights = area.grid(column, 44.5, 0.0)
    assert sum(weights) == pytest.approx(1.0)
    assert set(nodes) == {0, 1}
    single = Area(45.0, 0.0, 1, 1, 1.0, 1.0)
    assert area.grid(single, 45.0, 0.0) == ((0, 0, 0, 0), (0.0, 1.0, 0.0, 0.0))


def test_grid_outside_raises():
    with pytest.raises(ValueError):
        area.grid(FIGURE, 47.9, 12.0)
    assert k.grid(47.9, 12.0, *FIGURE)[0] == k.OUTSIDE


# ---- polynomial (8.11, 8.12; SF041, SF054) ---------------------------------------

C = (12.4, -1.6, 0.8, 0.05)  # ionosphere C00 TECU, C01 TECU/deg, C10 TECU/deg, C11 TECU/deg^2
T = (0.123, -0.004, 0.0025, 0.0003)  # troposphere T00 m, T01, T10 m/deg, T11 m/deg^2


@pytest.mark.parametrize("coefficients", [C, T], ids=["SF054", "SF041"])
@pytest.mark.parametrize("dlat,dlon", [(0.0, 0.0), (0.7, -1.3), (1.5, -2.25), (-3.0, 4.0)])
def test_polynomial_types_0_to_2(coefficients, dlat, dlon):
    c00, c01, c10, c11 = coefficients
    assert area.polynomial(0, coefficients[:1], dlat, dlon) == c00
    assert area.polynomial(1, coefficients[:3], dlat, dlon) == pytest.approx(
        c00 + c01 * dlat + c10 * dlon, abs=1e-14)
    assert area.polynomial(2, coefficients, dlat, dlon) == pytest.approx(
        c00 + c01 * dlat + c10 * dlon + c11 * dlat * dlon, abs=1e-14)


def test_polynomial_at_the_centroid_is_c00():
    for t in (0, 1, 2):
        assert area.polynomial(t, C, 0.0, 0.0) == 12.4


def test_c01_goes_with_latitude_and_c10_with_longitude():
    # SF041, SF054: "C00 + C01 (delta phi) + C10 (delta lambda)".
    assert area.polynomial(1, (0.0, 1.0, 0.0), 2.0, 5.0) == 2.0
    assert area.polynomial(1, (0.0, 0.0, 1.0), 2.0, 5.0) == 5.0


def test_unused_coefficients_are_not_read():
    nan = float("nan")
    s, value = k.polynomial(0, (3.0, nan, nan, nan), 1.0, 1.0)
    assert s == k.OK and value == 3.0
    s, value = k.polynomial(1, (3.0, 1.0, 1.0, nan), 1.0, 1.0)
    assert s == k.OK and value == 5.0


@pytest.mark.parametrize("t", [3, 4, 5, 6, 7, -1])
def test_undefined_equation_types(t):
    # SF041 and SF054 values 3 to 7 are "TBD".
    assert k.polynomial(t, C, 0.0, 0.0)[0] == k.EQUATION_TYPE
    with pytest.raises(ValueError):
        area.polynomial(t, C, 0.0, 0.0)


def test_too_few_coefficients():
    with pytest.raises(ValueError):
        area.polynomial(2, C[:3], 0.0, 0.0)


def test_polynomial_at_a_user_position_of_an_area():
    # The offsets are taken from the area centroid (8.11, 8.12): a type-2
    # polynomial written about the centroid, evaluated at a user position.
    a = Area(52.0, 4.0, 4, 5, 1.0, 1.5)        # centroid 50.5 N, 7.0 E
    dlat, dlon = area.offsets(a, 51.0, 5.5)     # +0.5, -1.5
    expected = 12.4 - 1.6 * 0.5 + 0.8 * -1.5 + 0.05 * 0.5 * -1.5
    assert area.polynomial(2, C, dlat, dlon) == pytest.approx(expected)


# ---- grid residual (8.11, 8.12) ---------------------------------------------------

def test_interpolate_is_the_weighted_sum():
    w = (0.1, 0.2, 0.3, 0.4)
    r = (1.0, -2.0, 0.5, 0.28)
    assert area.interpolate(w, r) == pytest.approx(0.1 - 0.4 + 0.15 + 0.112)
    r = (0.012, -0.008, 0.020, 0.004)
    assert area.interpolate((0.25, 0.25, 0.25, 0.25), r) == pytest.approx(sum(r) / 4.0, abs=1e-16)
    for node in range(4):
        w = [0.0] * 4
        w[node] = 1.0
        assert area.interpolate(w, r) == r[node]


def test_interpolate_at_a_user_position_from_area_weights():
    # Residuals defined on every node of a 3 x 4 grid as a bilinear field are
    # reproduced at any point of the area.
    a = Area(45.0, 10.0, 3, 4, 0.5, 0.5)
    residual = [0.04 * i - 0.2 for i in range(12)]   # linear in node index: bilinear in position
    nodes, weights = area.grid(a, 44.3, 11.1)
    got = area.interpolate(weights, [residual[n] for n in nodes])
    row, col = (45.0 - 44.3) / 0.5, (11.1 - 10.0) / 0.5
    assert got == pytest.approx(0.04 * (row * 4 + col) - 0.2)


def test_invalid_node_gives_no_grid_residual():
    # A7: any of the four nodes at the invalid value (SF052, SF053; SF064 to
    # SF067) leaves the polynomial alone, also a node of zero weight.
    assert area.interpolate((0.25, 0.25, 0.25, 0.25), (1.0, None, 1.0, 1.0)) is None
    assert area.interpolate((1.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, None)) is None
    assert k.interpolate((0.25,) * 4, (1.0,) * 4, (1, 1, 0, 1))[0] == k.INVALID_NODE
