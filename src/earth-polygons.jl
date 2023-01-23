"""
Utilities for handing regions of Earth.

Written by Gage Bonner September 2022
"""

###################################################################################################
###################################################################################################
###################################################################################################

"""
This is a dictionary of all relevant polygons, defined by a list of vertices (lon, lat) and edges which connect the vertices

north_atlantic =========> the entire North Atlantic between North and South America and Africa
africa_coast ===========> a single box off the coast of West Africa
sargasso_sea ===========> a single box in the Sargasso Sea
gulf_of_mexico_full ====> the entire Gulf Of Mexico, extending out roughly to Cuba
gulf_of_mexico_houston => a single box near the coast of Houston
nirvana_bottom =========> the bottom nirvana, roughly -9 degrees
nirvana_top ============> the top nirvana, roughly +39 degrees
"""

earth_polygons = Dict(
    "north_atlantic_verts" => [
        -82.5 39 ; -6.06 39 ; -5.9 24.9 ;
        15.3 21.9 ; 15.8 -9 ; -63.7 -9 ; 
        -75.67 5.47 ; -78.24 9.06 ; -80.2 9.65 ; 
        -81.8 9.14 ; -84.2 11.56 ; -84.47 14.12 ; 
        -90.1 15.88 ; -97.67 18.2 ; -100.2 26.08 ; 
        -95.5 31.95 ],
    "north_atlantic_edges" => [
        1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 
        6 7 ; 7 8 ; 8 9 ; 9 10 ; 10 11 ;
        11 12 ; 12 13 ; 13 14 ; 14 15 ; 15 16;
        16 1],

    "africa_coast_verts" => [-18.2 16.8 ; -18.2 17.2 ; -17.8 17.2 ; -17.8 16.8],
    "africa_coast_edges" => [1 2 ; 2 3 ; 3 4 ; 4 1],

    "sargasso_sea_verts" => [-55 27.8 ; -55 28 ; -54.8 28 ; -54.8 27.8],
    "sargasso_sea_edges" => [1 2 ; 2 3 ; 3 4 ; 4 1],

    "gulf_of_mexico_full_verts" => [
        -100 17; -100 32 ; -83 32 ; 
        -80 24.3 ; -88.8 19.7 ; -91.2 17],
    "gulf_of_mexico_full_edges" => [1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 6 1],
    
    "gulf_of_mexico_houston_verts" => [-95 27 ; -95 27.2 ; -94.8 27.2 ; -94.8 27],
    "gulf_of_mexico_houston_edges" => [1 2 ; 2 3 ; 3 4 ; 4 1],

    "nirvana_bottom_verts" => [-37.0 -12.0 ; -37 -9 ; 16 -9 ; 16 -12],
    "nirvana_bottom_edges" => [1 2 ; 2 3 ; 3 4 ; 4 1],

    "nirvana_top_verts" => [-75.0 42 ; -6 42 ; -6 39 ; -75 39],
    "nirvana_top_edges" => [1 2 ; 2 3 ; 3 4 ; 4 1],

    # direct from west coast of Africa to gulf of mexico
    "path_fast" => [-21.56 17.32;
    -23.88 17.48; 
    -26.41 17.48; 
    -28.41 17.48; 
    -30.97 17.44;
    -33.28 17.44; 
    -35.53 17.39; 
    -37.86 17.39; 
    -40.71 17.39; 
    -42.78 17.43; 
    -45.03 17.36; 
    -46.97 17.36;
    -49.5 17.39; 
    -52.52 17.35; 
    -54.2 17.29; 
    -56.59 17.25; 
    -58.82 17.4; 
    -59.17 15.17;
    -61.49 15.22; 
    -63.88 15.31; 
    -66 15.15; 
    -68.21 15.08;
    -70.09 15.04; 
    -72.97 14.89;
    -75.31 14.82; 
    -77.41 14.79; 
    -79.79 14.89; 
    -79.96 16.95; 
    -82.37 16.95; 
    -82.66 19.39;
    -84.48 19.38;
    -84.66 21.35; 
    -85.23 23.95; 
    -87.22 24.23; 
    -88.95 24.34],

    # from west africa to gulf of guinae back around to GoM
    "path_slow" => [-19.65 14.72;
    -19.68 13.19;
    -19.61 10.74;
    -16.92  9.807;
    -16.58 8.319;
    -15.14 8.033; 
    -14.8 6.272;
    -12.72 5.654; 
    -11.98 3.77;
    -9.928 3.532; 
    -8.008 3.362; 
    -6.125 3.284; 
    -3.294 2.941; 
    -3.078 0.5494; 
    -3.025 -1.847; 
    -4.762 -2.; 
    -4.965 -4.25;
    -6.961 -3.91;
    -9.258 -3.775;
    -12.48 -3.765; 
    -14.81 -3.816; 
    -17.06 -3.816; 
    -19.11 -3.857; 
    -21.67 -3.898; 
    -23.62 -3.893; 
    -25.95 -3.857;
    -28.92 -3.857;
    -30.79 -3.799; 
    -33.02 -3.751; 
    -35.12 -3.715; 
    -35.59 -1.793;
    -37.26 -1.553; 
    -40. -1.35; 
    -42.04 -1.362;
    -42.67 0.6899; 
    -45.08 0.847; 
    -47.13 0.9352;
    -47.38 2.484; 
    -49.08 2.834; 
    -49.58 4.686; 
    -49.73 7.544; 
    -51.98 7.809;
    -54.59 8.012; 
    -57.05  8.157; 
    -57.09 10.19; 
    -58.89 10.19; 
    -59.1 12.01; 
    -61.09 12.28; 
    -63.27 12.48; 
    -64.05 14.73; 
    -65.51 14.8; 
    -68.51 14.89; 
    -70.72 14.93; 
    -73.09 14.97; 
    -74.56 14.97; 
    -77.27 14.94;
    -79.2 14.9; 
    -79.94 16.9; 
    -82.72 17.35; 
    -82.95 19.47; 
    -84.44 19.61; 
    -84.69 21.19; 
    -84.76 23.33; 
    -86.52 23.89; 
    -88.96 24.21]
)



"""
Takes two regions and returns a dictionary defining a single region which is 
their union. Suitable for use with inpoly2() of PolygonInbounds.
"""

function combine_regions(verts1, edges1, verts2, edges2)
    region_verts = [verts1 ; verts2]
    region_edges = [edges1 ; edges2 .+ length(edges1[:,1])]

    return Dict(
        "verts" => region_verts,
        "edges" => region_edges
    )
end