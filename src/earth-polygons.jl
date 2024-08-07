module EarthPolygons
using UlamMethod


"""
    North_Atlantic_clipped

Vertices defining a polygon contained in the rectangle with lower left 
corner (-100, -9) and upper right corner (15, 39). This polygon clips away the 
land, leaving only the ocean in the rectangle.
"""
North_Atlantic_clipped = [
    -76.3 39; 
    -5.89 39; 
    -5.61 34.3; 
    -9.63 29.3; 
    -14.2 24.9; 
    -15.9 21.5; 
    -15.9 14.; 
    -8.15 5.09; 
    4.13 7.24; 
    5.96 4.95; 
    9.14 4.82; 
    10.9 2.79;
    9. -1.11;
    14.2 -9; 
    -35.5 -9; 
    -35.7 -5.8; 
    -40. -3.25; 
    -50.5 0.371; 
    -52. 3.89; 
    -56.7 5.65; 
    -63.3 10.; 
    -72.4 11.6; 
    -76.4 8.1; 
    -79.1 9.56; 
    -81.4 9.15; 
    -83.9 10.8; 
    -84.2 15.9; 
    -88.4 16.1; 
    -86.4 21.1; 
    -89.9 22.5; 
    -91.8 18.9; 
    -95.8 19.; 
    -98. 22.1; 
    -97.2 28.1; 
    -93.4 29.8; 
    -84. 29.9; 
    -82. 25.; 
    -79.3 25.8; 
    -81.4 31.7; 
    -76.5 34.9
] |> permutedims |> Boundary

"""
    North_Atlantic_box

Vertices defining a rectangle with lower left corner (-100, -9) and upper right corner (15, 39).
"""
North_Atlantic_box = [
    -100 -9;
    -100 39;
    15 39;
    15 -9
] |> permutedims |> Boundary

"""
    GoG_big

Vertices defining the entirety of the Gulf of Guinea.
"""
GoG_big = [
    -1.38 6.48;
    4.34 7.88; 
    10.1 4.52; 
    11.4 -0.871; 
    1.01 -1.49
] |> permutedims |> Boundary

"""
    GoG_small

Vertices defining the interior of the Gulf of Guinea closest to the shore.
"""
GoG_small = [
    6.96 4.82;
    6.82 0.0591;
    11.2 0.0827; 
    11.1 4.25
] |> permutedims |> Boundary

end # module