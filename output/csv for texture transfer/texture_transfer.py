import pyvista as pv
import numpy as np
from numpy import genfromtxt
import imageio

# MODERN THEME SETTING
pv.global_theme.background = 'white'

# Load mapping prediction (Using updated local absolute path)
#pred = genfromtxt(r'shrec16\1-17 full m=15 n350.csv', delimiter=',')
pred = genfromtxt(r'1-17 full m=15 n350.csv', delimiter=',')
pred_x = pred.astype(int)
prediction = pred_x

# Load Target and Source Meshes
a = pv.read(r'horse17_partial2.off') # target

b = pv.read(r'horse1.off') # source

# MODERN IMAGEIO V3 CALL
texture = pv.numpy_to_texture(imageio.v3.imread(r'texture1.png')) #large squares
#texture = pv.numpy_to_texture(imageio.v3.imread(r'texture2.png')) #small squares

# Generate texture coordinates on target
a.texture_map_to_plane(inplace=True)

# MODERN TEXTURE COORDINATE TRANSFER (Via point_data)
b.active_texture_coordinates = a.active_texture_coordinates[prediction]

# Display
p = pv.Plotter(notebook=0, shape=(1,2)) # Automatically uses global white theme
p.add_mesh(a, texture=texture)
p.subplot(0,1)
p.add_mesh(b, texture=texture)

p.link_views()
p.show()