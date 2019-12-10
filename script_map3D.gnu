# Format i nom de la imatge
set term pngcairo enhanced font 'Verdana,9'
set output "P9-1920-fig-map3D.png"

# Plot 
sp "aux3.dat" w pm3d
