# Format i nom de la imatge
set term pngcairo enhanced font 'Verdana,9'
set output "P9-1920-fig-map-2D.png"

# Títol del gràfic
set title "Mapa de Temperatura"

# Plot 
plot "aux3.dat" w image
