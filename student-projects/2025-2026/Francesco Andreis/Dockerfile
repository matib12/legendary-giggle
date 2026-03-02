FROM rootproject/root:6.34.00-ubuntu24.10

WORKDIR /app

# Copia tutto il contenuto della repository nel container
COPY . /app

# Esegui il programma compilato con ACLiC
CMD ["root", "-l", "-b", "-q", "geometry_and_physics_aau_15.cc+"]
