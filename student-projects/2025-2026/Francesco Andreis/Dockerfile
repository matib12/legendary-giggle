# ROOT 6.34 su Ubuntu 24.10
FROM rootproject/root:6.34.00-ubuntu24.10

WORKDIR /app

# Copia il codice nel container
COPY geometry_and_physics_aau_15.cpp .

# Compila con ROOT
RUN g++ -O2 geometry_and_physics_aau_15.cpp -o sim \
    $(root-config --cflags --libs)

# Esegui il programma
CMD ["./sim"]
