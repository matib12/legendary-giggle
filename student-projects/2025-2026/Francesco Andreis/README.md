## Contenuto repository
Questa repository contiene la simulazione ROOT: geometry_and_physics_aau_15.cc
Il codice è stato testato con l'ultima versione di ROOT 6.34.08

## Funzionamento programma
Il programma è una simulazione di scattering elastico di 4He su 197Au a 15 MeV. Il programma costruisce la geometria dell'esperimento, simula gli eventi, randomizza la posizione all'interno del pixel del detector e tiene conto dei contributi dovuti alle perdite di energia nel target e negli spessori morti. Successivamente ogni evento è pesato con la sezione d'urto Rutherford per l'angolo a cui viene rivelato.
I file di uscita contengono un TTree che per ogni evento specifica energia e angoli in coordinate sferiche e anche una flag per determinare in quale rivelatore l'evento sia finito e degli istogrammi con lo spettro in energia per ogni pixel dei rivelatori.

## Run with Docker (ROOT 6.34.08)

Build the image:
docker build -t aau_sim .

Run the simulation:
docker run --rm -v $PWD:/app aau_sim

## Output
Il programma produce:
- sim_aau_15.root
- h_sim_aau_15.root
I file di output compariranno nella stessa cartella

## Note
Il tempo con cui il programma funziona è proporzionale agli eventi ninc che si vogliono simulare. Il problema è che per ottenere statistica sufficiente, soprattutto per i rivelatori lontani (BU,BD,BR,BL) è necessario che il parametro ninc sia "alto", ossia circa 1e09, il che rende il compile time molto lungo.
