# Project name

   Database for Probe Station Measurements

# Problem description

   Attualmente ogni misura realizzata alla probe station in camera pulita è rappresentata da un file .txt, la cui struttura interna dipende dal tipo di misura realizzata (I-V, C-V, C-f, ...), e dunque per fare il plotting delle misure è necessario caricare tutti i 
file di interesse, risultando in una procedura time-consuming. Nella cartella data sono riportati degli esempi dei file .txt risultati da diverse tipologie di misura.

   L'idea del progetto è quella di realizzare una nuova forma di storage più efficiente dei dati, al fine di poter fare il loro retrieval più velocemente per il successivo plotting e analisi dati. 

   La soluzioni migliore a cui sono attualmente approdato è quella di realizzare un file .parquet per ciascuna tipologia di misura, quindi ci sarà ad esempio un IV.parquet e un CV.parquet. Il file .parquet è una tabella in cui ciascuna misura è rappresentata da una 
riga e le informazioni relative a quella misura sono disponibili nelle colonne. Ad esempio il CV.parquet per ciascuna misura ha le seguenti colonne: Produzione, Wafer1, Wafer2, Shot1, Shot2, Struttura, Temperatura, Frequenza, ACampl, C e V. Le ultime due colonne, 
C e V, sono dei vettori che contengono le capacità misurate ai diversi bias point applicati. Questo file CV.parquet è anche esso presente nella cartella data.



