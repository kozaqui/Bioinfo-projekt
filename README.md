#Bioinfo-projekt
Projekt na zajęcia z Bioinformatyki na UJ w roku 2015/2016/ Poszukiwanie struktur ze zmodyfikowanym RNA w PDB

Plik z programem to *modRNAinPDB.py*. W obecnej wersji do działania potrzebuje pliku *Components-pub.cif*, umieszczonego w working directory. Jest to plik w formacie mmCIF opisujący monomery z bazy PDB. Jest on zbyt duży jak na GitHub (~200 Mb), więc zamieszczam do niego [link](http://ligand-expo.rcsb.org/dictionaries/Components-pub.cif). Jest to pierwszy plik ze [strony](http://ligand-expo.rcsb.org/ld-download.html) Ligand Expo Downloads (Chemical component dictionaries).

Program tworzy lokalną bazę danych **modRNA.db** za pomocą SQLite (czyli plik z bazą danycn można otworzyć np. za pomocą sqlitebrowser)