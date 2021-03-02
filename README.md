# Velocity limiter
Pythonban írt kód, mely meghatározott struktúrájú (.csv) fájlt beolvas, gyorsulás - értékeket\n
(külön gyorsulás és lassulás) limitál, ezzel megváltoztatva a sebesség - értékeket.

Fájl megnyitása/beolvasása az "Open file" gomb megnyomása után felugró ablakban történik. (Előre beállított (.csv) filter)
A gyorsulás/lassulás értékek 0 - 5 [m/s^2] között változhatnak, értéküket csúszkán állíthatjuk be.

"Generate!" gomb megnyomásával a program kirajzoltatja az aktuálisan megnyitott fájl adatait a beállított limitekkel.

Ha a beállított értékek megfelelnek számunkra, akkor a "Save file" gomb megnyomásával menthetjük tetszőleges néven, tesztőleges helyre
a módosított (.csv) fájlt.

GUI létrehozásához a PyQt5 modult használtam.

Diagramok rajzolásához pedig a Matplotlib modult.
