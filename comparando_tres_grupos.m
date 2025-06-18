% Cargar los 3 grupos de datos:
grupoA = load('grupo_A.txt');
grupoB = load('grupo_B.txt');
grupoC = load('grupo_C.txt');

Prueba_KS(grupoA, grupoB)
Prueba_KS(grupoA, grupoC)
Prueba_KS(grupoB, grupoC)

Prueba_WMW(grupoA, grupoB)
Prueba_WMW(grupoA, grupoC)
Prueba_WMW(grupoB, grupoC)

Prueba_Bootstrap(grupoA, grupoB)
Prueba_Bootstrap(grupoA, grupoC)
Prueba_Bootstrap(grupoB, grupoC)
