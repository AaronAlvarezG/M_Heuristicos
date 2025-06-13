% prueba de Kolmogorov-Smirnov (K-S) 
% Generar datos de ejemplo
data1 = randn(100, 1); % Datos de la primera muestra
data2 = randn(100, 1) + 0.5; % Datos de la segunda muestra con desplazamiento

funcion_FDA = @(x) normcdf(x, 6, 2); % Función anónima

% Paso 1: Ordenar la muestra y unirlos en un solo vector
valores = sort([data1; data2]);

% Paso 2: Calcular las ECDF
ECDF1 = arrayfun(@(x) sum(data1 <= x)/length(data1), valores);
ECDF2 = arrayfun(@(x) sum(data2 <= x)/length(data2), valores);

% Paso 3: Calcular la CDF teórica (o ECDF de la segunda muestra)
%x = linspace(min(data2), max(data2), 100); % Rango para la CDF teórica
cdfTeorica = funcion_FDA(valores); % Calcular la CDF teórica usando la función anónima

% Paso 4: Calcular las diferencias
diferencias_vs_dis_norm = abs(ECDF1 - cdfTeorica); % Calcular las diferencias absolutas entre ECDF y CDF teórica
diferencias_vs_dist = abs(ECDF1 -ECDF2); % Calcular las diferencias absolutas las ECDF de las dos distribusciones 

%  Paso 5: Obtener el estadístico D
D1 = max(diferencias_vs_dis_norm);
D2 = max(diferencias_vs_dist);% Calcular el estadístico D para la segunda muestra

% Paso 6: Calcular el valor crítico o p-valor 
alpha = 0.05; % Nivel de significancia

valorCritico_1 = 1.36/sqrt(length(data1));
valorP_1 = 1 - normcdf(D1, 0, sqrt(1/length(data1))); % Calcular el p-valor

valorCritico_2 = 1.36/sqrt(length(data2));
valorP_2 = 1 - normcdf(D2, 0, sqrt(1/length(data2))); % Calcular el p-valor

disp(' ');

% Comparando el Conjunto 1 con una distribución normal
disp('Conjunto 1 VS Distribución normal');
% Paso 6: Tomar decisión mediante valorP
if valorP_1 < alpha
    disp('    Prueba valorP: H0 Rechazada');
else
    disp('    Prueba valorP: H0 Aceptada');
end

% Paso 6: Tomar decisión mediante valor crírico
if D1 > valorCritico_1
    disp('    Prueba Valor crítico: Rechazada');
else
    disp('    Prueba Valor crítico: Aceptada');
end

% Comparando el Conjunto 1 con ek Conjunto 2
disp('Conjunto 1 VS Conjunto 2')
% Paso 6: Tomar decisión mediante valorP
if valorP_2 < alpha
    disp('    Prueba valorP: H0 Rechazada');
else
    disp('    Prueba valorP: H0 Aceptada');
end

% Paso 6: Tomar decisión mediante valor crírico
if D2 > valorCritico_2
    disp('    Prueba Valor crítico: Rechazada');
else
    disp('    Prueba Valor crítico: Aceptada');
end