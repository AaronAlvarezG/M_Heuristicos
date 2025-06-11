% prueba de Kolmogorov-Smirnov (K-S) 
% Generar datos de ejemplo
data1 = randn(100, 1); % Datos de la primera muestra
data2 = randn(100, 1) + 0.5; % Datos de la segunda muestra con desplazamiento

funcion_FDA = @(x) normcdf(x, 6, 2) % Función anónima

% Ordenar los datos
data1 = sort(data1)
data2 = sort(data2)

%




