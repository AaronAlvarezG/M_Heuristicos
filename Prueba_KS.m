% prueba de Kolmogorov-Smirnov (K-S) 
% Generar datos de ejemplo
% data1 = randn(100, 1); % Datos de la primera muestra
% data2 = randn(100, 1) + 0.5; % Datos de la segunda muestra con desplazamiento

function KS = Prueba_KS(data1, data2)

    mu = 6;
    sigma = 2;
    
    funcion_FDA = @(x) normcdf(x, mu, sigma); % Función anónima
    
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
    n1 = length(data1);
    n2 = length(data2);
    
    valorCritico_1 = 1.36/sqrt(n1);
    % valorP_1 = 1 - normcdf(D1, 0, sqrt(1/length(data1))); % Calcular el p-valor
    
    valorCritico_2 = 1.36*sqrt((n1 + n2)/(n1*n2));
    % valorP_2 = 1 - normcdf(D2, 0, sqrt(1/length(data2))); % Calcular el p-valor
    
    
    % Paso 5: Calcular p-valores (aproximación asintótica)
    
    % % Prueba de una muestra
    % n1 = length(data1);
    % lambda1 = (sqrt(n1) + 0.12 + 0.11/sqrt(n1)) * D1;
    % valorP_1 = 2 * exp(-2 * lambda1^2);
    % 
    % % Prueba de dos muestras
    % n2 = length(data2);
    % n_eff = sqrt(n1 * n2 / (n1 + n2));
    % lambda2 = (n_eff + 0.12 + 0.11/n_eff) * D2;
    % valorP_2 = 2 * exp(-2 * lambda2^2);
    % 
    % % Paso 6: Calcular valores críticos aproximados (opcional)
    % alpha = 0.05;
    % valorCritico_1 = 1.36 / sqrt(n1);
    % valorCritico_2 = 1.36 * sqrt((n1 + n2) / (n1 * n2));
    
    disp('======================|| PRUEBA DE KOLMOGOR-SMIRNOV ||============================================= ');
    
    % Comparando el Conjunto 1 con una distribución normal
    disp('Conjunto 1 VS Distribución normal');
    
    % Paso 6: Tomar decisión mediante valorP
    % fprintf('    Prueba valorP: ¿Es %d menor que %d ?',valorP_1,alpha);
    % if valorP_1 < alpha
    %     disp('    SI -> H0 Rechazada');
    % else
    %     disp('    No -> H0 Aceptada');
    % end
    
    % Paso 6: Tomar decisión mediante valor crírico
    fprintf('    Prueba Valor crítico: ¿Es %.4f mayor que %.4f ?',D1,valorCritico_1);
    
    if D1 > valorCritico_1
        disp('    SI -> H0 Rechazada');
    else
        disp('    NO -> H0 Aceptada');
    end
    
    % Comparando el Conjunto 1 con ek Conjunto 2
    disp('Conjunto 1 VS Conjunto 2')
    
    % Paso 6: Tomar decisión mediante valorP
    % fprintf('    Prueba valorP: ¿Es %d menor que %d ?',valorP_2,alpha);
    % if valorP_2 < alpha
    %     disp('    SI -> H0 Rechazada');
    % else
    %     disp('    NO -> H0 Aceptada');
    % end
    
    % Paso 6: Tomar decisión mediante valor crírico
    fprintf('    Prueba Valor crítico: ¿Es %.4f menor que %.4f ?',D2,valorCritico_2);
    if D2 > valorCritico_2
        disp('    SI -> H0 Rechazada');
    else
        disp('    NO -> H0 Aceptada');
    end
    
    x = sort(data1);
    Fx =  normcdf(x, mu, sigma);
    [h_oficial1, p_oficial1, ksstat_oficial1] = kstest(data1, 'CDF', [x, Fx]);   
    [h_oficial, p_oficial, ksstat_oficial] = kstest2(data1, data2);
    
    % Comparando
    fprintf('Resultado Normal: oficial: %.4f\n', h_oficial1);
    fprintf('Resultado 2 Conjuntos: oficial: %.4f\n', h_oficial);
    
    fprintf('Normal: oficial: %.4f\n', ksstat_oficial1);
    fprintf('Normal: manual:  %.4f\n', D1);
    
    fprintf('2 Conjuntos: oficial: %.4f\n', ksstat_oficial);
    fprintf('2 Conjuntos: manual:  %.4f\n', D2);
    
    % fprintf('P-valor oficial: %.4f\n', p_oficial);
    % fprintf('P-valor manual:  %.4f\n', valorP_2);
end