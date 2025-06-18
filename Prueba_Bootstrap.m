% Prueba de Bootrstrap
% Paso 1: Definir dos conjuntos de datos aleatorios
% n = 100;
% grupo_a = randn(n,1);
% grupo_b = randn(n,1);

% Paso 2: Elegir el estadístico a comparar
% Media
% Varianza
% Desviación estándar
%

function PB = Prueba_Bootstrap(grupo_a, grupo_b)
    na = length(grupo_a);
    nb = length(grupo_b);

    fprintf('\n');
    disp('===============|| PRUEBA DE BOOTSTRAP ||==================================');
    
    disp('Seleccione el estadístico a comparar:');
    disp('  1) Media');
    disp('  2) Mediana');
    disp('  3) Desviación Estándar');
    
    opcion = input('Ingresa tu seleccion:');
    
    if opcion ==  1
        estadistico = @(x) mean(x);
        disp('  Estadístico seleccionado: media');
    elseif opcion == 2
        estadistico = @(x) median(x);
        disp('  Estadístico seleccionado: mediana');
    elseif opcion == 3
        estadistico = @(x) std(x); 
        disp('  Estadístico seleccionado: desviacione estándar');
    else
        error('  Error: Opción iválida');
    end
    
    % Paso 3: Calcular la diferencia observada
    estd_original_a = estadistico(grupo_a);
    estd_original_b = estadistico(grupo_b);
    
    diferencia_observada = estd_original_a - estd_original_b;
    
    fprintf('\n  Estadístico observado para el grupo A: %.4f', estd_original_a);
    fprintf('\n  Estadístico observado para el grupo B: %.4f', estd_original_b);
    
    fprintf('\n  Diferencia observada: %.4f', diferencia_observada);
    % Paso 4: incializar parametros de Bootrstrap
    tamano_muestra_B = 1000;
    diferencias_B = zeros(tamano_muestra_B,1);
    
    % Paso 5 Bucle de Bootstrap:
    
    for i = 1:tamano_muestra_B
        muestra_A = grupo_a(randi(na, na,1));
        muestra_B = grupo_b(randi(nb, nb,1));
        
        estadistico_temp_a = estadistico(muestra_A);
        estadistico_temp_b = estadistico(muestra_B);
        
        diferencias_B(i) = estadistico_temp_a - estadistico_temp_b;
    end
    % Paso 6: Calcular el intervalo de confianza
    IC = prctile(diferencias_B, [2.5, 97.5]); % intervalo de confinza del 95%
    
    % Paso 7: Evaluar si hay evidencia de diferencia
    % Si 0 pertence al Intervalo: no hay evidencia
    % otro: Si hay evidencia significativa
    
    fprintf('\n\n  Resultados:');
    fprintf('\n  El intervalo de confianza es: [%.4f, %.4f]', IC(1), IC(2));
    
    if IC(1) <= 0 && IC(2) >= 0
        disp(' por lo que NO hay evidencia suficiente de diferencia para rechazar H0');
    else
        disp(' por lo que existe evidencia significactiva para rechazar H0')
    end
    
    % (Opcional) Paso 8: Calcular el p-valor empírico
    % Calcualar la proporcion de diferencias de bootstrap 
    %   cuya magnitud sea mayor o igual a la diferencia observada. 
    
    conteo = sum(abs(diferencias_B) >= abs(diferencia_observada));
    p_valor = conteo / tamano_muestra_B;
    
    if p_valor <= 0.05
        fprintf('\n  p-valor= %.4f es menor que 0.05, por lo tanto se rechaza H0', p_valor);
    else
        fprintf('\n  p-valor= %.4f es mayor que 0.05, por lo tanto NO se rechaza H0', p_valor);
    end
    fprintf('\n');
    % Paso 9 (Opcional): Visualización
    % Generar un histograma de las diferecias de bootstrap
    % Dibujar una línea vertical con la diferencia observada
    % Resaltar el intervalo de confianza.
    hold on;
    histogram(diferencias_B);
    xline(diferencia_observada);
    
    xline(IC(1), '--r', 'IC inferior');
    xline(IC(2), '--r', 'IC superior');
    xline(0, '--k', 'Cero');
    legend('Distribución bootstrap', 'Diferencia observada', 'IC 95%');
    title('Distribución de diferencias bootstrap');
    xlabel('Diferencia de estadístico');
    ylabel('Frecuencia');
    hold off;
end