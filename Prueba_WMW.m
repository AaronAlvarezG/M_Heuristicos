% Prueba Wilcoxon-Mann-Whitney

% Paso 1: Inicialización de datos
% Crear o cargar lo dos conjuntos de datos independientes
% grupo_a = randi(20,100,1);
% grupo_b = randi(20,100,1) + 10;

function WMW = Prueba_WMW(grupo_a, grupo_b)
    % Paso 2: Unión y ordenaminto
    grupo_a_ordenado = sort(grupo_a);
    grupo_b_ordenado = sort(grupo_b);
    
    union = [grupo_a_ordenado; grupo_b_ordenado];
    
    % Paso 3: Asignación de rangos
    rangos = tiedrank(union);
    
    % Paso 4: Sepración de rangos por grupo
    rangos_a = rangos(1:length(grupo_a));
    rangos_b = rangos(length(grupo_a)+1:end);
    
    % Paso 5: Cálculo del estadístico U
    R_a = sum(rangos_a);
    R_b = sum(rangos_b);
    
    na = length(grupo_a);
    nb = length(grupo_b);
    
    Ua = na*nb+(na*(na+1)/2)-R_a;
    Ub = na*nb+(nb*(nb+1)/2)-R_b;
    
    U_observado = min(Ua, Ub);
    % Paso 6: Comparación con valor crítio y p_valor
    mu_U = (na*nb)/2;
    sigma_U = sqrt((na * nb * (na + nb + 1)) / 12);
    z = (U_observado - mu_U) / sigma_U;
    
    % calcular p_valor
    p_valor = 2 * (1 - normcdf(abs(z)));
    
    % Paso 7: Interpretación
    fprintf('\n');
    disp('===============|| PRUEBA DE Wilcoxon-Mann-Whitney ||==================================');
    
    fprintf('\rResultados:');
    fprintf('   \nSuma de los rangos del gurpo A: %.4f', R_a);
    fprintf('   \nSuma de los rangos del gurpo B: %.4f', R_b);
    fprintf('   \nValor p: %.4f\n', p_valor);
    fprintf('   \nEstadístico U observado: %.4f', U_observado);
    
    if p_valor < 0.05
        fprintf('\nSe rechaza la hipótesis nula: hay diferencias significativas entre los grupos.\n');
    else
        fprintf('\nNo se rechaza la hipótesis nula: no hay diferencias significativas entre los grupos.\n');
    end
    
    % Paso 8: Visualización de los resultados
    figure;
    boxplot([grupo_a, grupo_b], 'Labels', {'Grupo A', 'Grupo B'});
    title('Comparación de los Grupos A y B');
    ylabel('Valores');
end
