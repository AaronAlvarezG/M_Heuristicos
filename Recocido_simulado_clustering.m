% Recocido simulado para resilver clustering con optimización de k
datos = readtable('winequality-white.csv');
% Normalizamos los valores de vinos
X = table2array(normalize(datos)); % Normalizamos los valores de vinos
n = size(X,1);

% Rango de k (número de clústeres)
k_min = 2;
k_max = 10;

% Parámetros para Recocido Simulado
% T_inicial = 2000;
% T_final = 1e-2;
% alpha = 0.95; % factor de enfriamiento
% no_iteraciones_por_temperatura = 50;
% no_vecinos = 5;
% no_corridas = 5;
% no_iteraciones_sin_mejora = 1000;

opcion = 2;  % Opciones: 'rapido', 'equilibrado', 'preciso'
switch lower(opcion)
        case 1 %'rapido'
            disp('rapido');
            T_inicial = 500;
            T_final = 1e-1;
            alpha = 0.85;
            no_iteraciones_por_temperatura = 10;
            no_vecinos = 3;
            no_corridas = 3;
            no_iteraciones_sin_mejora = 100;

        case 2 %'balanceado'
            disp('balanceado');
            T_inicial = 2000;
            T_final = 1e-2;
            alpha = 0.95;
            no_iteraciones_por_temperatura = 50;
            no_vecinos = 5;
            no_corridas = 5;
            no_iteraciones_sin_mejora = 1000;

        case 3 %'preciso'
            disp('preciso')
            T_inicial = 5000;
            T_final = 1e-4;
            alpha = 0.99;
            no_iteraciones_por_temperatura = 100;
            no_vecinos = 10;
            no_corridas = 10;
            no_iteraciones_sin_mejora = 3000;

        otherwise
            error('Opción no reconocida. Usa: ''rapido'', ''balanceado'' o ''preciso''.');
end
% Incializar estructuras para evaluación
eval_por_corrida_y_k = zeros(no_corridas, k_max - k_min + 1);
mejor_objetivo_por_k = zeros(k_max - k_min + 1, 1);
asignaciones = cell(k_max - k_min + 1,1);
tiempos_por_k = zeros(k_max - k_min + 1,1);

total_iteraciones = (k_max - k_min + 1) * no_corridas;
contador_global = 0;

disp('Avance: ');
hora_actual = datetime('now', 'Format', 'HH:mm:ss');
fprintf(' [%s]: %5.1f %% (k = %2d, corrida = %2d de %2d)\n', hora_actual, 100 * contador_global / total_iteraciones, k, corrida, no_corridas);
for k=k_min: k_max
    col_k = k - k_min + 1;
    tic;

    puntajes_corrida = zeros(no_corridas, 1); % Inicializar el vector de objetivos para cada corrida
    
    mejor_asignacion = [];

    for corrida=1:no_corridas
        % generar una solucion inial
        asignacion_actual = randi(k,n,1);
        puntajes_corrida(corrida,1)=evaluar(X, asignacion_actual, k); % Evaluar la calidad de esa solución
        eval_por_corrida_y_k(corrida,col_k) = eval_por_corrida_y_k(corrida,col_k) + 1;
        
        mejor_sol_corrida = asignacion_actual;
        mejor_valor_corrida = puntajes_corrida(corrida, 1);
        contador_sin_mejora = 1;

        % Busca meejor vecino
        T = T_inicial;
        while T > T_final
            for c = 1 : no_iteraciones_por_temperatura
                puntajes_vecinos=zeros(no_vecinos,1);
                vecinos = cell(no_vecinos, 1);
    
                for v=1:no_vecinos
                    vecino_candidato = asignacion_actual;
                    idx = randi(n);
                    nueva_etiqueta = randi(k);
                    while nueva_etiqueta == vecino_candidato(idx)
                        nueva_etiqueta = randi(k);
                    end
                    vecino_candidato(idx) = nueva_etiqueta;
                    
                    puntajes_vecinos(v)=evaluar(X, vecino_candidato, k);
                    vecinos{v} = vecino_candidato;
                    eval_por_corrida_y_k(corrida,col_k) = eval_por_corrida_y_k(corrida,col_k) + 1;
                end
                   
                % Buscar mínimo objetivo_vecino
                [mejor_valor, idx_mejor]=min(puntajes_vecinos);
                if mejor_valor <= puntajes_corrida(corrida)
                    asignacion_actual = vecinos{idx_mejor};
                    puntajes_corrida(corrida)=mejor_valor;
                else
                    delta = mejor_valor - puntajes_corrida(corrida);
                    prob = exp(-delta / T);
    
                    if rand() < prob
                        asignacion_actual = vecinos{idx_mejor};
                        puntajes_corrida(corrida)=mejor_valor;                        
                    end
                end    

                if puntajes_corrida(corrida) < mejor_valor_corrida
                    mejor_valor_corrida = puntajes_corrida(corrida);
                    mejor_sol_corrida = asignacion_actual;
                    contador_sin_mejora = 0;
                else
                    contador_sin_mejora = contador_sin_mejora + 1;
                end

                if contador_sin_mejora > no_iteraciones_sin_mejora
                    asignacion_actual = mejor_sol_corrida;
                    puntajes_corrida(corrida) = mejor_valor_corrida;

                    contador_sin_mejora = 0;
                end


            end
            T = T*alpha;
        end
        contador_global = contador_global + 1;
        hora_actual = datetime('now', 'Format', 'HH:mm:ss');
        fprintf(repmat('\b', 1, 50));
        fprintf(' [%s]: %5.1f %% (k = %2d, corrida = %2d de %2d)\n', hora_actual, 100 * contador_global / total_iteraciones, k, corrida, no_corridas);
        
    end
    
    mejor_objetivo_por_k(col_k) = min(puntajes_corrida); % Guardar el mejor objetivo para el k actual
    asignaciones{col_k} = mejor_sol_corrida;
    
    tiempos_por_k(col_k) = toc;
end


function [J] = evaluar(X, etiquetas, k)
    J = 0;
    % Inicializar los centroides
    centroide = cell(k, 1);    
    % Para cada j hasta k en soluciones extraer todos lo puntos asignados al cluster j
    for j = 1:k
        puntos_del_cluster = X(etiquetas == j, :);
        % Si no hay puntos penalizar J con un valor muy grande
        if isempty(puntos_del_cluster)
            J = J + 1e6; % Penalizar si no hay puntos en el cluster
        else
            centroide{j} = mean(puntos_del_cluster);
            diferencias = centroide{j} - puntos_del_cluster;
            J = J + sum(sum(diferencias.^2)); % Calcular la suma de distancias
        end
    end
end

%% Muestra el tiempo total de ejecución
disp('Tiempo total (minutos):');
disp(sum(tiempos_por_k)/60);

%% Gráfica: Desempeño del Recocido Simulado para distintos valores de k
figure;
plot(k_min:k_max, mejor_objetivo_por_k, '-o')
xlabel('Número de Clústeres (k)')
ylabel('Mejor objetivo encontrado')
title('Desempeño del Recocido Simulado para distintos valores de k')

%% Gráfica: Tiempo por valor de k
figure;
bar(k_min:k_max, tiempos_por_k)
xlabel('Número de Clústeres (k)')
ylabel('Tiempo de ejecución (segundos)')
title('Tiempo por valor de k')

%% Comparación visual para múltiples valores de k
% Compara en subplots cómo se agrupan los datos para diferentes valores de k

k_visual = 4 - k_min + 1 ;
[~, X_pca] = pca(X);              % PCA: todas las componentes
X_reducido = X_pca(:, 1:2);       % Nos quedamos con las dos principales

valores_k_visuales = [2, 3, 4, 5, 6, 7, 8, 9, 10];   % Puedes ajustar los valores según los resultados
figure;

for i = 1:length(valores_k_visuales)
    k = valores_k_visuales(i);
    subplot(3, 3, i);  % Mosaico 2x2
    gscatter(X_reducido(:,1), X_reducido(:,2), asignaciones{i});
    title(['k = ', num2str(k)]);
    xlabel('PC1'); ylabel('PC2');
    axis tight;
    grid on;
end

sgtitle('Comparación de clustering con diferentes valores de k');