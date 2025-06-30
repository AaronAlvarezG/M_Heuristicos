% Asenso a la montaña para resolver clustering con optimización de k
datos = readtable('winequality-white.csv');
% Normalizamos los valores de vinos
X = table2array(normalize(datos)); % Normalizamos los valores de vinos
n = size(X,1);

% Rango de k (número de clústeres)
k_min = 2;
k_max = 10;

% Parámetros de asenso en la montaña
max_inter_sin_mejora = 1000;
epsilon = .000000001; 
no_vecinos = 20;
no_corridas = 30;

% Incializar estructuras para evaluación
no_evaluaciones_por_corrida_y_k = zeros(no_corridas, k_max - k_min + 1);
mejores_puntajes_por_k = zeros(k_max - k_min + 1, 1);
mejores_asignaciones_por_k = cell(k_max - k_min + 1,1);
tiempos_por_k = zeros(k_max - k_min + 1,1);

mejores_asignaciones_por_corrida = cell(k_max - k_min + 1, no_corridas);
historial_de_puntajes_por_corrida = cell(k_max - k_min + 1, no_corridas);

% === Progreso general del proceso ===
total_iteraciones = (k_max - k_min + 1) * no_corridas;
contador_global = 0;
disp('Avance: ');
hora_actual = datetime('now', 'Format', 'HH:mm:ss');
fprintf(' [%s]: %5.1f %% (k = 00, corrida = 00 de %2d)\n', hora_actual, 100 * contador_global / total_iteraciones, no_corridas);

for k=k_min: k_max
    col_k = k - k_min + 1;
    tic;

    puntajes_corrida = zeros(no_corridas, 1); % Guardará el mejor puntaje por corrida
    
    mejor_puntaje = inf; % Inicializa con un valor muy alto
    mejor_asignacion = []; % Variable para guardar la mejor asignación de clústeres
    puntajes_iteracion = [];    

    for corrida=1:no_corridas
        % generar una solucion inicial
        asignacion_actual = randi(k,n,1);
        puntaje_actual = evaluar(X, asignacion_actual, k); % Evaluar la calidad de esa solución
        puntajes_corrida(corrida,1) = puntaje_actual;
        puntajes_iteracion(end + 1) = puntaje_actual;

        no_evaluaciones_por_corrida_y_k(corrida,col_k) = no_evaluaciones_por_corrida_y_k(corrida,col_k) + 1;
        
        % Busca mejor vecino
        contador_sin_mejora = 0;
        while contador_sin_mejora < max_inter_sin_mejora
            
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
                no_evaluaciones_por_corrida_y_k(corrida,col_k) = no_evaluaciones_por_corrida_y_k(corrida,col_k) + 1;
            end
            % Buscar mínimo objetivo_vecino
            [mejor_valor, idx_mejor]=min(puntajes_vecinos);
            if mejor_valor<puntajes_corrida(corrida)
                puntajes_corrida(corrida) = mejor_valor;
                asignacion_actual = vecinos{idx_mejor};
                puntajes_iteracion(end + 1) = mejor_valor;

                if puntajes_corrida(corrida) - mejor_valor < epsilon
                    contador_sin_mejora = contador_sin_mejora + 1;
                else
                    contador_sin_mejora=0;
                end 

                if mejor_valor < mejor_puntaje
                    mejor_puntaje = mejor_valor;
                    mejor_asignacion = asignacion_actual;
                end
            else
                contador_sin_mejora = contador_sin_mejora + 1;
            end             
        end
        % Almacenar el mejor puntaje y asignación de la corrida actual
        mejores_asignaciones_por_corrida{col_k} = mejor_asignacion;
        historial_de_puntajes_por_corrida{col_k, corrida} = puntajes_iteracion;

        contador_global = contador_global + 1;
        hora_actual = datetime('now', 'Format', 'HH:mm:ss');
        fprintf(repmat('\b', 1, 50));
        fprintf(' [%s]: %5.1f %% (k = %2d, corrida = %2d de %2d)\n', hora_actual, 100 * contador_global / total_iteraciones, k, corrida, no_corridas);
    end
    
    mejores_puntajes_por_k(col_k) = min(puntajes_corrida); % Guardar el mejor objetivo para el k actual
    mejores_asignaciones_por_k{col_k} = mejor_asignacion;
    
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

%% Visualización: Método del codo
% Esta gráfica muestra cómo disminuye la distorsión intra-cluster a medida que aumenta k.
% Permite identificar el "codo", que sugiere un buen valor de k (trade-off entre complejidad y calidad).

valores_k = k_min:k_max;

figure;
plot(valores_k, mejores_puntajes_por_k, '-o', 'LineWidth', 2);
xlabel('Número de clústeres (k)');
ylabel('Distorsión intra-clúster');
title('Método del codo para selección de k');
grid on;

%% Visualización PCA con última asignación (de k = 4, por ejemplo)
% Proyectamos los datos a 2D con PCA para visualizar la asignación de clústeres.

k_visual = 4 - k_min + 1 ;
[~, X_pca] = pca(X);              % PCA: todas las componentes
X_reducido = X_pca(:, 1:2);       % Nos quedamos con las dos principales

figure;
gscatter(X_reducido(:,1), X_reducido(:,2), mejores_asignaciones_por_k{k_visual}); % asignación generada para último k evaluado
xlabel('Componente principal 1');
ylabel('Componente principal 2');
title(['Visualización PCA para k = ', num2str(k_visual)]);
legend('show');
grid on;

%% Comparación visual para múltiples valores de k
% Compara en subplots cómo se agrupan los datos para diferentes valores de k

valores_k_visuales = [2, 3, 4, 5, 6, 7, 8, 9, 10];   % Puedes ajustar los valores según los resultados
figure;

for i = 1:length(valores_k_visuales)
    k = valores_k_visuales(i);
    subplot(3, 3, i);  % Mosaico 2x2
    gscatter(X_reducido(:,1), X_reducido(:,2), mejores_asignaciones_por_k{i});
    title(['k = ', num2str(k)]);
    xlabel('PC1'); ylabel('PC2');
    axis tight;
    grid on;
end

sgtitle('Comparación de clustering con diferentes valores de k');

%% Evaluaciones promedio por valor de k
% Esta gráfica muestra cuántas veces se evaluó la función por k en promedio

promedio_eval_por_k = mean(no_evaluaciones_por_corrida_y_k, 1);  % Promedio por columna (corridas)

figure;
plot(k_min:k_max, promedio_eval_por_k, '-o', 'LineWidth', 2);
xlabel('Número de clústeres (k)');
ylabel('Evaluaciones promedio de la función');
title('Evaluaciones necesarias por valor de k');
grid on;

%% Graficar el tiempo de ejecución por cada valor de k:

valores_k = k_min:k_max;

figure;
plot(valores_k, tiempos_por_k, '-o', 'LineWidth', 2);
xlabel('Número de clústeres (k)');
ylabel('Tiempo de ejecución (segundos)');
title('Tiempo de ejecución por valor de k');
grid on;

%% Visualizar el historial de puntajes para alguna corrida

k_ejemplo = 10; corrida_ejemplo = 1;

figure;
puntajes = historial_de_puntajes_por_corrida{k_ejemplo - k_min + 1, corrida_ejemplo};
plot(puntajes, '-o');
xlabel('Iteración');
ylabel('Función objetivo');
title(sprintf('Convergencia (k = %d, corrida = %d)', k_ejemplo, corrida_ejemplo));
