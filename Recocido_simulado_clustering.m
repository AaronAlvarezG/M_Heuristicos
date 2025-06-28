% Recocido simulado para resilver clustering con optimización de k
datos = readtable('winequality-white.csv');
% Normalizamos los valores de vinos
X = table2array(normalize(datos)); % Normalizamos los valores de vinos
n = size(X,1);

% Rango de k (número de clústeres)
k_min = 2;
k_max = 10;

% Parámetros para Recocido Simulado
T_inicial = 2000;
T_final = 1e-3;
alpha = 0.95; % factor de enfriamiento

% Parámetros de asenso en la montaña
max_inter_sin_mejora = 100;
epsilon = .000000001;
no_vecinos = 5;
no_corridas = 5;

% Incializar estructuras para evaluación
eval_por_corrida_y_k = zeros(no_corridas, k_max - k_min + 1);
mejor_objetivo_por_k = zeros(k_max - k_min + 1, 1);
asignaciones = cell(k_max - k_min + 1,1);
tiempos_por_k = zeros(k_max - k_min + 1,1);

for k=k_min: k_max
    col_k = k - k_min + 1;
    tic;

    puntajes_corrida = zeros(no_corridas, 1); % Inicializar el vector de objetivos para cada corrida
    
    mejor_puntaje = inf;
    mejor_asignacion = [];

    for corrida=1:no_corridas
        % generar una solucion inial
        asignacion_actual = randi(k,n,1);
        puntajes_corrida(corrida,1)=evaluar(X, asignacion_actual, k); % Evaluar la calidad de esa solución
        eval_por_corrida_y_k(corrida,col_k) = eval_por_corrida_y_k(corrida,col_k) + 1;
        
        % Busca meejor vecino
        T = T_inicial;
        while T < T_final
            
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
                
                if abs(puntajes_corrida(corrida) - mejor_valor) < epsilon
                    mejor_puntaje = puntajes_corrida(corrida);
                    mejor_asignacion = asignacion_actual; % Actualizar la mejor asignación encontrada

                    mejora = mejora + 1;
                else
                    mejora=0;
                end
                puntajes_corrida(corrida)=mejor_valor;
            else
                mejora=1+mejora;
            end

            T = T*alpha;
        end
    end
    
    mejor_objetivo_por_k(col_k) = min(puntajes_corrida); % Guardar el mejor objetivo para el k actual
    asignaciones{col_k} = mejor_asignacion;
    
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

%% Visualización: Método del codo
% Esta gráfica muestra cómo disminuye la distorsión intra-cluster a medida que aumenta k.
% Permite identificar el "codo", que sugiere un buen valor de k (trade-off entre complejidad y calidad).

valores_k = k_min:k_max;

figure;
plot(valores_k, mejor_objetivo_por_k, '-o', 'LineWidth', 2);
xlabel('Número de clústeres (k)');
ylabel('Distorsión intra-clúster');
title('Método del codo para selección de k');
grid on;

%% Visualización PCA con última asignación (de k = 4, por ejemplo)
% Proyectamos los datos a 2D con PCA para visualizar la asignación de clústeres.

k = 4;
k_visual = k - k_min + 1 ;
[~, X_pca] = pca(X);              % PCA: todas las componentes
X_reducido = X_pca(:, 1:2);       % Nos quedamos con las dos principales

figure;
gscatter(X_reducido(:,1), X_reducido(:,2), asignaciones{k_visual}); % asignación generada para último k evaluado
xlabel('Componente principal 1');
ylabel('Componente principal 2');
title(['Visualización PCA para k = ', num2str(k)]);
legend('show');
grid on;

%% Comparación visual para múltiples valores de k
% Compara en subplots cómo se agrupan los datos para diferentes valores de k

valores_k_visuales = [3, 4, 5, 6];   % Puedes ajustar los valores según los resultados
figure;

for i = 1:length(valores_k_visuales)
    k = valores_k_visuales(i);
    subplot(2, 2, i);  % Mosaico 2x2
    gscatter(X_reducido(:,1), X_reducido(:,2), asignaciones{k});
    title(['k = ', num2str(k)]);
    xlabel('PC1'); ylabel('PC2');
    axis tight;
    grid on;
end

sgtitle('Comparación de clustering con diferentes valores de k');

%% Evaluaciones promedio por valor de k
% Esta gráfica muestra cuántas veces se evaluó la función por k en promedio

promedio_eval_por_k = mean(eval_por_corrida_y_k, 1);  % Promedio por columna (corridas)

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

