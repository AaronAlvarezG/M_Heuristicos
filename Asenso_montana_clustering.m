% Asenso a la montaña para resilver clustering con optimización de k
datos = readtable('winequality-white.csv');
% Normalizamos los valores de vinos
X = table2array(normalize(datos)); % Normalizamos los valores de vinos
n = size(X,1);

% Rango de k (número de clústeres)
k_min = 2;
k_max = 10;

% Parámetros de asenso en la montaña
max_inter_sin_mejora = 100;
epsilon = .000000001;
no_vecinos = 5;
no_corridas = 5;

% Incializar estructuras para evaluación
eval_por_corrida_y_k = zeros(no_corridas, k_max - k_min + 1);
mejor_objetivo_por_k = zeros(k_max - k_min + 1, 1);

for k=k_min: k_max
    col_k = k - k_min + 1;
    puntajes_corrida = zeros(no_corridas, 1); % Inicializar el vector de objetivos para cada corrida

    for corrida=1:no_corridas
        % generar una solucion inial
        asignacion_actual = randi(k,n,1);
        puntajes_corrida(corrida,1)=evaluar(X, asignacion_actual, k); % Evaluar la calidad de esa solución
        eval_por_corrida_y_k(corrida,col_k) = eval_por_corrida_y_k(corrida,col_k) + 1;
        
        % Busca meejor vecino
        mejora=0;
        while mejora < max_inter_sin_mejora
            
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
            if mejor_valor<puntajes_corrida(corrida)
                asignacion_actual=vecinos{idx_mejor};
                
                if abs(puntajes_corrida(corrida) - mejor_valor) < epsilon
                    mejora = mejora + 1;
                else
                    mejora=0;
                end
                puntajes_corrida(corrida)=mejor_valor;
            else
                mejora=1+mejora;
            end
        end
    end
    
    mejor_objetivo_por_k(col_k) = min(puntajes_corrida); % Guardar el mejor objetivo para el k actual
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

%% Visualización: método del codo

valores_k = k_min:k_max;

figure;
plot(valores_k, mejor_objetivo_por_k, '-o', 'LineWidth', 2);
xlabel('Número de clústeres (k)');
ylabel('Distorsión intra-clúster');
title('Método del codo para selección de k');
grid on;

k = 4;
[coef, X_pca] = pca(X);  % Reduce a todas las componentes principales
X_reducido = X_pca(:, 1:2);  % Tomamos solo las dos primeras para graficar

figure;
gscatter(X_reducido(:,1), X_reducido(:,2), asignacion_actual);
xlabel('Componente principal 1');
ylabel('Componente principal 2');
title(['Visualización PCA para k = ', num2str(k)]);
legend('show');
grid on;



[~, X_pca] = pca(X);             % Reducimos con PCA
X_reducido = X_pca(:, 1:2);      % Usamos las 2 primeras componentes

valores_k = [3, 4, 5, 6];  % Puedes modificar esta lista
figure;

for i = 1:length(valores_k)
    k = valores_k(i);
    subplot(2, 2, i);  % Organiza en una cuadrícula 2x2 (ajustable)
    gscatter(X_reducido(:,1), X_reducido(:,2), asignaciones{k});
    title(['k = ', num2str(k)]);
    xlabel('PC1'); ylabel('PC2');
    axis tight;
    grid on;
end

sgtitle('Comparación de clustering con diferentes valores de k');

promedio_eval_por_k = mean(eval_por_corrida_y_k, 1);
figure;
plot(k_min:k_max, promedio_eval_por_k, '-o', 'LineWidth', 2);
xlabel('Número de clústeres (k)');
ylabel('Evaluaciones promedio de la función');
title('Evaluaciones necesarias por valor de k');
grid on;