% Implemtanción de GRASP para obtener el óptimo k en clusterisación

% Paso 1 : Construcción - se genera una solución incial de manera greedy aleatorizada
% Paso 2 : Búsqueda local - se mejora la solución inicial con un método local .
% Paso 3 : Selección - si la solución mejorada es mejor que la global, se actualiza.

% Fase de contrucción: 
%  a) empieza con una solución vacía
%  b) calcula todos los candidatos posibles para añadir
%  c) evalua el costo de cada candidato
%  d) construye una lista restringida de candidatos (RCL) con los mejores candidatos según un umbral
%  e) escoge aleatoriamente uno de la RCL y lo añade a la solución
%  f) repite hasta completar la solución

% Fase de búsqueda local:
%  - Explorar el vecindario de la solución inicial
%  - Si se encuentra una solución mejor se mueve a ella y se sigue buscando
%  - Cuando no hay mejoras se detiene -> óptimo local

% Criterio de parada:
%  - número máximo de iteraciones
%  - Tiempo máximo de ejecución
%  - No hay mejora después de varias iteraciones


% ---------- CONFIGURACIÓN ----------
archivo = 'winequality-white.csv';
k_min = 2; k_max = 10;
num_corridas = 10;         % <= número de corridas independientes
semillas = 100 + (1:num_corridas);  % o define tus propias semillas

% Parámetros GRASP 
grasp_iters         = 20;      % iteraciones GRASP por k
alpha_rcl           = 0.30;
max_iter_local      = 20;
mejora_minima_local = 1e-6;
penalizacion_vacio  = 1e6;

% ---------- CARGA Y NORMALIZACIÓN ----------
T = readtable(archivo);
X = table2array(normalize(T));
valores_k = k_min:k_max;
nK = numel(valores_k);

% Matriz de resultados: filas = corridas, cols = k
sse_por_corrida = nan(num_corridas, nK);

% ---------- CORRIDAS ----------
for run = 1:num_corridas
    rng(semillas(run), 'twister');
    sse_por_corrida(run, :) = grasp_una_corrida_k_sweep( ...
        X, k_min, k_max, grasp_iters, alpha_rcl, ...
        max_iter_local, mejora_minima_local, penalizacion_vacio);
    fprintf('Corrida %d/%d terminada.\n', run, num_corridas);
end

% ---------- ESTADÍSTICAS ----------
mediana_sse = median(sse_por_corrida, 1, 'omitnan');
q1_sse      = quantile(sse_por_corrida, 0.25, 1);
q3_sse      = quantile(sse_por_corrida, 0.75, 1);

% ---------- GRÁFICO (mediana + IQR) ----------
figure;
plot(valores_k, mediana_sse, '-o', 'LineWidth', 2); hold on; grid on;
% Banda IQR
x_fill = [valores_k, fliplr(valores_k)];
y_fill = [q1_sse, fliplr(q3_sse)];
fill(x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(valores_k, mediana_sse, '-o', 'LineWidth', 2);
uistack(findobj(gca,'Type','line','-and','LineWidth',2),'top');

xlabel('Número de clústeres (k)');
ylabel('Distorsión intra-clúster (SSE)');
title('GRASP con múltiples corridas (mediana e IQR)');
legend({'IQR (Q1–Q3)', 'Mediana'}, 'Location','northeast');

% Parámetros en la figura
txt = sprintf('\\alpha_{RCL}=%.2f, itersLocal=%d, mejora=%.1e, pena=%.1e, iters/k=%d, corridas=%d', ...
    alpha_rcl, max_iter_local, mejora_minima_local, penalizacion_vacio, grasp_iters, num_corridas);
text(valores_k(1), mediana_sse(1), ['\leftarrow ', txt], ...
     'FontSize', 9, 'Interpreter','tex', 'VerticalAlignment','bottom');

% Guardar EPS 300 dpi
print(gcf, '-depsc', '-r300', 'grasp_multicorridas_elbow.eps');

% ---------- GUARDAR TABLA ----------
tabla = table(valores_k(:), mediana_sse(:), q1_sse(:), q3_sse(:), ...
    'VariableNames', {'k','SSE_mediana','SSE_Q1','SSE_Q3'});
writetable(tabla, 'grasp_multicorridas_resumen.csv');
disp('Guardados: grasp_multicorridas_elbow.eps y grasp_multicorridas_resumen.csv');


% --------------------------------------------------------------
function mejores_J_por_k = grasp_una_corrida_k_sweep(X, k_min, k_max, grasp_iters, alpha_rcl, max_iter_local, mejora_min, pen_vacio)
    valores_k = k_min:k_max;
    mejores_J_por_k = inf(1, numel(valores_k));

    pos = 1;
    for k = valores_k
        mejor_J_k = inf;

        for it = 1:grasp_iters
            % Construcción tipo k-means++ con RCL
            centroides_ini = construir_centroides_grasp(X, k, alpha_rcl);
            etiquetas = asignar_a_centroides(X, centroides_ini);

            % Búsqueda local
            [etiquetas, J] = busqueda_local_labels(X, etiquetas, k, max_iter_local, mejora_min, pen_vacio);

            if J < mejor_J_k
                mejor_J_k = J;
            end
        end

        mejores_J_por_k(pos) = mejor_J_k;
        pos = pos + 1;
    end
end

% ==== Auxiliares: mismas que ya usas ====
function centroides = construir_centroides_grasp(X, k, alpha_rcl)
    [n, d] = size(X);
    centroides = zeros(k, d);
    idx = randi(n); centroides(1,:) = X(idx,:);
    dist_min = inf(n,1);
    for c = 2:k
        dist_min = min(dist_min, sum((X - centroides(c-1,:)).^2, 2));
        dmin = min(dist_min); dmax = max(dist_min);
        umbral = dmin + alpha_rcl * (dmax - dmin);
        cand = find(dist_min >= umbral); if isempty(cand), cand = 1:n; end
        idx = cand(randi(numel(cand)));
        centroides(c,:) = X(idx,:);
    end
end

function etiquetas = asignar_a_centroides(X, centroides)
    k = size(centroides,1); n = size(X,1);
    etiquetas = ones(n,1);
    for i = 1:n
        d2 = sum((centroides - X(i,:)).^2, 2);
        [~, j] = min(d2); etiquetas(i) = j;
    end
end

function [etiquetas, J] = busqueda_local_labels(X, etiquetas, k, max_iter_local, mejora_min, pen_vacio)
    J = f_objetivo_SSE(X, etiquetas, k, pen_vacio);
    n = size(X,1);
    for it = 1:max_iter_local
        hubo_mejora = false;
        for i = 1:n
            etiqueta_actual = etiquetas(i);
            mejor_delta = 0; mejor_cluster = etiqueta_actual;
            for c = 1:k
                if c == etiqueta_actual, continue; end
                if sum(etiquetas == etiqueta_actual) == 1, break; end % no vaciar
                etiquetas_trial = etiquetas; etiquetas_trial(i) = c;
                J_trial = f_objetivo_SSE(X, etiquetas_trial, k, pen_vacio);
                delta = J_trial - J;
                if delta < mejor_delta
                    mejor_delta = delta; mejor_cluster = c;
                end
            end
            if mejor_cluster ~= etiqueta_actual
                etiquetas(i) = mejor_cluster;
                J = J + mejor_delta; hubo_mejora = true;
            end
        end
        if ~hubo_mejora || abs(mejor_delta) < mejora_min, break; end
    end
end

function J = f_objetivo_SSE(X, etiquetas, k, pen_vacio)
    etiquetas = etiquetas(:);
    J = 0;
    for c = 1:k
        pts = X(etiquetas == c, :);
        if isempty(pts)
            J = J + pen_vacio;
        else
            cen = mean(pts,1);
            dif = pts - cen;
            J = J + sum(sum(dif.^2));
        end
    end
end
