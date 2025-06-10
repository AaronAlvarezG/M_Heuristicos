% Problema del conjunto mínimo cubridor - Ascenso en la montaña

% Matriz de adyacencia simétrica A (grafo no dirigido)
A = [0 1 0 0 1
     1 0 1 0 0
     0 1 0 1 0
     0 0 1 0 1
     1 0 0 1 0];

[a1, ~] = size(A); % Número de nodos
penalizacion = 1000;

% Generar una solución aleatoria: 1 si el nodo es seleccionado
sol_1 = round(rand(1, a1));

% Evaluar cuántos nodos están cubiertos
contador = 0;
for i = 1:a1
    cubierto = sol_1(i); % está cubierto si se selecciona a sí mismo
    for j = 1:a1
        if A(i,j) == 1 && sol_1(j) == 1
            cubierto = 1; % está cubierto por un vecino
            break;
        end
    end
    if cubierto == 0
        contador = contador + 1; % nodo i no está cubierto
    end
end

% Función objetivo: minimizar nodos seleccionados + penalización
numero_nodos_seleccionados = sum(sol_1);
funcion_objetivo = numero_nodos_seleccionados + penalizacion * contador;

disp(numero_nodos_seleccionados)
disp(funcion_objetivo)

