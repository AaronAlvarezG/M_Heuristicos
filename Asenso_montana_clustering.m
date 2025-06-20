% Asenso a la montaña para resilver clustering
% Cargar el dataset de vinos
vinos = readtable('wine_dataset.csv');
% Normalizamos los valores de vinos
vinos_normalizada = normalize(vinos); % Normalizamos los valores de vinos

% Definir los parámetros globlaes
k_min = 2;
k_max = 10;
max_repeticiones_k = 200;

% Definir los parámetros de asenso en la montaña
n_iteraciones_cambio=1000;
epsilon=.000000001;
vecinos=20;
distancia_max=1;
corridas=30;
NEFO=zeros(corridas,1);
for k=k_min: k_max

    for corrida=1:corridas
        %%%% generar una solucion inial
        sol=randi(k,height(vinos),1);
        
        objetivo(corrida,1)=evaluar(vinos, sol, k);
        NEFO(corrida,1)=NEFO(corrida,1)+1;
        %%%%%
        % Busca meejor vecino
        mejora=0;
        while mejora < n_iteraciones_cambio
            sol_vecionos=[];
            objetivo_vecino=[];
            for v=1:vecinos
                sol_vecionos(v,:)=sol(corrida,:);
                seleccionar=round(1+rand()*(dimenciones-1));
                cambio=-distancia_max+rand()*(2*distancia_max);
                sol_vecionos(v,seleccionar)=sol_vecionos(v,seleccionar)+cambio;
                if sol_vecionos(v,seleccionar)<-100
                    sol_vecionos(v,seleccionar)=-100;
                end
                if sol_vecionos(v,seleccionar)>100
                    sol_vecionos(v,seleccionar)=100;
                end
                objetivo_vecino(v,1)=evaluar(sol_vecionos(v,:),dimenciones,a,b,c);
                NEFO(corrida,1)=NEFO(corrida,1)+1;
            end
            [a1,a2]=min(objetivo_vecino);
            if a1<objetivo(corrida,1)
                %%% remplazo 
                sol(corrida,:)=sol_vecionos(a2,:);
                objetivo(corrida,1)=objetivo_vecino(a2,1);
                if objetivo(corrida,1)-a1<epsilon
                    mejora=mejora+1;
                else
                    mejora=0;
                end
            else
               mejora=1+mejora;
            end
        end
    end
end

function [J]=evaluar(datos, soluciones,k)
    J = 0;
    % Inicializar los centroides
    centroide = cell(k, 1);    
    % Para cada j hasta k en soluciones extraer todos lo puntos asignados al cluster j
    for j = 1:k
        clusterPoint = datos(soluciones == j, :);
        % Si no hay puntos penalizar J con un valor muy grande
        if isempty(clusterPoint)
            J = J + 1e6; % Penalizar si no hay puntos en el cluster
        else
            centroide{j} = mean(clusterPoint);
            J = J + sum(pdist2(clusterPoint, centroide{j})); % Calcular la suma de distancias
        end
        
    end

    
end

