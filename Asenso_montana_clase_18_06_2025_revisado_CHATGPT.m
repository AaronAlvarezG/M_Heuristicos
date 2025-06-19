%%%%
%%%%Asenso_montaña ackley
%%%
%%% entrada
%%%
%% Problema
a=20;
b=.2;
c=2*pi;
dimenciones=2;
%%% parametros metodo 
n_iteraciones_cambio=1000;
epsilon=.000000001;
vecinos=20;
distancia_max=1;
corridas=30;
NEFO=zeros(corridas,1);
for corrida = 1:corridas
    sol = -100 + rand(1, dimenciones) * 200;
    objetivo(corrida,1) = evaluar(sol, dimenciones, a, b, c);
    NEFO(corrida,1) = 1;
    mejora = 0;

    while mejora < n_iteraciones_cambio
        sol_vecinos = zeros(vecinos, dimenciones);
        objetivo_vecino = zeros(vecinos, 1);

        for v = 1:vecinos
            sol_vecinos(v,:) = sol;
            seleccionar = randi(dimenciones);
            cambio = -distancia_max + rand() * (2 * distancia_max);
            sol_vecinos(v,seleccionar) = min(max(sol_vecinos(v,seleccionar) + cambio, -100), 100);
            objetivo_vecino(v) = evaluar(sol_vecinos(v,:), dimenciones, a, b, c);
            NEFO(corrida,1) = NEFO(corrida,1) + 1;
        end

        [min_val, idx_mejor] = min(objetivo_vecino);
        if min_val < objetivo(corrida,1)
            anterior = objetivo(corrida,1);
            sol = sol_vecinos(idx_mejor,:);
            objetivo(corrida,1) = min_val;
            mejora = abs(anterior - min_val) < epsilon ? mejora + 1 : 0;
        else
            mejora = mejora + 1;
        end
    end
end


function [ob]=evaluar(sol,dimenciones,a,b,c)
      ob=-a*(exp(-b*sqrt((1/dimenciones)*sum(sol.^2))))-exp((1/dimenciones)*sum(cos(c*sol)))+a+exp(1);
    
end