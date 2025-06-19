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
for corrida=1:corridas
    %%%% generar una solucion inial
    for j=1:dimenciones
        sol(corrida,j)=-100+rand()*200;
    end
    objetivo(corrida,1)=evaluar(sol(corrida,:),dimenciones,a,b,c);
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

function [ob]=evaluar(sol,dimenciones,a,b,c)
      ob=-a*(exp(-b*sqrt((1/dimenciones)*sum(sol.^2))))-exp((1/dimenciones)*sum(cos(c*sol)))+a+exp(1);
    
end