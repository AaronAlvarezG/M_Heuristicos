%%%%
%%%%Asenso_montaña ackley
%%%
%%% entrada
%%%
%% Problema


contador_evaluaciones = 0;

clear all
close all
a=20;
b=.2;
c=2*pi();
dimenciones=2;
%%% parametros metod
ISC=1000;
NIACT=100;
TI=3000;
TF=.01;
alpha=.99;
vecinos=5;
distancia_max=2;
corridas=30;
for corrida=1:corridas
    TA=TI;
    %%%% generar una solucion inial
    for j=1:dimenciones
        sol(corrida,j)=-100+rand()*200;
    end
    contador=0;
    contador_1=0;
    objetivo(corrida,1)=evaluar(sol(corrida,:),dimenciones,a,b,c);
    %%% mejo_sol_visitada
    ms(corrida,:)=sol(corrida,:);
    mso(corrida,:)=objetivo(corrida,1);
    while TA>TF
        for c=1:NIACT
            for v=1:vecinos
                sol_vecionos(v,:)=sol(corrida,:);
                cambio=-distancia_max+rand(1,dimenciones)*(2*distancia_max);
                sol_vecionos(v,:)=sol_vecionos(v,:)+cambio;
                for k=1:dimenciones
                    if sol_vecionos(v,k)<-100
                        sol_vecionos(v,k)=-100;
                    end
                    if sol_vecionos(v,k)>100
                        sol_vecionos(v,k)=100;
                    end
                end
                objetivo_vecino(v,1)=evaluar(sol_vecionos(v,:),dimenciones,a,b,c);
            end
            %%%% Criterio metropoli
            [ob_v, pos_v]=min(objetivo_vecino);
            if ob_v <= objetivo(corrida,1)
                objetivo(corrida,1)=ob_v;
                sol(corrida,:)=sol_vecionos(pos_v,:);
            else
                probabilidad=exp(-((ob_v-objetivo(corrida,1))/TA));
                des_1=rand();
                if des_1 < probabilidad
                    objetivo(corrida,1)=ob_v;
                    sol(corrida,:)=sol_vecionos(pos_v,:);
                end
            end
            if objetivo(corrida,1)< mso(corrida,:)
                ms(corrida,:)=sol(corrida,:);
                mso(corrida,:)=objetivo(corrida,1);
                contador_1=0;
            else
                contador_1=contador_1+1;
            end
            if contador_1>ISC
                sol(corrida,:)=ms(corrida,:);
                objetivo(corrida,1)=mso(corrida,:);
            end
            contador=contador+1;
            soluciones(contador,corrida)=objetivo(corrida,1);
        end
        TA=alpha*TA;
    end
end
function [ob]=evaluar(sol,dimenciones,a,b,c, contador_evaluaciones)

contador_evaluaciones = contador_evaluaciones + 1;

ob=-a*(exp(-b*sqrt((1/dimenciones)*sum(sol.^2))))-exp((1/dimenciones)*sum(cos(c*sol)))+a+exp(1);
end

% Estadísticas de resultados:
fprintf('Promedio: %.4f\n', mean(objetivo));
fprintf('Mínimo: %.4f\n', min(objetivo));
fprintf('Máximo: %.4f\n', max(objetivo));
fprintf('Desviación estándar: %.4f\n', std(objetivo));

histogram(objetivo)
xlabel('Valor de función Ackley')
ylabel('Frecuencia')
title('Distribución de soluciones en 30 corridas')

corridas_fallidas = find(objetivo >= 20);
fprintf('Corridas que no mejoraron: %s\n', mat2str(corridas_fallidas));

fprintf('Cantidad de llamadas a la funcion objetivo: %.4f\n', );