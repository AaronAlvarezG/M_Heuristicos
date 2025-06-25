%%%%
%%%%Asenso_montaña ackley
%%%
%%% entrada
%%%
%% Problema
clear all
close all
a=20;
b=.2;
c=2*pi();
dimenciones=2;
%%% parametros metod
ISC=1000; % Iteraciones sin mejora permitidas antes de reiniciar a la mejor solución conocida.
NIACT=100;% Iteraciones internas por temperatura
TI=2000; % Temperatura inicial
TF=.01; % Temperatura final
alpha=.99; % Factor de enfriamiento
vecinos=5; % Número de vecinos generados por iteración
distancia_max=4; % Máximo desplazamiento aleatorio por dimensión
corridas=30; % Número de ejecuciones independientes
cont_eval = 0; % Cantidad de llamadas a la funcion evaluar

evaluaciones = zeros(corridas, 1);

for corrida=1:corridas
    TA=TI;
    %%%% generar una solucion inial
    for j=1:dimenciones
        sol(corrida,j)=-100+rand()*200;
    end
    contador=0;
    contador_1=0;
    [cont_eval, objetivo(corrida,1)]=evaluar(sol(corrida,:),dimenciones,a,b,c, cont_eval);
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
                [cont_eval, objetivo_vecino(v,1)]=evaluar(sol_vecionos(v,:),dimenciones,a,b,c, cont_eval);
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
    evaluaciones(corrida) = cont_eval;
    cont_eval = 0;
end
function [cont, ob]=evaluar(sol,dimenciones,a,b,c, contador)
    contador = contador + 1;
    ob=-a*(exp(-b*sqrt((1/dimenciones)*sum(sol.^2))))-exp((1/dimenciones)*sum(cos(c*sol)))+a+exp(1);
    cont = contador;
end
disp('============================================================')
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

fprintf('Evaluaciones promedio por corrida: %.2f\n', mean(evaluaciones)/100)

scatter(evaluaciones, objetivo)
xlabel('Evaluaciones de la función')
ylabel('Valor final de Ackley')
title('Eficiencia por corrida')

% 18,825,030
% 18,825,030
% 627,501