Distancia=[0 21.71931168 40.43651815 26.7211246 3.6138622 7.486360932 4.003598381 6.765153361 10.82819006 2.199386278 10.77398719
21.71931168 0 18.90154756 21.00564686 18.18714106 28.84084777 22.88317504 15.81425306 30.42388864 22.22471597 11.19402073
40.43651815 18.90154756 0 28.65859906 36.9807842 47.28491514 41.78284816 34.71510478 48.19058933 40.70181937 30.08896808
26.7211246 21.00564686 28.65859906 0 25.18333775 30.01701684 30.31469776 25.76436687 28.05080035 25.37502709 23.07995017
3.6138622 18.18714106 36.9807842 25.18333775 0 11.06576703 5.428590977 3.354355378 13.99418808 4.874105046 7.164719115
7.486360932 28.84084777 47.28491514 30.01701684 11.06576703 0 9.018425583 14.24406192 4.962307931 6.617439082 18.16366703
4.003598381 22.88317504 41.78284816 30.31469776 5.428590977 9.018425583 0 7.069943423 13.40624108 6.095088186 11.69562311
6.765153361 15.81425306 34.71510478 25.76436687 3.354355378 14.24406192 7.069943423 0 17.34548933 8.225825187 4.638156962
10.82819006 30.42388864 48.19058933 28.05080035 13.99418808 4.962307931 13.40624108 17.34548933 0 9.120142543 20.57368465
2.199386278 22.22471597 40.70181937 25.37502709 4.874105046 6.617439082 6.095088186 8.225825187 9.120142543 0 11.68587609
10.77398719 11.19402073 30.08896808 23.07995017 7.164719115 18.16366703 11.69562311 4.638156962 20.57368465 11.68587609 0];

Tamano_lista_tabu=3;
Numero_iteraciones=10;
Corridas=30;
[puntos,~]=size(Distancia);
historial_costos = zeros(Corridas, Numero_iteraciones);

for corrida = 1: Corridas
    % solicion_inicial
    visitdos=zeros(1,puntos);
    solucion=zeros(1,puntos);
    solucion(1,1)=1;
    visitdos(1,1)=1;
    for i=2:puntos
        seleccion=round(2+rand()*(puntos-2));
        while visitdos(1,seleccion)==1
            seleccion=round(2+rand()*(puntos-2));
        end
        solucion(1,i)=seleccion;
        visitdos(1,seleccion)=1;
    end
    costo_sol=objetivo(solucion,Distancia,puntos);
    %%%%% Movimientos posibles
    movimientos=0;
    lista=[];
    for i=2:puntos-1
        for j=1+i:puntos
            movimientos=movimientos+1;
            lista(movimientos,1)=i;
            lista(movimientos,2)=j;
            lista(movimientos,3)=0;
        end
    end
        
    for iteraciones=1:Numero_iteraciones
        % Se generan todos los vecinos candidatos
        Candidatos=[];
        c_candidatos=0;
        for i=2:puntos-1
            for j=1+i:puntos
                c_candidatos=c_candidatos+1;
                Candidatos(c_candidatos,:)=solucion(1,:);
                Candidatos(c_candidatos,i)=solucion(1,j);
                Candidatos(c_candidatos,j)=solucion(1,i);
            end
        end
        
        % Se evalua el costo de cada candidato con la función objetivo
        costo_sol_candidatos=ones(c_candidatos,1)*inf;
        for movimiento=1:c_candidatos
            if lista(movimiento,3)==0
                costo_sol_candidatos(movimiento,1)=objetivo(Candidatos(movimiento,:),Distancia,puntos);
            end
        end
    
        % Se elige el mejor candidato que no esté en la lista tabú y se
        % acualiza la solución si mejora
        [valor_mc,mejor_candito]=min(costo_sol_candidatos);
        if valor_mc < costo_sol
            costo_sol=valor_mc;
            solucion=Candidatos(mejor_candito,:);
            lista(mejor_candito,3)=Tamano_lista_tabu;
        end
        % Criterio aspiracional: permite aceptar un candidato tabú si mejora la
        % solución actual
        costo_sol_candidatos_p=ones(c_candidatos,1)*inf;
        for movimiento=1:c_candidatos
            if lista(movimiento,3) > 0
                costo_temp = objetivo(Candidatos(movimiento,:),Distancia,puntos);
                if costo_temp < costo_sol
                    costo_sol_candidatos_p(movimiento) = costo_temp;
                end
            end
        end
        [valor_mct,mejor_candito_t]=min(costo_sol_candidatos_p);
        if valor_mct < costo_sol
            costo_sol=valor_mct;
            solucion=Candidatos(mejor_candito_t,:);
            lista(mejor_candito_t,3)=Tamano_lista_tabu;
        end
    
        % Actualizar lista tabú (disminuir duración)
        for m = 1:size(lista,1)
            if lista(m, 3) > 0
                lista(m, 3) = lista(m, 3) - 1; 
            end
        end
        historial_costos(corrida, iteraciones) = costo_sol;
    end
end

%%% 
function [costo]=objetivo(solucion,Distancia,puntos)
    solucion=[solucion,1];
    costo=0;
    for i=1:puntos
        a=solucion(1,i);
        b=solucion(1,i+1);
        costo=costo+Distancia(a,b);
    end
end

%% Graficar la solucion

% Calcular promedio y mejor corrida por iteración
promedio = mean(historial_costos, 1);
mejor_corrida = min(historial_costos, [], 1);

% Graficar
figure;
plot(1:Numero_iteraciones, promedio, 'b-', 'LineWidth', 2); hold on;
plot(1:Numero_iteraciones, mejor_corrida, 'r--', 'LineWidth', 2);
xlabel('Iteraciones');
ylabel('Costo de la solución');
title('Evolución del costo por iteración');
legend('Promedio de corridas', 'Mejor corrida');
grid on;