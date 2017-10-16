function [TRANS_MOD] = matrixAdjust(estrato, goalPow, aparato, TRANSref, devPower, delta, length)
maxIter = 1000; % iteraciones por cada matrix de transmision
maxTries = 100; %maximo numero de variaciones en las matrices
EMISref = [1 0; 0 1]; % es la misma porque solo son dos estados y la emisión es el estado
tol = 0.05;
%T = TRANSref;

% if((aparato=='B')||(aparato=='O')||(aparato=='F'))
%     cambios=1; % cambian duracion y frecuencia
% else
%     if ((aparato=='C'))
%         cambios = 2; % cambia solo duracion
%     else
%         cambios = 3; % cambia solo frecuencia - menos o mas veces
%     end
% end

cambios = 2;

TRANStry = zeros(2);



for t=1:maxTries
    
    if t==1
        T = TRANSref;
        disp(TRANSref);
        % comprobar si hay que aumentar o disminuir el consumo
        if estrato > 4 % aumenta consumo
            subir=true;
        else
            if estrato <4 %bajar consumo
                subir = false;
            else
                TRANS_MOD = T;
                break
            end
        end       
    else
        T = TRANStry;
        disp(t)
        disp(TRANStry);
    end
    
    
    
    % comprobar si hay que aumentar o disminuir el consumo
    if subir==true % aumenta consumo
        switch cambios
            case 1
                TRANStry = [T(1)+ 0.5*delta  T(3)- 0.5*delta;T(2)-0.5*delta T(4) + 0.5*delta];
            case 2
                if (T(2)==1)||(T(4)==1)
                    disp('La matriz no puede ajustarse')
                    TRANS_MOD = T;
                    break
                else
                    TRANStry = [T(1)  T(3);T(2)-delta T(4)+delta];
                end
            otherwise
                TRANStry = [T(1)+ delta  T(3)- delta;T(2) T(4)];
        end
    else
        if subir==false %bajar consumo
            switch cambios
                case 1
                    TRANStry = [T(1)- 0.5*delta  T(3)+ 0.5*delta;T(2)+0.5*delta T(4)-0.5*delta];
                case 2
                    if (T(2)==1)||(T(4)==1)
                        disp('La matriz no puede ajustarse')
                        TRANS_MOD = T;
                        break
                    else
                        TRANStry = [T(1)  T(3);T(2)+delta T(4)-delta];
                    end
                otherwise
                    TRANStry = [T(1)+delta  T(3)-delta;T(2) T(4)];
            end
        else
            TRANS_MOD = T;
            break
        end
    end
    testDevicePow = zeros(1,maxIter);
    for i = 1:maxIter
        [seq, ~] = hmmgenerate(length,TRANStry,EMISref);
        %disp(unique(seq))
        test_normalizado = (seq - 1);
        testDevicePow(i) = devPower*sum(test_normalizado)/(60*1000);
    end
    disp(goalPow);
    disp(mean(testDevicePow));
    % comprobación de cumplimiento, permitiendo una tolerancia
    if abs(goalPow - mean(testDevicePow)) <= tol
        disp('Matrix found!')
        TRANS_MOD = TRANStry;
        break
    end
    
    if subir==true %viene aumentando el consumo
        if goalPow > mean(testDevicePow)
            subir=true; %sube consumo en la proxima iteración para alcanzar la meta
        else
            if goalPow < mean(testDevicePow)
                subir = false; %baja consumo en la proxima iteracion porque esta volado
                delta = 0.5*delta; % reduce el paso para acercarse mas lento
            end
        end        
    else %subir es false y viene reduciendo consumo
        if goalPow < mean(testDevicePow)
            subir=false;
        else
            if goalPow > mean(testDevicePow) % se pasó y redujo mucho
                subir = true;
                delta = 0.5*delta; % reduce el paso para acercarse mas lento
            end
        end
    end            
    
    if (t==maxTries)&&((goalPow - mean(testDevicePow)) > tol)
        disp('Sin convergencia');
        TRANS_MOD = TRANStry;
    end

end
