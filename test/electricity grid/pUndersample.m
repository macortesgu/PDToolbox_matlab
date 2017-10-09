function Prefs = pUndersample(PrefsOrig,factor,strategies)

T_ = size(PrefsOrig,2)/factor;
Tperf = size(PrefsOrig,2);
N = size(PrefsOrig,1);
P = size(PrefsOrig,3);
Prefs = zeros(N,T_,P);

if(mod(Tperf,factor)~=0)
    error('La cantidad de periodos seleccionada no es m�ltiplo entero del perfil original');
else
    mean_ocurr = zeros(N,P);
    for dev = 1:N
        for u = 1:P
            mean_ocurr(dev,u) = length(nonzeros(PrefsOrig(dev,:,u)))/Tperf;
        end
    end
    
    occurrences = zeros(N,T_,P);
    for u = 1:P
        for dev = 1:N
            for t = 0:(T_-1)
                occurrences(dev,t+1,u) = length(nonzeros(PrefsOrig(dev,(t*factor+1):(t*factor+factor),u)));
                Prefs(dev,t+1,u) = (round(occurrences(dev,t+1,u)/factor * strategies) / strategies);               
            end
        end
    end
end
end