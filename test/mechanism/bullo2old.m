function [X]=bullo2old(s,N,A,theta,tMax,w, c_eps, varargin)
% c_ij: value the user i gives to j neighbours opinion.
% c_eps: error tolerance for reaching equilibrium.

NORM_TYPE = Inf;
wantPlot = false;

% Parse input
if (~isempty(varargin))
    for c=1:length(varargin)
        switch varargin{c}
            case {'plot'}
                wantPlot = true;
            otherwise
                error(['Invalid optional argument, ', varargin{c}]);
        end % switch
    end % for
end % if

%%var init
X = zeros(N,tMax+1);
X(:,1) = s;
%c_ij = 1;

theta(1,1)

%%Implementación del modelo combinado
for t=1:tMax
    for i=1:N
        vec=find(A(i,:));
        sel=randi(max(size(vec)),1,1);
        X(i, t+1)=(1-theta(i))*X(i,1)+theta(i)*w(i)*X(i,t)+theta(i)*(1-w(i))*X(vec(1,sel),t);
    end
     % Check if we have reached an equilibrium
    opinionDiff = diff(X(:,t:t+1),1,2);
    if (norm(opinionDiff, NORM_TYPE) < c_eps)
        disp(['[Bullo stubborness model] Reached equilibrium after ' num2str(t) ' rounds.']);
        break;
    end
end

X(:,t+2:end) = [];

  
    % Plot
if (wantPlot)
    plotOpinions2(X);
end

end

