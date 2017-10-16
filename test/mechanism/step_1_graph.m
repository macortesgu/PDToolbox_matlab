function Aw = step_1_graph(P)

%%graph part
%% opinion  dynamics parameters

netType = 'WS';
K = 4;

switch netType
    
    case 'WS'
        % Watts-Strogatz Graph Model Parameters
        beta = 0.6; %link probability
        Graph = WattsStrogatz(P,K,beta); %%graph with small world properties
        A = adjacency(Graph);
        % % %
    case 'BA'
        %Barabasi-Albert Graph Model Parameters
        m = 2*K;% Same K as in WS model
        mo = m+1; %seed
        A = sparse(BAgraph_dir(P,mo,m));
        % % %
    otherwise
        % Watts-Strogatz Graph Model Parameters
        beta = 0.6; %link probability
        Graph = WattsStrogatz(P,K,beta); %%graph with small world properties
        A = adjacency(Graph);
        % % %
end


Aw = addWeights(A); %add weights to the adjacency matrix, weiths are agents valuation of neighbours'respective opinion.

%% Plot weighted graph

Gw = digraph(Aw);
figure(20)
pp = plot(Gw);
colormap(spring)
layout(pp,'auto')
%view(3)

Gw.Nodes.NodeColors = indegree(Gw); %or outdegree
pp.NodeCData = Gw.Nodes.NodeColors;
colorbar

Gw.Edges.LWidths = 7*Gw.Edges.Weight/max(Gw.Edges.Weight);
pp.LineWidth = Gw.Edges.LWidths;

name_Graph = strcat('graph',netType,num2str(P));

%savefig(1, name_Graph, 'compact');


% activefig = figure(1);
% activefig.PaperOrientation = 'landscape';
% activefig.PaperType = 'a3';
% print -painters -dpdf -r600 -bestfit opinionEvolution-N60.pdf

%print -painters -dpdf -r600 -bestfit WS-20.pdf