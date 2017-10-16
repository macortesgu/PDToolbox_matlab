function [Aw] = addWeights(A)

%parse input

if(size(A,1)~=size(A,2))
    error('Invalid non-square adjacency matrix of dim: %dx%d',size(A,1),size(A,2));
end

Aw = zeros(length(A));    

for i = 1:size(A,2)
    
    col_indexes = find(A(i,:)); %take row by row
    
    c_acum = 0;
    a = 0;
    
    num_neighbours = size(col_indexes,2);
    weights = zeros(1,num_neighbours);
   
    for j = 1:num_neighbours
        
        if (j == num_neighbours)
             weights(j) = 1 - sum(weights(1:num_neighbours-1));
             break;
        end
        
        b = 1 - c_acum;
        weights(j) = a + (b-a)*rand(1,1); 
        c_acum = c_acum + weights(j);         
        
    end    
    Aw(i, col_indexes) = weights;
        
end
Aw = sparse(Aw);
if(size(Aw,1)~=size(Aw,2))
    error('Corrupter output. Invalid non-square adjacency matrix of dim: %dx%d',size(A,1),size(A,2));
end
end