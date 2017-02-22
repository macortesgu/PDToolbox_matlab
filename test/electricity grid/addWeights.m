function [Aw] = addWeights(A)

%parse input

if(size(A,1)~=size(A,2)
    error(['Invalid non-square adjacency matrix of dim: %dx%d'],size(A,1),size(A,2));
end

    

for i = 1:size(A,2)
    
    col_indexes = find(A(i,:); %take row by row
    
    

end
