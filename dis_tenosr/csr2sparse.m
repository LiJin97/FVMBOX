function A = csr2sparse(val, row_ptr, col_ind) 

 

 
 % Determine indexing base. 
 
 
 % Compute full row index. 
 m = length(row_ptr) - 1; 
 row_ind = zeros(length(col_ind), 1); 
 for i = 1:m 
   row_ind(row_ptr(i)+1:row_ptr(i+1)) = i; 
 end 
 
 
 % Construct matrix. 
% [row_ind,col_ind,]
 A = sparse(row_ind, col_ind , val); 
 
 
 end 
