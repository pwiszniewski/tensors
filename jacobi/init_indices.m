function [ indB B ] = init_indices (A)
	 
	 ## usage: [ indx B ]= init_indices (A)
	 ## 
	 ## initialize vector indx and a matrix B
	 ## B(k,l) contains sums of squares 
	 ## of an element (k,l) in the set of matrices A
	 ##
	 ## indx are column indices of the maximal elements in each
         ## row of B.
	 
	 
  matsz = size(A,1);
  nmat = size(A,2)/matsz;

  B = zeros(matsz);

  for j=1:matsz
    for k=1:matsz
      indA = [k:matsz:nmat*matsz];
      B(j,k) = norm(A(j,indA),'fro')^2;
    end
  end
  
  indB = zeros(1,matsz);
  
  for j = 1:matsz
    indB(j) = maxcol(B, j);
  end

endfunction
