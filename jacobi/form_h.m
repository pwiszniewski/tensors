function h = form_h (A, k, l)
	 
	 ## usage: h = form_h (A, k, l)
	 ## 
	 ##        (k,l) - indices of the maximal element 
	 ## form row vector h used in joint diagonalization
	 ## A is an array of matrices indexed by the last index.
	 ##

  
  nmat = size(A,3);
  h = zeros(nmat,3);

  for j = 1:nmat
      h(j,:)= [A(k,k,j) - A(l,l,j), A(k,l,j), A(l,k,j)] ;
  end

endfunction
