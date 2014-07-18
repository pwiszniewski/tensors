function A = rotate_set (A, k, l, c, s)
	 
	 ## usage: A = rotate_set (A, k, l, c, s)
	 ## 
	 ## apply a Jacobi rotation (c,s) to a set of matrices A. 
	 ##  
	 ##   A = J * A * J^h  

  nmat = size(A,3);
  matsz = size(A,1);
  
  for n = 1:nmat

    ## left multiply by R
      
    for j = 1:matsz
      [A(k,j,n), A(l,j,n)] = rotate(A(k,j,n),A(l,j,n),c,s);
    end
    
    ## right multiply by R^h
    
    for j = 1:matsz
      [A(j,k,n), A(j,l,n)] = rotate(A(j,k,n),A(j,l,n),conj(c),conj(s));
    end
    
  end

endfunction
