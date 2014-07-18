function EI = csehub (dim, cplx)
	 
	 ## usage: ei = csehub (dim)
	 ## 
         ## build CSE Hamiltonian for 
         ## Hubbard model (mulliken)
         ## This function needs to be revised!!!!!!!!!!!!!!!!!
  
  t = 1;
  dopbc = 1;

  
  h = oneeh(dim, t, dopbc);
  EI = initrandhub(dim, cplx, 0);
  h = h ./ (dim - 1);

  for i = 1 : dim
    for j = 1 : dim
      for k = 1 : dim 
	for l = 1 : dim

	  if (k == l)
	    EI(i,j,k,l) += h(i,j);
	  endif
	  
	  if (i == j)
	    E(i,j,k,l) += h(k,l);
	  endif

	end
      end
    end
  end
  
  
endfunction
