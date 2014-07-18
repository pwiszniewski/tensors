function EI = repack (ET)
	 
	 ## usage: EI = repack (ET)
	 ## symmetrizes integrals in the Mulliken or Dirak ordering 
         ## 
	 
  dim = size(ET,1);
  EI = zeros(dim,dim,dim,dim);


  ## repack <ik|jl> = (ij|kl)

    for i = 1 : dim
       for j = 1 : dim
	for k = 1 : dim
	  for l = 1 : dim
	    EI(i,k,j,l) = ET(i,j,k,l);		
	  end
	end
      end
    end
    
 
endfunction
