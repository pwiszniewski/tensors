function EI = dirtocum (ET)
	 
	 ## usage: EI = dirtocum (ET)
	 ## Reorders Dirak integrals to cumulants 
         ## Note that this is not correct in general, as 
	 ## only complex cumulants have the same number 
         ## of symmetries as real Dirak integrals
	 
  dim = size(ET,1);
  EI = zeros(dim,dim,dim,dim);


  ## repack ij|lk = (ij|kl)

    for i = 1 : dim
       for j = 1 : dim
	for k = 1 : dim
	  for l = 1 : dim
	    EI(i,j,l,k) = ET(i,j,k,l);		
	  end
	end
      end
    end
    
 
endfunction
