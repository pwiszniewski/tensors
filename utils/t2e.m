function EI_U = t2e (EI, U, packing)
	 
	 ## usage: EI_U = t2e (EI, U, packing)
	 ## 
	 ## packing = 1 Dirak
	 ##        <> 1 Mulliken 
	 ##
	 ## Transform a 4-index tensor EI of dimension 
         ## dim x dim x dim x dim to new basis given by U

  dim  = size(U,1);

  if ( packing == 0 )

    TMP = zeros(dim,dim,dim,dim);

    for i = 1 : dim
      for j = 1 : dim
	for k = 1 : dim
	  for s = 1 : dim
	    for l = 1 : dim
	      TMP(i,j,k,s) += EI(i,j,k,l) * U(l,s);
	    end
	  end
	end
      end
    end
    
    EI_U = zeros(dim,dim,dim,dim);
    
    for i = 1 : dim
      for j = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for k = 1 : dim
	      EI_U(i,j,r,s) += TMP(i,j,k,s) * conj(U(k,r));
	    end
	  end
	end
      end
    end
    
    TMP = zeros(dim,dim,dim,dim);
    
    for i = 1 : dim
      for q = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for j = 1 : dim
	      TMP(i,q,r,s) += EI_U(i,j,r,s) * U(j,q);
	    end
	  end
	end
      end
    end
    
    EI_U = zeros(dim,dim,dim,dim);
    
    for p = 1 : dim
      for q = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for i = 1 : dim
	      EI_U(p,q,r,s) += TMP(i,q,r,s) * conj(U(i,p));
	    end
	  end
	end
      end
    end
  
  elseif (packing == 1)
      
    TMP = zeros(dim,dim,dim,dim);

    for i = 1 : dim
      for j = 1 : dim
	for k = 1 : dim
	  for s = 1 : dim
	    for l = 1 : dim
	      TMP(i,j,k,s) += EI(i,j,k,l) * U(l,s);
	    end
	  end
	end
      end
    end
    
    EI_U = zeros(dim,dim,dim,dim);
    
    for i = 1 : dim
      for j = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for k = 1 : dim
	      EI_U(i,j,r,s) += TMP(i,j,k,s) * U(k,r);
	    end
	  end
	end
      end
    end
    
    TMP = zeros(dim,dim,dim,dim);
    
    for i = 1 : dim
      for q = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for j = 1 : dim
	      TMP(i,q,r,s) += EI_U(i,j,r,s) * conj(U(j,q));
	    end
	  end
	end
      end
    end
    
    EI_U = zeros(dim,dim,dim,dim);
    
    for p = 1 : dim
      for q = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for i = 1 : dim
	      EI_U(p,q,r,s) += TMP(i,q,r,s) * conj(U(i,p));
	    end
	  end
	end
      end
    end
    
  elseif (packing == 2)

    TMP = zeros(dim,dim,dim,dim);

    for i = 1 : dim
      for j = 1 : dim
	for k = 1 : dim
	  for s = 1 : dim
	    for l = 1 : dim
	      TMP(i,j,k,s) += EI(i,j,k,l) * U(l,s);
	    end
	  end
	end
      end
    end
    
    EI_U = zeros(dim,dim,dim,dim);
    
    for i = 1 : dim
      for j = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for k = 1 : dim
	      EI_U(i,j,r,s) += TMP(i,j,k,s) * conj(U(k,r));
	    end
	  end
	end
      end
    end
    
    TMP = zeros(dim,dim,dim,dim);
    
    for i = 1 : dim
      for q = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for j = 1 : dim
	      TMP(i,q,r,s) += EI_U(i,j,r,s) * conj(U(j,q));
	    end
	  end
	end
      end
    end
    
    EI_U = zeros(dim,dim,dim,dim);
    
    for p = 1 : dim
      for q = 1 : dim
	for r = 1 : dim
	  for s = 1 : dim
	    for i = 1 : dim
	      EI_U(p,q,r,s) += TMP(i,q,r,s) * U(i,p);
	    end
	  end
	end
      end
    end

  endif

  endfunction
  
