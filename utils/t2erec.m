function EI_U = t2erec (EI, U)
	 
	 ## usage: EI_U = t2erec (EI, U)
	 ## 
	 ## transform Dirak ordered integrals with rectangular 
	 ## matrices

  dim  = size(U,2);
  rank = size(U,1);

  TMP = zeros(rank,rank,rank,rank);

  for j = 1 : dim
    for k = 1 : dim
      for l = 1 : dim
	for p = 1 : rank
	  for i = 1 : dim
	    TMP(p,j,k,l) += EI(i,j,k,l) * U(p,i);
	  end
	end
      end
    end
  end
    
  EI_U = zeros(rank,rank,rank,rank);
  
  for p = 1 : rank
    for k = 1 : dim
      for l = 1 : dim
	for q = 1 : rank
	  for j = 1 : dim
	    EI_U(p,q,k,l) += TMP(p,j,k,l) * U(q,j);
	  end
	end
      end
    end
  end
    
  TMP = zeros(rank,rank,rank,rank);
  
  for p = 1 : rank
    for q = 1 : rank
      for l = 1 : dim
	for r = 1 : rank
	  for k = 1 : dim
	    TMP(p,q,r,l) += EI_U(p,q,k,l) * conj(U(r,k));
	  end
	end
      end
    end
  end
  
  EI_U = zeros(rank,rank,rank,rank);
  
  for p = 1 : rank
    for q = 1 : rank
      for r = 1 : rank
	for s = 1 : rank
	  for l = 1 : dim
	    EI_U(p,q,r,s) += TMP(i,q,r,l) * conj(U(s,l));
	  end
	end
      end
    end
  end
  
  endfunction
  
