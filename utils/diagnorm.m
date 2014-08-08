function [ normd normof ] = diagnorm (EI, isym )
	 
	 ## usage: [ normd normof ] = bidiagnorm (EI)
	 ## 
	 ## 
	 
  dim = size(EI,1);
  
  normf = norm(vec(EI),'fro');

  normd = 0;

  ## Any ordering

  for j = 1:dim
      normd += EI(j,j,j,j)*conj(EI(j,j,j,j));
  end

  normof = normf^2 - normd;

endfunction
