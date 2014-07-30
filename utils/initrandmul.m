function [T, U, lambda] = initrandmul (dimin, rnk, cplx)
	 
	 ## usage: [T, U, lambda] = initrandrec (dimin, rnk, cplx)
	 ## 
	 ## Builds 4-index tensor of a specified rank in Mulliken 
	 ## notation

	 T = zeros(dimin,dimin,dimin,dimin);
	 
	 lambda = rand(rnk,1);

#	 if (mod(rank,2) == 0 )
#	   lai = i*rand(rank/2,1);
#	   lambda(1:rank/2)    = lambda(1:rank/2) + lai;
#	   lambda(rank/2+1:rank) = conj(lambda(1:rank/2));
#	 else
#	   lai = i*rand((rank-1)/2,1);
#	   lambda(1:(rank-1)/2)      = lambda(1:(rank-1)/2) + lai;
#	   lambda((rank-1)/2+1:rank-1) = conj(lambda(1:(rank-1)/2));
     
#	 end

	 if ( rnk >= dimin )
	   U = rand(rnk,dimin); 
	   if (cplx == 1)
	     U += i*rand(rnk,dimin);
	   end
	   U = orth(U);
	   U = U';
	 else
	   U = rand(dimin); 
	   if (cplx == 1)
	     U += i*rand(dimin);
	   end
	   U = orth(U);
	   U = U(:,1:rnk);
	 end

	 for j = 1:dimin
	   for k = 1:dimin
	     for l = 1:dimin
	       for m = 1:dimin
		 for d = 1:rnk
		   T(j,k,l,m) += lambda(d) * conj(U(j,d)) * U(k,d) * \
				 conj(U(l,d)) * U(m,d);
		 end
	       end
	     end
	   end
	 end
	     
	     
endfunction
