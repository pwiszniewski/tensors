function [T, U, lambda] = initrandcum (dimin, rnk, cplx)
	 
	 ## usage: [T, U, lambda] = initrandrec (dimin, rnk, cplx)
	 ## 
	 ## Builds 4-index tensor of a specified rank as a cumulant

	 T = zeros(dimin,dimin,dimin,dimin);
	 
	 lambda = rand(rnk,1);
#	 lambda = zeros(rnk,1);
#	 for j = 1:rnk
#	     lambda(j) = j;
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
		   T(j,k,l,m) += lambda(d) * U(j,d) * conj(U(k,d)) * conj(U(l,d)) * U(m,d);
		 end
	       end
	     end
	   end
	 end
	     
endfunction
