function [T, U, lambda] = initranddir (dimin, rnk, cplx)
	 
	 ## usage: [T, U, lambda] = initrandrec (dimin, rnk, cplx)
	 ## 
	 ## Builds 4-index tensor of a specified rank in Mulliken 
	 ## notation

	 T = zeros(dimin,dimin,dimin,dimin);
	 
	 lambda = rand(rnk,1);

	 if ( rnk >= dimin )
	   U = rand(dimin,rnk); 
	   if (cplx == 1)
	     U += i*rand(dimin,rnk);
	   end
	 else
	   U = rand(dimin); 
	   if (cplx == 1)
	     U += i*rand(dimin);
	   end
	   U = orth(U);
	 end

	 for j = 1:dimin
	   for k = 1:dimin
	     for l = 1:dimin
	       for m = 1:dimin
		 for d = 1:rnk
		   T(j,k,l,m) += lambda(d) * conj(U(j,d)) * conj(U(k,d)) * \
				 U(l,d) * U(m,d);
		 end
	       end
	     end
	   end
	 end
	     
	     
endfunction
