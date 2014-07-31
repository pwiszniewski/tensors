function [ T ] = buildcum (U, lambda)
	 
	 ## usage: [ T ] = buildcum (U, lambda)
	 ## 
	 ## Builds 4-index tensor of a specified rank as a cumulant

	 dimin = size(U,1);
	 rnk   = size(U,2);

	 if (rnk != length(lambda))
	    error ('incompartible size\n');
	 end

	 T = zeros(dimin,dimin,dimin,dimin);
	 
	 for j = 1:dimin
	   for k = 1:dimin
	     for l = 1:dimin
	       for m = 1:dimin
		 for d = 1:rnk
		   T(j,k,l,m) += lambda(d) * U(j,d) * conj(U(k,d)) * \
				 conj(U(l,d)) * U(m,d);
		 end
	       end
	     end
	   end
	 end
	     
	     
endfunction
