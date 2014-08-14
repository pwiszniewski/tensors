function [ T ] = buildten (U, lambda, isym)
	 
	 ## usage: [ T ] = buildten (U, lambda, isym)
	 ## 
	 ## Builds 4-index tensor of a specified rank and symmetry
	 ##
	 ## isym = 0 Cumulants
	 ##      = 1 Mulliken
	 ##      = 2 Dirak

	 dimin = size(U,1);
	 rnk   = size(U,2);

	 if (rnk != length(lambda))
	    error ('incompartible size\n');
	 end

	 T = zeros(dimin,dimin,dimin,dimin);
	 
	 if ( isym == 0 )
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
	   
	 elseif ( isym == 1 )     
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

	 elseif ( isym == 2 )
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
	   
	 end		

endfunction
