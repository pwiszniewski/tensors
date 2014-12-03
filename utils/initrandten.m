function [T, U, lambda] = initrandten (dimin, rnk, isym, cplx)
	 
	 ## usage: [T, U, lambda] = initrandrec (dimin, rnk, isym, cplx)
	 ## 
	 ## Builds 4-index tensor of a specified rank with the 
         ## symmetry requested by isym
         ## dimin - dimension
	 ##   rnk - rank, rank =< dimin^2
	 ##  isym - symmetry 
         ##       = 0 cumulants
	 ##       = 1 Mulliken
	 ##       = 2 Dirak

	 T = zeros(dimin,dimin,dimin,dimin);
	 
#	 lambda = rand(rnk,1);
	 lambda = ones(rnk,1);
	 lambda = zeros(rnk,1);
	 for j = 1:rnk
	     lambda(j) = j;
	 end

	 if ( rnk > dimin )
	   U = rand(dimin,rnk); 
	   if (cplx == 1)
	     U += i*rand(dimin,rnk);
	   end

	   for j = 1:rnk
	       U(:,j) = U(:,j)/norm(U(:,j));
	   end

	 else
	   U = rand(dimin); 
	   if (cplx == 1)
	     U += i*rand(dimin);
	   end
	   U = orth(U);
	 end
	 
	 if ( isym == 0 )
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
	   
	 elseif ( isym == 1 ) 
	   for j = 1:dimin
	     for k = 1:dimin
	       for l = 1:dimin
		 for m = 1:dimin
		   for d = 1:rnk
		     T(j,k,l,m) += lambda(d) * conj(U(j,d)) * U(k,d) * conj(U(l,d)) * U(m,d);
		   end
		 end
	       end
	     end
	   end
	   
	 elseif (isym == 2)
	   for j = 1:dimin
	     for k = 1:dimin
	       for l = 1:dimin
		 for m = 1:dimin
		   for d = 1:rnk
		     T(j,k,l,m) += lambda(d) * conj(U(j,d)) * conj(U(k,d)) * U(l,d) * U(m,d);
		   end
		 end
	       end
	     end
	   end
		
	 end
endfunction
