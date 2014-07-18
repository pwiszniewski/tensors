function [ EI, U ] = initrandlrnk (dim, cplx, packing)
	 
	 ## usage: EI U = initrandhub (dim, cplx, packing)
	 ## 
	 ## 
	 ## packing = 1 Dirak
	 ##        <> 1 Mulliken 
	 ##
	 ## Initialize a random 4-index tensor EI of dimension 
         ## dim x dim x dim x dim and rank dim 

 EI_on_site = zeros(dim,dim,dim,dim);

 for j =  1 : dim
     EI_on_site(j,j,j,j) = rand(1,1);
     if (mod(j,2) != 0)
	EI_on_site(j,j,j,j) *= -1;
     endif
 end

 U = rand(dim) + i*rand(dim);
 U = orth(U);

 EI = t2e(EI_on_site, U, packing);

endfunction
