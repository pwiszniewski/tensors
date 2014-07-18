function [ EI, U ] = initrandhub (dim, cplx, packing)
	 
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
     EI_on_site(j,j,j,j) = j;
 end

 U = rand(dim);
 if (cplx == 1)
   U += i*rand(dim);
 end
 
 U = orth(U);
     
 EI = t2e(EI_on_site, U, packing);

endfunction
