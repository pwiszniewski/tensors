function [c, s] = find_angle (h, minmax)
	 
	 ## usage: [c, s] = find_angle_max (h, minmax)
	 ## 
	 ## h - row vector
	 ## minmax - 0 - diagonalize
         ##          1 - zero-diagonalize 
	 ##
	 ## Finds the angle which diagonalizes or zero-diagonalizes
	 ## a general set of matrices.
	 ## h_k is a 1x3 vector  defined as 
	 ## h_k = [ A^k_ii  - A^k_jj, A_ij, A_ji ]
	 ## The factor G is sum_k [ h_k^H * h_k ] 

 C = [ 1 0 0; 0 1 -i; 0 1 i ];
 G = real(C' * h' * h * C);

 [v d] = eig(G);
 [dump, perm] = sort(diag(d));

 if ( minmax == 0 )
   v = v(:,perm(3));
 else
   v = v(:,perm(1));
 endif

 if (v(1) < 0 )
    v = -v;
 endif
 
 r = norm(v);

 c = sqrt( ( v(1) + r ) / 2);
 s = ( v(2) - i*v(3) ) / sqrt( 2 * r * ( v(1) + r ) );

endfunction
