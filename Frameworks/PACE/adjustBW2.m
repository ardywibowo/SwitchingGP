function [bw_xcov] = adjustBW2(kernel,bw_xcov, npoly, nder,regular,verbose)

  if strcmp(kernel,'gauss')
     if regular == 2
         bwxcov_fac = [1.1 0.8 0.8];
     else
         bwxcov_fac = [1.1 1.2 2];
     end
     if nder > 2
        facID = 3;
     elseif nder >= 0 && nder <= 2
        facID = nder+1;
     else
        facID = 1;
     end
     bw_xcov = bw_xcov.*bwxcov_fac(facID);
  end

  if strcmp(kernel,'epan')
     if regular == 2
         bwxcov_fac = [1.1 1.0 1.1];
     else
         bwxcov_fac = [1.1 1.2 1.5];
     end
     if nder > 2
        facID = 3;
     elseif nder >= 0 && nder <= 2
        facID = nder+1;
     else
        facID = 1;
     end
     bw_xcov = bw_xcov.*bwxcov_fac(facID);
  end

end
