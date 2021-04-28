function [x, t_x, y, t_y, isYFun, newx, new_tx, newy, new_ty, invalid] = pre(x, t_x, y, t_y, isNewSub)

   invalid=0;
   newx = []; new_tx= []; newy = 0; new_ty = [];
   isYFun = checkY(y);  %check whether y is a function
   if length(x)~=length(t_x)        % check whether x and t_x have the same length
       fprintf(1,'Error: invalid input, x and t_x should have equal length!\n');
       invalid=1;
       return;
   end       
   
   if isempty(isNewSub) || all(isNewSub == 0)
      if length(y)~=length(x)       % check whether y and x have the same length
          fprintf(1,'Error: invalid input, x and y should have equal length!\n');
          invalid=1;
          return;
      end  
      if isYFun == 1
          if length(t_y)~=length(y)
              fprintf(1,'Error: invalid input, y and t_y should have equal length!\n');
              invalid=1;
              return;
          end
      end
      newx = []; new_tx = []; newy = []; new_ty = [];
          
   elseif length(isNewSub) == 1 && isNewSub > 0            %when isNewSub is scalar, the last number 
      if isNewSub == length(x)
         fprintf(1,'Error: invalid isNewSub, at least one or more subjects are not for predicton!\n');
         invalid=1;
         return;
      end
      if length(y)>length(x)-isNewSub             % check length(y)
          fprintf(1,'Warning: length(y) should be length(x)-isNewSub. Reset it to be so now!\n');
          y = y(1:(length(x)-isNewSub));
      elseif length(y)<length(x)-isNewSub
          fprintf(1,'Error: invalid input of y. Please check length(y)!\n');
          invalid=1;
          return;
      end
      newx = x((end-isNewSub+1):end);                       %of subjects ("isNewSub") are used for prediction
	  new_tx = t_x((end-isNewSub+1):end);
      if isYFun == 1
          if length(t_y)~=length(t_x)
              fprintf(1,'Error: invalid input, t_y and t_x should have equal length!\n');
              invalid=1;
              return;
          end
          new_ty = t_y((end-isNewSub+1):end);
          t_y = t_y(1:(end-isNewSub)); 
      else
         new_ty = [];
      end
      x = x(1:(end-isNewSub));
      t_x = t_x(1:(end-isNewSub));              
      
 
   elseif length(isNewSub) == length(x)             %when isNewSub is an indicator,1 : new subject , 0: observed subject
     if all(isNewSub == 1)
        fprintf(1,'Error: invalid isNewSub, at least one or more subjects are not for predicton!\n');
        invalid=1;
        return;
     end
     if length(y)==length(x)
         fprintf(1,'Warning: length(y) should be length(x)-sum(isNewSub). Reset it to be so now!\n');
         y = y(isNewSub==0);
     elseif length(y)~=length(x)-sum(isNewSub)
         fprintf(1,'Error: invalid input of y. Please check length(y)!\n');
         invalid=1;
         return;
     end
     newx = x(isNewSub == 1);
     new_tx = t_x(isNewSub == 1);
      
      if isYFun == 1
          if length(t_y)~=length(t_x)
              fprintf(1,'Error: invalid input, t_y and t_x should have equal length!\n');
              invalid=1;
              return;
          end
          new_ty = t_y(isNewSub == 1);
          t_y = t_y(isNewSub == 0);
      else
         new_ty = [];
      end
      x = x(isNewSub == 0);
      t_x = t_x(isNewSub == 0);
      
 
   else
      fprintf(1,'Warning: isNewSub must be a positive integer or a 0-1 indicator vector!Reset isNewSub = [] now!\n');
      newx = []; new_tx = []; newy = []; new_ty = []; 
   end