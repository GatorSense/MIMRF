function   [indexRowCol] = rowcol(indexlist,Xim)

% USAGE 
%    This function yields the row and col number given index number in Xim
% WRITTEN BY 
%    X. Du 12/17/2013

indexRowCol(:,1) =  mod(indexlist,size(Xim,1));
indexRowCol(:,2) = min( size(Xim,2), floor(indexlist/size(Xim,1))+1);
 for j = 1:size(indexRowCol(:,1))
   if indexRowCol(j,1)==0
       indexRowCol(j,2) = min( size(Xim,2), floor(indexlist(j)/size(Xim,1)));
   end
       
end

 indexRowCol(find(indexRowCol(:,1)==0))=size(Xim,1);
end