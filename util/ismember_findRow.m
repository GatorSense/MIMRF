function [row_id] = ismember_findRow(A,B)
%%% This function returns row number if A is a row of B
%%% Same as [~,row_id] = ismember(A,B,'rows'), may be faster.
%%% E.g., A=[3 4]; B=[1 2; 3 4; 5 6]; THEN row_id=2
%%% Written by: X. Du 08/24/2015

Atemp = repmat(A,[size(B,1),1]);
diff = Atemp-B;
alltemp = all(diff==0,2); %find rows that has all zero difference values
row_id = find(alltemp==1);

end