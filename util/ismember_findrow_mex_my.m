function [aaa,bbb,ccc] = ismember_findrow_mex_my(A,B)
% Find the index of which row in matrix B fully contains A vector
% Input:
%  Example: A=[3 4];B=[1 2 3;1 2 4;1 3 4; 2 3 4];
% Output:
% Example: aaa = [1,1] (logical); bbb = [2,3,0,...] (12x1 vector, the first
% 3&4 element in B); ccc=[3;4] (THIS IS THE INDEX FOR ROW NUMBER!)
% Uses the "ismember_findrow_mex.c" code. 
% Written by: X. Du 09/10/2015

input2 =reshape(B',[size(B,1)*size(B,2),1]);
numrow = size(B,2);
[aaa,bbb,ccc]=ismember_findrow_mex(A',input2,numrow);
ccc(ccc==0)=[];

end