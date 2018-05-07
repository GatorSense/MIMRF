/* Find row index if vector A is part of a row in Vector B
 * Used in the function [aaa,bbb,ccc] = ismember_findrow_mex_my(A,B)
 *Notes: This code expects two sorted arrays!!! Has to be asending order in both arrays
 * Can deal with repetative arrays but only output the first element that is a member
 * INPUT
 *      - Vector A 
 *      - Vector B 
 * OUTPUT
 *      -ans: logical
 *      -outIndex: vector of index in the second vector
 *      -count: number of index in the second vector that is a member
 *
 *
 *  Example:
 *      Input A=[3 4];B=[1 2 3;1 2 4;1 3 4; 2 3 4];
 *      Output:
 *          aaa = [1,1] (logical); 
 *          bbb = [2,3,0,...] (12x1 vector, the first 3&4 element in B); 
 *          ccc=[3;4] (THIS IS THE INDEX FOR ROW NUMBER!)
 *
 *
 **  Written by: X. Du 09/09/2015
*
* This product is Copyright (c) 2016 University of Missouri
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
* 3. Neither the name of the University nor the names of its contributors
* may be used to endorse or promote products derived from this software
* without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF MISSOURI AND
* CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
* OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
 */

#include <matrix.h>
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *in1,*in2full;
    mxLogical *out;
    int N,N1,N2,count,i,j,buflen,buflen2,ii,jj,Nout,inout;
    double *loc;
    double *outIndex;
    double *outcount;
    double tempb;
    double *in3;
    double NumRows;
    int NumRowsInt,NumColsInt;
    double *in1p, *in2p;
    double *in2;  
    
    in1 = mxGetPr(prhs[0]);
    in2full = mxGetPr(prhs[1]);
    in3 =  mxGetPr(prhs[2]);
            
    N = (int)mxGetNumberOfElements(prhs[0]); 
    N1= (int)mxGetNumberOfElements(prhs[0]);  // size of A vector
    N2 = (int)mxGetNumberOfElements(prhs[1]); // size of B vector
    Nout = (int)mxGetNumberOfElements(prhs[0]); 
   
 
    NumRows = in3[0];   
    NumRowsInt = (int)NumRows;  
    NumColsInt = N2/NumRowsInt;
   
            
    plhs[0] = mxCreateLogicalMatrix( N,1 );
    out = mxGetLogicals( plhs[0] );
    
    plhs[1] = mxCreateDoubleMatrix(N2,1,mxREAL);
    outIndex = mxGetPr( plhs[1] );
    
    plhs[2] = mxCreateDoubleMatrix(NumColsInt,1,mxREAL);
    outcount = mxGetPr( plhs[2] );
  
     // double in2[NumRowsInt];

        for (j=0;j<NumColsInt;j++) {
              //  in2 = mxMalloc(NumRowsInt+1);
  buflen = mxGetNumberOfElements(prhs[1]) + 1;
  buflen2 = mxGetNumberOfElements(prhs[1]) + mxGetNumberOfElements(prhs[1]) + 2;
  in2 = mxCalloc(buflen, sizeof(double));
  in2p = mxCalloc(buflen2, sizeof(double));
  
                for(i=0;i<NumRowsInt;i++) { //for loop terminates if count>n
        in2[i] = in2full[i+(j*NumRowsInt)]; 
        in2p[i]=in2full[i+(j*NumRowsInt)];
        }
   
  
      in1p = mxGetPr(prhs[0]);
    count = 0;
    N = (int)mxGetNumberOfElements(prhs[0]); 
    jj=0;
    ii=0;
    for( N+=NumRowsInt; N>0; N--){
    
        while((jj<NumRows)&&(ii<N1)){
            
        if( in1p[ii] <= in2p[jj] ){ //is or is not member
              out[ii]=(in1p[ii] == in2p[jj]); 
           if(in1p[ii] == in2p[jj]){     
                outIndex[count] = jj+1;
                count += 1;

            }     
           
                ii= ii+1;
               
        }
        
        else{  //keep looking     
            jj= jj+1;
        }
                         
    }}  
    
    
    
    
    
    
    if (count==N1){
    *outcount++ = j+1; 
    }
        
    mxFree(in2);
    mxFree(in2p);
    }
   
}