function [ omegaMatrix ] = omegaMatJPL( e )
%OMEGAMATJPL Summary of this function goes here
%   Detailed explanation goes here
   
    
    e0=e(1);
    e1=e(2);
    e2=e(3);
     
    omegaMatrix=[0,e2,-e1,e0;-e2,0,e0,e1;e1,-e0,0,e2;-e0,-e1,-e2,0];
end


