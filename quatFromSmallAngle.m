function [ quatr ] = quatFromSmallAngle( e )
%QUATFROMSMALLANGLE Summary of this function goes here
%   Detailed explanation goes here
quatr=[1 0 0 0];
e1=e(1);
e2=e(2);
e3=e(3);
q_squared=e1*e1+e2*e2+e3*e3;
if(q_squared<1)
   quatr(1)=(1-q_squared)^0.5;
   quatr(2)=e1*0.5;
   quatr(3)=e2*0.5;
   quatr(4)=e3*0.5;
else
    w=1.0/((1+q_squared)^0.5);
    f=w*0.5;
    quatr(1)=w;
    quatr(2)=e1*f;
    quatr(3)=e2*f;
    quatr(4)=e3*f;
end
quatr=quatnormalize(quatr);

end

