function [powx,powy,powz]=powmat

   vecx=[0 1 0 1 2 0];    vecy=[0 0 1 1 0 2];
   for ip=1:6
   powx(ip,:)=vecx+vecx(ip);
   powy(ip,:)=vecy+vecy(ip);
   end

