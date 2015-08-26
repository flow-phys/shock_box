clear all;
clc;

load density;

line = 0;

 figure(2);clf;

 xmax = max(max(x_c));
 ymax = max(max(max(RHO)));
 ymin = min(min(min(RHO(:,:,2:end))));

 for tt = 2: size(RHO,3)
     
 [nx,ny,time] = size(RHO);    
 halfx = nx/2 + 1;
 
     if (line == 1)
     hold off;
     plot(x_c(halfx,:),RHO(halfx,:,tt),'b*');
     hold on;
     plot(x_c(1,:),RHO(1,:,tt),'b*');
     axis([-xmax xmax ymin ymax]);
     drawnow;
     pause(.05);
     end
         
         if (line==0)
            
            surf(x_c,y_c,RHO(:,:,tt));
            axis([-xmax xmax -xmax xmax ymin ymax]);
            %view(0,90);
            view(45,60);
            drawnow;
            pause
         end 
 end