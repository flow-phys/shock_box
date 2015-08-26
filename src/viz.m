function viz(plotter)
global rho u v p x_c y_c Dxy nA nB probname perx pery

%  Legend for plot types
%  1 - 1D density only plot
%  2 - 1D all variables plot via subplot
%  3 - 2D Psuedo color plot of pressure and density 
%  4 - 2D Contour color plot of pressure and density 
%  5 - 2D Psuedo color plot of mag(U) and 2D Contour color plot of pressure
%  6 - 2D Velocity Vector plot

% switch probname 
%     case{'advect1d'}
%         plotter(1) = 1;           
%              
%     case{'advect2d_cart'}
%         plotter(1) = 4;
%                 
%     case{'advect2d_wavy'}
%         plotter(1) = 4;
%         
%     case{'shocktube'}
%         plotter(1) = 2;
%         
%     case{'cylinder'}
%         plotter(1) = 5;
%             
%     case{'bullet'}
%         plotter(1) = 5;
%         
%     case{'gauss2d_cart'}
%         plotter(1) = 5;
%     
%     case{'gauss2d_wavy'}
%         plotter(1) = 4;
%         
%     case{'nozzle'}
%         plotter(1) = 5;
%     
%     case{'ramp'}
%         plotter(1) = 5;
% end

for i=1: max(size(plotter))
    switch plotter(i)
        % 1D density only plot
        case{1}
            figure(i);
            plot(x_c(:,:), rho(:,:) );drawnow;
       
        % 1D all variables plot via subplot
        case{2}                
            figure(i);
            subplot(3,1,1); plot(x_c(:,:), p(:,:) ,'*');drawnow;
            subplot(3,1,2); plot(x_c(:,:), rho(:,:) ,'*');drawnow;
            subplot(3,1,3); plot(x_c(:,:), u(:,:) ,'*');drawnow;
            
        % 2D Psuedo color plot of pressure and density 
        case{3}
            figure(i);
            subplot(2,1,1);
            surf(x_c,y_c,p,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');view(0,90);axis image;drawnow;
            subplot(2,1,2);
            surf(x_c,y_c,rho,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');view(0,90);axis image;drawnow;
            
        % 2D Contour color plot of pressure and density 
        case{4}
            figure(i);
            subplot(1,2,1);
            contour(x_c,y_c,p,30);axis image;drawnow;
            subplot(1,2,2);
            contour(x_c,y_c,rho,30);axis image;drawnow;
            
        % 2D Psuedo color plot of mag(U) and 2D Contour color plot of
        % pressure
        case{5}
            mag = sqrt ( u.*u + v.*v );
            %[tmp,d1] = grad(u,perx,pery);
            %[d2,tmp] = grad(v,perx,pery);
            %vort = d1 + d2;
            figure(i);
            subplot(1,2,1);
            surf(x_c,y_c,mag,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');view(0,90);axis image;drawnow;
            subplot(1,2,2);
            contour(x_c,y_c,p,30);axis image;drawnow;
            
        %  2D Vector plot of velocity (quiver)
        case{6}
            figure(i);
            quiver(x_c,y_c,u,v);drawnow;
            
    end
end
