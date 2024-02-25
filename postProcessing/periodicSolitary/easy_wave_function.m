    
   x=linspace(0,20,100);
   y=cos(x/2.5);
   y1=0.6*cos(x/3+3);
   
   plot(x,y)
   
   hold on
   plot(x,y1)
   legend('reference','simulation');
   set(gca,'XTick',[])
   set(gca,'YTick',[])