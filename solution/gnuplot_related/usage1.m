% example illustrating USAGE 1 (see the readme file)
b=gpt('foo.plt');
x1=b.set0.x;y1=b.set0.y;
x2=b.set1.x;y2=b.set1.y;
x3=b.set2.x;y3=b.set2.y;
plt=plot(x1,y1,'b-',x2,y2,'k-',x3,y3,'r-');
set(plt(2),'linewidth',2)
a=axis;
axis([a(1) a(2) -1.2 1.2]);
legend('cos(\theta)','sin(\theta)/\theta','\theta');
xlabel('\theta  \rightarrow');
ylabel('f(\theta)  \rightarrow');
title('DEMO PLOT FOR GPT (USAGE 1)');