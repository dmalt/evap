script

In=Result(:,1);
Te=Result(:,3);

plot(In,Te,'LineWidth',3,'Color','r');
axis([-0.1 1.1 120. 170.]);
xlabel('In');
ylabel('Tw');
