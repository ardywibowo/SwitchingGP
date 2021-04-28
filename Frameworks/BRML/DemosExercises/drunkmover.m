close all
clear all
import brml.*

Gx=50; Gy=50;

T=100; % number of timesteps
P=500; % number of people

X{1}=zeros(Gy,Gx);X{2}=zeros(Gy,Gx);
for p=1:P
    x{p}(1)=randi([3 Gx-3]); y{p}(1)=randi([3 Gy-3]);
    if p==1
        x{p}(2)=x{p}(1)+2; y{p}(2)=y{p}(1)+2;
    else
        x{p}(2)=x{p}(1)+sign(rand-0.5); y{p}(2)=y{p}(1)+sign(rand-0.5);
    end
end

for p=P:-1:1
    for t=1:2
        if p==1
            X{t}(y{p}(t),x{p}(t))=2; % use =2 if you want to see the fast mover
        else
            X{t}(y{p}(t),x{p}(t))=1;
        end
        
    end
end


rn=0.99;
for t=3:T
    pause(0.5)
    X{t}=zeros(Gy,Gx);
    
    p=1; % person 1 is dangerous fast mover
    x{p}(t)=x{p}(t-1)+ 2*sign(rand-0.5);
    y{p}(t)=y{p}(t-1)+ 2*sign(rand-0.5);
    
    for p=2:P
        if rand<rn
            x{p}(t)=x{p}(t-1)+ sign(x{p}(t-1)-x{p}(t-2));
        else
            x{p}(t)=x{p}(t-1)- sign(x{p}(t-1)-x{p}(t-2));
        end
        
        if rand<rn
            y{p}(t)=y{p}(t-1)+ sign(y{p}(t-1)-y{p}(t-2));
        else
            y{p}(t)=y{p}(t-1)- sign(y{p}(t-1)-y{p}(t-2));
        end
    end
    for p=P:-1:1
        if x{p}(t)<=Gx && x{p}(t)>0 && y{p}(t)<=Gy && y{p}(t)>0 % the person is in the grid
            if p==1
                X{t}(y{p}(t),x{p}(t))=2;  % use =2 if you want to see the fast mover
            else
                X{t}(y{p}(t),x{p}(t))=1;
            end
        else % the person is outside the grid
            x{p}(t)=randi([3 Gx-3]); y{p}(t)=randi([3 Gy-3]);
            if p==1
                X{t}(y{p}(t),x{p}(t))=2;  % use =2 if you want to see the fast mover
            else
                X{t}(y{p}(t),x{p}(t))=1;
            end
        end
    end
    imagesc(X{t}); drawnow;
end

[x{1}(:) y{1}(:)] % list of the drunk's position

Xtrue=X;
for t=1:T
    X{t}=Xtrue{t}>0;
end
if 1==0
    save drunkproblemX X
    save drunkproblemXtrue Xtrue
end

