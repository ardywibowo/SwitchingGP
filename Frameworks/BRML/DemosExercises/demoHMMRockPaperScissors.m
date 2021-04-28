function demoHMMRockPaperScissors(showfiltered)
%DEMOHMMROCKPAPERSCISSORS another HMM inference demo
close all
import brml.*
H = 4; % number of Hidden states
V = 3; % number of Visible states
T = 20; % length of the time-series

[rock paper scissors]=assign(1:3);
s{1}='rock'; s{2}='paper'; s{3}='scissors';

c(1)=randi(3); % random first computer move
for currenttime=1:T
    
    done=false;
    while ~done
        in=input(['Game ',num2str(currenttime), ': Your Go! Rock(R), Paper (P) or Scissors (S)'],'s');
        switch lower(in)
            case 'r'
                o(currenttime)=rock; done=true;
            case 'p'
                o(currenttime)=paper; done=true;
            case 's'
                o(currenttime)=scissors; done=true;
            otherwise
                disp('false input: must be r, p, or s. Try again.')
                done=false;
        end
    end
    
    disp(['you play ' s{o(currenttime)},  ' , computer plays ' s{c(currenttime)}]);
    
    w(currenttime)=RPSwinner(o(currenttime),c(currenttime));
    switch w(currenttime)
        case 1
            st='you win!';
        case 0
            st='draw!';
        case -1
            st='you lose!';
    end
    disp([st ' You won ' num2str(100*mean(w>0)) ' percent of the games up to now. Computer won '  num2str(100*mean(w<0)) ' percent'])
    
    % Do filtering in the HMM to estimate the strategy the human is playing
     ph1 = condp(ones(H,1)); % could be any strategy at timestep 1
    
    for t=1:(currenttime+1)
        phghm(:,:,t) = condp(ones(H,H)+15*eye(H));% transition distribution p(h(t)|h(t-1)). Tend to stay in the same strategy
        pvgh(:,:,t)=zeros(V,H)+10e-10; 
        if t>1
            pvgh(:,1,t)=1; % random play strategy
            pvgh(:,2,t)=zeros(V,1); pvgh(o(t-1),2,t)=1; % do what you did last time
            pvgh(:,3,t)=zeros(V,1); pvgh(c(t-1),3,t)=1; % do what computer did last time
            pvgh(:,4,t)=zeros(V,1); pvgh(setdiff(1:3,union(c(t-1),o(t-1))),4,t)=1; % do different to human or computer
        end
    end
    pvgh=condp(pvgh+0.00001);
    [alpha,loglik]=HMMforwardTimeDependent(o(1:currenttime),phghm(:,:,1:currenttime),ph1,pvgh(:,:,1:currenttime));
    
    alpha(:,1)=condp(ones(H,1)); probnextmove(:,1)=condp(ones(V,1));
    % one step ahead prediction:
    t=currenttime+1;
    probnextmove(:,t)=pvgh(:,:,t)*phghm(:,:,t)*alpha(:,t-1);
    %human(t)=randgen(probnextmove(:,t)); % sample next move by the human
    human(t)=brml.argmax(probnextmove(:,t)); %  most likely next move by the human
    c(t)=human(t)+1; % computer plays a move to beat the human
    if c(t)==4; c(t)=1;end
   
    
    if showfiltered
        subplot(2,1,1)
        imagesc(alpha); title('filtered distribution');
        set(gca,'YTickLabel',{'Random','You Last Time','Computer last time', 'Different'}, 'YTick',1:4);
        set(gca,'XTick',1:t);
        subplot(2,1,2)
        imagesc(probnextmove); title('next move distribution');
        set(gca,'YTickLabel',{'Rock','Paper','Scissors'}, 'YTick',1:3);
        set(gca,'XTick',1:t);
    end
    
    
end