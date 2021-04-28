function demoRockPaperScissorsMarkovHuman(showprediction)
% play a game with a human opponent
close all
import brml.*


T=20;
[rock paper scissors]=assign(1:3);
s{1}='rock'; s{2}='paper'; s{3}='scissors';


for t=1:20 % number of games
    done=false;
    while ~done
        disp(' ')
        in=input(['Game ',num2str(t), ': Your Go! Rock(R), Paper (P) or Scissors (S)'],'s');
        switch lower(in)
            case 'r'
                o(t)=rock; done=true;
                
            case 'p'
                o(t)=paper; done=true;
            case 's'
                o(t)=scissors; done=true;
            otherwise
                disp('false input: must be r, p, or s. Try again.')
                done=false;
        end
    end
    
    if t<=2
        c(t)=randi(3);
    else
        % fit a Markov Model
        P=zeros(3,3);
        for tt=2:t-1
            P(o(tt),o(tt-1))=P(o(tt),o(tt-1))+1; % count the number of transitions
        end
        P=condp(P); % normalise the counts to get a probability
        prediction(:,t)=P(:,o(t-1));
        if showprediction
            imagesc(prediction); title('prediction')
            set(gca,'ytick',[1 2 3]);
            set(gca,'YTickLabel',{'Rock','Paper','Scissors'})
        end
        h(t)=randgen(P(:,o(t-1))); % sample next move by the human
        c(t)=h(t)+1; % computer plays a move to beat the human
        if c(t)==4; c(t)=1;end
    end
    
    disp(['you play ' s{o(t)},  ' , computer plays ' s{c(t)}]);
    
    w(t)=RPSwinner(o(t),c(t));
    switch w(t)
        case 1
            st='you win!';
        case 0
            st='draw!';
        case -1
            st='you lose!';
    end
    disp([st ' You won ' num2str(100*mean(w>0)) ' percent of the games up to now. Computer won '  num2str(100*mean(w<0)) ' percent'])
    
    
end
