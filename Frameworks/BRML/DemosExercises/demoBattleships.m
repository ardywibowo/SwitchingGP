function occ=demoBattleships(standardplay,occ)
% occ is the occupancy state of the board
% To play a first game with a human opponent, use
% occ=demoBattleships(true,[])
% To then play the same game with a computer opponent, use
% demoBattleships(false,occ)

if isempty(occ)
    generateship=true;
else
    generateship=false;
end

import brml.*
D=10; S=D*D; M=reshape(1:S,D,D);

ship(1).x=[0 1 2 3 4]; ship(1).y=[0 0 0 0 0];
ship(2).x=[0 0 0 0 0]; ship(2).y=[0 1 2 3 4];

if generateship
    shiporigin(1).x=1; shiporigin(1).y=8;
    shiporigin(2).x=1; shiporigin(2).y=2;
    
    notok=true;
    while notok
        shiporigin(1).x=randi(10); shiporigin(1).y=randi(10);
        shiporigin(2).x=randi(10); shiporigin(2).y=randi(10);
        
        occ=false(D,D);
        notok=false;
        for s=1:2
            for i=1:length(ship(s).x)
                if validgridposition(shiporigin(s).x+ship(s).x(i),shiporigin(s).y+ship(s).y(i),D,D)
                    if  occ(shiporigin(s).x+ship(s).x(i),shiporigin(s).y+ship(s).y(i))
                        notok=true;
                    else
                        occ(shiporigin(s).x+ship(s).x(i),shiporigin(s).y+ship(s).y(i))=true;
                    end
                else
                    notok=true;
                end
            end
        end
    end
    
    %figure(2); imagesc(occ); pause
end


data=zeros(D,D);

for attempt=1:25
    
    figure(1);
    
    
    count=0; prior=zeros(S^2,1); post=zeros(S^2,1);
    poss=zeros(D,D,2); prob=zeros(D,D); postocc=zeros(D,D);
    for s1=1:S
        mask1=false(D,D);
        [x1,y1]=find(M==s1);
        onboard1=true;
        for i=1:length(ship(1).x)
            if ~validgridposition(x1+ship(1).x(i),y1+ship(1).y(i),D,D)
                onboard1=false;
            else
                mask1(x1+ship(1).x(i),y1+ship(1).y(i))=true;
            end
        end
        
        onboard2=true;
        if onboard1
            for s2=1:S
                onboard2=true;
                mask2=false(D,D);
                [x2,y2]=find(M==s2);
                count=count+1;
                
                for i=1:length(ship(2).x)
                    if ~validgridposition(x2+ship(2).x(i),y2+ship(2).y(i),D,D)
                        onboard2=false;
                    else
                        mask2(x2+ship(2).x(i),y2+ship(2).y(i))=true;
                    end
                end
                
                if onboard1 && onboard2
                    if ~any( mask1(:) & mask2(:)) % no clashing ships
                        prior(count)=1;
                        mask=mask1+mask2;
                        if ~( any(mask(:)==1 & data(:)==-1) || any(data(:)==1 & mask(:)==0)) % compatible with data
                            post(count)=1;
                            postocc=postocc+mask;
                        end
                    end
                end
            end
        end
    end
    
    postocc=(~standardplay)*postocc/sum(post(:));
    imagesc(postocc); colorbar; colormap bone
    for i=1:10
        line([0 11],[i-0.5 i-0.5],'color',[0 0 0])
        line([i-0.5 i-0.5],[0 11],'color',[0 0 0])
    end
    title(['attempt ' num2str(attempt)])
    
    [yy xx]=find(data==1);
    for i=1:length(xx)
        text(xx(i)-0.1,yy(i),'o','color',[0,1,0],'fontsize',20);
    end
    [yy xx]=find(data==-1);
    for i=1:length(xx)
        text(xx(i)-0.1,yy(i),'x','color',[1,0,0],'fontsize',20);
    end
    p=postocc.*(data==0);
    [g1,g2]=find(p==max(p(:)));
    
    if ~standardplay
        %fprintf(1,'posterior is highest at %d, %d\n',g2(1),g1(1));
        text(g2(1)-0.1,g1(1),'*','color',[0,0,1],'fontsize',20);
    end
    guess=fliplr(round(ginput(1)));
    if occ(guess(1),guess(2))
        data(guess(1),guess(2))=1;
    else
        data(guess(1),guess(2))=-1;
    end
    
end