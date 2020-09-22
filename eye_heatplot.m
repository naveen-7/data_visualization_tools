
function FULLMAT1 = makeeyegram_n(eyex,eyey)

res =2;
clear GRIDX GRIDY FULLMAT
[GRIDX,GRIDY] = meshgrid(-5:1/(10^res):5,-5:1/(10^res):5);
GRIDX = round(GRIDX,res); GRIDY = round(GRIDY,res);
FULLMAT = zeros(size(GRIDX));



eyex = EYE_POP_S{52,1}.EYE_Xn;
eyey = EYE_POP_S{52,1}.EYE_Yn;

for i=1:length(eyex)
    
    if ~isempty(eyex(i,1))
        
        i
        
        TEMPX = eyex{i,1}; TEMPY = eyey{i,1};
        TEMPX  = round(TEMPX,res); TEMPY  = round(TEMPY,res);
        if ~(any(TEMPX>5) | any(TEMPX<-5) | any(TEMPY>5) | any(TEMPY<-5))
            
            %        figure();
            %        subplot(2,2,1)
            %        hold on;
            %        plot(eyex{i,1},eyey{i,1});
            %        plot(TEMPX,TEMPY);
            %        xlim([-5 5]); ylim([-5 5]);
            %
            for ii=1:length(TEMPX)
                [~,COL] =  find(GRIDX==TEMPX(ii)); COL = COL(1);
                [ROW,~] =  find(GRIDY==TEMPY(ii)); ROW = ROW(1);
                
                FULLMAT(ROW,COL) = 1+FULLMAT(ROW,COL);
            end
        end
    end
end



FULLMAT = flipud(FULLMAT);
% figure();imagesc(FULLMAT); 

FULLMAT1 = mat2gray(FULLMAT);
% figure();imagesc(FULLMAT1); colormap('hot');

end


