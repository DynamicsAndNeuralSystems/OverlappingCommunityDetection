function[F1] = F1_calc(actual,predicted)
%Calculate binary F1 score 

idx=(actual()==1);
pos=length(actual(idx));
neg=length(actual)-pos;

tp=sum(actual(idx)==predicted(idx));
tn=sum(actual(~idx)==predicted(~idx));
fp=neg-tn;
fn=pos-tp;

prec=tp/(tp+fp);
rec=tp/pos;

if tp == 0
% If there are no true positives, F1 score is undefined. Avoiding that by
% replacing with a 0 instead. 
    F1 = 0;
else
    F1=2*((prec*rec)/(prec+rec));
end
