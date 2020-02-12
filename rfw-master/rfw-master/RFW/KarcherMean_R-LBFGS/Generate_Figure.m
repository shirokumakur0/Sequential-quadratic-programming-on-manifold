function Generate_Figure
ImportGROPT;

% specify which Figure do you want to generate
Fig = 2;
position = 'top';

data(Fig, position); % generate test set
testAH;
dataAveAH; % take average
Plot;
end