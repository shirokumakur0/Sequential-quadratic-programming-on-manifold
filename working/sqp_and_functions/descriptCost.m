function descriptCost(meritproblem, xCur, deltaXast)
iter_num = 100;
stepwidth = 1e-8;

data.step = zeros(iter_num+10, 1);
data.cost = zeros(iter_num+10, 1);

step = 0;

data.step(1) = step;
data.cost(1) = meritproblem.cost(xCur);

for i = 1:iter_num
    step = step + stepwidth;
    newx = meritproblem.M.retr(xCur, deltaXast, step);
    data.step(i+1) = step;
    data.cost(i+1) = meritproblem.cost(newx);
end

figure;
semilogy([data.step], [data.cost], '.-');
xlabel('Step');
ylabel('Cost of the merit');
fprintf("cost initial:%f     cost next:%f \n",data.cost(1), data.cost(2));
end