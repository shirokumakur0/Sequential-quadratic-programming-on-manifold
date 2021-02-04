        % qp_mineigval_pow = 1;
        % while true
            [coeff, ~, qpexitflag, ~, Lagmultipliers] = quadprog(qpinfo.H, qpinfo.f,...
                qpinfo.A, qpinfo.b, qpinfo.Aeq, qpinfo.beq, [], [], [], qpoptions);
                     
        %    if options.qp_loop <= 0
        %        fprintf("break, qp_loop <= 0 \n")
        %        break;
        %    end
        %    if qp_mineigval_pow >= options.qp_mineigval_max_pow 
        %        fprintf("break, pow max");
        %        break; 
        %    end
        %    if qpexitflag == -6 && (strcmp(options.modify_hessian, 'mineigval_matlab') || strcmp(options.modify_hessian, 'mineigval_manopt') )       
        %        fprintf("continue, qp_loop <= 0 \n")               
        %        options.mineigval_correction = options.mineigval_correction * 10;
        %        qpinfo.diagcoeff = max(options.mineigval_threshold,...
        %            abs(qpinfo.mineigval)) + options.mineigval_correction;
        %        qpinfo.H = H + qpinfo.diagcoeff * eye(qpinfo.n);
        %        qpinfo.H = 0.5 * (qpinfo.H.' + qpinfo.H);
        %        qp_mineigval_pow = qp_mineigval_pow + 1;
        %    else
        %        fprintf("break , anyway with %d\n", qp_mineigval_pow);
        %        break;
        %    end    
        % end
        % if qp_mineigval_pow ~= 1
        %    fprintf('afterwards%d, %d\n',qpexitflag,qp_mineigval_pow)
        % end
        
                % if strcmp(options.modify_hessian, "eye")
        %     df0 = meritproblem.M.inner(xCur, deltaXast, deltaXast);
        % elseif strcmp(options.modify_hessian, 'mineigval_matlab')
        % elseif strcmp(options.modify_hessian, 'mineigval_manopt')
        %     if qpinfo.mineigval < 0
        %         df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
        %             mus, lambdas) + qpinfo.mineigval_diagcoeff * deltaXast, deltaXast);
        %     else 
        %         df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
        %             mus, lambdas), deltaXast);
        %     end
        % else % Standard setting
        %     df0 = meritproblem.M.inner(xCur, hessLagrangian(xCur, deltaXast,...
        %             mus, lambdas), deltaXast);
        % end
        
        
        % qpsolnorm = stepsize * problem0.M.norm(xCur, deltaXast);
        
                %elseif qpsolnorm <= options.tolqpnorm
        %    fprintf('QP solution norm with stepsize tolerance reached\n');
        %    options.reason = "Legrangian Gradient norm tolerance reached";
        %    stop = true;
        
        % if strcmp(options.modify_hessian, 'mineigval_matlab') ||...
        %         strcmp(options.modify_hessian, 'mineigval_manopt')
        %    stats.qp_pow = NaN;
        % end

        