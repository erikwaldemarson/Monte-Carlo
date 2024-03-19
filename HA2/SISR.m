function c_n = SISR(d,n,N)
    %Same SISR as in 1.5

    %%Start of SISR
    X = zeros(n, d); %initial X, each row represents a coordinate in the random walk
    X(:,:, N) = 0; %each page represents a different particle
    w = ones(N, n); %weights, particle per row and step per column

    for k = 1:n
        for i = 1:N
            %selection
            z = 1:N;
            W = w(:,k)./sum(w(:,k));
            x_idx = randsample(z,1,true,W);

            %mutation
            X(:, :, i) = X(:, :, x_idx); 

            %self avoiding walk g
            isAllowed = ones(2,d); %0 if step is not allowed

            for j = 1:d
                for l = 1:k
                    step = zeros(1,d);
                    step(1,j) = 1;

                    if X(l, :, i) == X(k, :, i) + step
                        isAllowed(1,j) = 0;
                    end

                    if X(l, :, i) == X(k, :, i) - step 
                        isAllowed(2,j) = 0;
                    end
                end
            end

            nr_of_allowed = sum(isAllowed, 'all');
            set_of_allowed = cell(1, nr_of_allowed);

            %create set of all allowed steps
            x = 1;
            for l = 1:2
                for j = 1:length(isAllowed)
                    if isAllowed(l,j) == 1
                        set_of_allowed{x} = [l,j];
                        x = x + 1;
                    end
                end
            end

            if nr_of_allowed == 0 %if no steps are allowed then X_k+1 = X_k
                X(k+1, :, i) = X(k, :, i);
                w(i, k+1) = w(i, k);
            else
                step = zeros(1, d); %random step
                pos = randi([1,nr_of_allowed]); %choose random index

                idx = set_of_allowed{pos}(2);

                step(idx) = 1;

                if set_of_allowed{pos}(1) == 2 
                    step = -1*step;
                end

                X(k+1, :, i) = X(k, :, i) + step;
                w(i, k+1) = nr_of_allowed; 

            end

        end
    end
    %%end of SISR
    c_n = 1;
    for i = 2:n+1
        c_n = mean(w(:,i))*c_n; %c_n = average of weights
    end
end

