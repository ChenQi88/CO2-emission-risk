tic

load = load1;
zerocost = 7504.1;%0cost generators including nuclear hydro and bio:ercot = 7504.1 pjm = 40667
netloadnw = load - zerocost;
netload = netloadnw - wexp;

h = size(load);
h = h(1);
n = size(mc);
n = n(1);

result = [];
p_value = 0.5;
corr = 0.99;

sum_co2w = [];
sum_generate = [];
sum_co2_per = [];
sum_wind = [];

wexp_original = wexp;
wsigma_original = wsigma;
netloadnw_original = netloadnw;

wind_ratio = [];
co2 = [];
generate = [];
co2_per_ge = [];
wind = [];
IVA = [];
a = [];

for w_ratio = 0.05 :0.05 :0.55

	w_ratio 
    
    sratio = w_ratio*sum(load(1:h,1))/sum(wexp_original);
	wexp = sratio*wexp_original;
	wsigma = sratio* wsigma_original;
 	solution = zeros(n,h);
	solution = zeros(n,h);
    solutionnw = zeros(n,h);
    
    result_co2 = [];
    result_generate = [];
    result_co2_per_ge = [];
    result_wind = [];
  
    wind_ratio = [wind_ratio; sum(wexp)/sum(load(1:h,1))];
    
    for random = 1:1000
        random
        RTload(1) = normrnd(wexp(1),wsigma(1));
        for i = 2:h
            wexp_new = wexp(i) + corr* wsigma(i)/wsigma(i-1)*(RTload(i-1)-wexp(i-1));
            wsigma_new = ((1-corr*corr)*wsigma(i)*wsigma(i))^0.5;
            RTload(i) = normrnd(wexp_new, wsigma_new);
        end

        RTload(find(RTload<0))=0;
        %adjust_coefficient = sum(wexp)/sum(RTload); %RPS scenario
        %RTload = adjust_coefficient*RTload;

        solutionRT = solution;
        solution_limit = 0.5*capacity;

        for i = 1:h
            resid = zeros(n,h); 
            resid(:,i) = min(solution_limit+ramp,capacity);
            if RTload(i) >= netloadnw_original(i)
                solution_limit = solutionRT(:,i);
            else
                r = netloadnw_original(i) - RTload(i);
                if r > sum(resid(:,i))
                    a = [a; r-sum(resid(:,i))];
                    r = sum(resid(:,i))-0.1;
                end
                j=1;
                while r - resid(j,i)> 0 %if the demand needed to be balanced is larger than the available capaicity of generator i in hour j in the RT market, the step width is the available capacity of generator i
                    solutionRT(j,i) = solutionRT(j,i) + resid(j,i);
                    r = r - resid(j,i);
                    resid(j,i) = 0;
                    j = j+1;
                end
                solutionRT(j,i) =  solutionRT(j,i) + r;
                solution_limit = solutionRT(:,i);
                resid(j,i) = 0;
                r = 0; 
            end
        end
        co2w = cer'*solutionRT;
        result_co2 = [result_co2 ; sum(co2w)];
        result_co2_per_ge = [result_co2_per_ge ;sum(co2w)/ sum(sum(solutionRT))];
        result_generate =  [result_generate ; sum(sum(solutionRT))];
        result_wind = [result_wind ; sum(RTload)];
    end
    sum_co2w = [sum_co2w, result_co2];
    sum_generate = [sum_generate, result_generate];
    sum_co2_per = [sum_co2_per, result_co2_per_ge];
    sum_wind = [sum_wind, result_wind];
    co2 = [co2; (prctile(sort(result_co2),(100-p_value))-prctile(sort(result_co2),(p_value)))/mean(result_co2)];
    generate = [generate; (prctile(sort(result_generate),(100-p_value))-prctile(sort(result_generate),(p_value)))/mean(result_generate)];
    co2_per_ge = [co2_per_ge; (prctile(sort(result_co2_per_ge),(100-p_value))-prctile(sort(result_co2_per_ge),(p_value)))/mean(result_co2_per_ge)];
    wind = [wind; (prctile(sort(result_wind),(100-p_value))-prctile(sort(result_wind),(p_value)))/mean(result_wind)];
    IVA = [IVA, std(result_wind)/mean(result_wind)];
end

result = [wind_ratio, co2, generate, co2_per_ge, wind  ];        
toc        


%Resampling 
N = 200;
e_co2 = [];
for i = 1:11
    new_co2 = sum_co2w(:,i);
    emission_co2 = [];
    for j = 1:200
        renew_co2 = randperm(numel(new_co2));
        a_co2 = new_co2(renew_co2(1:N));
        emission_co2 = [emission_co2; (prctile(sort(a_co2),(100-p_value))-prctile(sort(a_co2),(p_value)))/mean(a_co2)];
    end
    e_co2 = [e_co2, emission_co2];
end

        
        
        
